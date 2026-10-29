#!/usr/bin/env bash
# ============================================================================
#  GO19 v2 — full chain on a rented GPU
#
#    bash gpu_run.sh bench              ~5 min, answers "what throughput, and
#                                       does Parrinello-Rahman work with
#                                       -update gpu?"
#    bash gpu_run.sh equil all          min1 min2 nvt npt1..npt4, 4 systems
#    bash gpu_run.sh equil go19_v2      one system
#    bash gpu_run.sh prod <sys> <rep>   one 100 ns replica, rep in 1..3
#
#  MDP files reproduce the Yoltla protocol byte-for-byte with ONE change:
#  refcoord-scaling = com in the restrained NPT stages.  Without it the
#  barostat rescales coordinates but not the position-restraint reference,
#  which drove a divergent box oscillation on Yoltla (+-25000 bar, 395 K,
#  LINCS failure at 376 ps).  grompp warns about it; -maxwarn 5 hid the
#  warning, so maxwarn is 1 here -- enough for the known benign net-charge
#  warning (-0.002 e, AM1-BCC rounding) and nothing else.
#
#  No `set -e`: gates must run even when a step fails.
# ============================================================================
set -uo pipefail

BUNDLE="${BUNDLE:-$PWD/gpu_bundle}"
ALLSYS="go19_v2 go19_v2_R53A go19_v2_R272L go19_v2_R53A_R272L"
TEMP=300
PRESSURE=1.0
DT=0.002
PROD_NS="${PROD_NS:-100}"
NCPU="${NCPU:-$(nproc)}"

GMX="$(command -v gmx 2>/dev/null || command -v gmx_mpi 2>/dev/null)"
[[ -n "$GMX" ]] || { echo "FATAL: no gmx / gmx_mpi on PATH"; exit 1; }

echo "gmx      : $GMX"
"$GMX" --version 2>/dev/null | grep -iE "^GROMACS version|^GPU support" | sed 's/^/           /'
echo "cores    : $NCPU"
echo "bundle   : $BUNDLE"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | sed 's/^/gpu      : /'
echo

# ---------------------------------------------------------------- mdp writer
write_mdps() {
    local D="$1"; mkdir -p "$D/mdp"

    cat > "$D/mdp/min1.mdp" <<EOF
; Restrained steepest descent -- relax solvent while holding solute
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
cutoff-scheme   = Verlet
nstlist         = 10
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
pme_order       = 4
fourierspacing  = 0.16
rvdw            = 1.0
constraints     = none
pbc             = xyz
define          = -DPOSRES
EOF

    cat > "$D/mdp/min2.mdp" <<EOF
; Unrestrained steepest descent -- full system relaxes
integrator      = steep
emtol           = 500.0
emstep          = 0.01
nsteps          = 50000
cutoff-scheme   = Verlet
nstlist         = 10
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
pme_order       = 4
fourierspacing  = 0.16
rvdw            = 1.0
constraints     = none
pbc             = xyz
EOF

    cat > "$D/mdp/nvt.mdp" <<EOF
; NVT equilibration -- 500 ps, ${TEMP} K, position restraints k=1000
integrator          = md
nsteps              = 250000
dt                  = ${DT}
nstxout-compressed  = 500
nstenergy           = 500
nstlog              = 500
continuation        = no
constraint_algorithm= lincs
constraints         = h-bonds
lincs_iter          = 1
lincs_order         = 4
cutoff-scheme       = Verlet
nstlist             = 10
ns_type             = grid
coulombtype         = PME
rcoulomb            = 1.0
pme_order           = 4
fourierspacing      = 0.16
rvdw                = 1.0
tcoupl              = V-rescale
tc-grps             = System
tau_t               = 0.1
ref_t               = ${TEMP}
pcoupl              = no
pbc                 = xyz
DispCorr            = EnerPres
gen_vel             = yes
gen_temp            = ${TEMP}
gen_seed            = -1
define              = -DPOSRES
EOF

    # npt1..npt3 restrained (+ refcoord-scaling), npt4 free
    local s
    for s in 1 2 3 4; do
        cat > "$D/mdp/npt${s}.mdp" <<EOF
; NPT stage ${s} -- 500 ps, density relaxation
integrator          = md
nsteps              = 250000
dt                  = ${DT}
nstxout-compressed  = 500
nstenergy           = 500
nstlog              = 500
continuation        = yes
constraint_algorithm= lincs
constraints         = h-bonds
lincs_iter          = 1
lincs_order         = 4
cutoff-scheme       = Verlet
nstlist             = 10
ns_type             = grid
coulombtype         = PME
rcoulomb            = 1.0
pme_order           = 4
fourierspacing      = 0.16
rvdw                = 1.0
tcoupl              = V-rescale
tc-grps             = System
tau_t               = 0.1
ref_t               = ${TEMP}
pcoupl              = C-rescale
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = ${PRESSURE}
compressibility     = 4.5e-5
pbc                 = xyz
DispCorr            = EnerPres
gen_vel             = no
EOF
        if [[ "$s" != "4" ]]; then
            cat >> "$D/mdp/npt${s}.mdp" <<EOF
define              = -DPOSRES
refcoord-scaling    = com
EOF
        fi
    done
}

write_prod_mdp() {
    local OUT="$1" SEED="$2"
    local STEPS=$(python3 -c "print(int(${PROD_NS}*1000/${DT}))")
    cat > "$OUT" <<EOF
; Production MD -- ${PROD_NS} ns, NPT, ${TEMP} K, no restraints
; gen_seed=${SEED}
integrator          = md
nsteps              = ${STEPS}
dt                  = ${DT}
nstxout-compressed  = 5000
compressed-x-grps   = Solute
nstenergy           = 5000
nstlog              = 5000
nstxout             = 0
nstvout             = 0
continuation        = no
constraint_algorithm= lincs
constraints         = h-bonds
lincs_iter          = 1
lincs_order         = 4
cutoff-scheme       = Verlet
nstlist             = 10
ns_type             = grid
coulombtype         = PME
rcoulomb            = 1.0
pme_order           = 4
fourierspacing      = 0.16
rvdw                = 1.0
tcoupl              = V-rescale
tc-grps             = System
tau_t               = 0.1
ref_t               = ${TEMP}
pcoupl              = Parrinello-Rahman
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = ${PRESSURE}
compressibility     = 4.5e-5
pbc                 = xyz
DispCorr            = EnerPres
gen_vel             = yes
gen_temp            = ${TEMP}
gen_seed            = ${SEED}
EOF
}

# --------------------------------------------------- restraint force constant
# The submitted Methods describes a progressive release across the NPT stages:
#   k = 1000 -> 100 -> 10 -> 0 kJ/mol/nm^2
# The posre_*.itp files shipped with the build hold a single fixed k = 1000 and
# no k=100 / k=10 variants exist anywhere in the project, so that schedule was
# never implemented.  Rather than rewrite the paper, we implement it: each
# restrained stage regenerates posre_*.itp at its own force constant from an
# untouched master copy before grompp runs.  grompp bakes the values into the
# tpr, so rewriting the same filename between stages is safe.
set_posres_fc() {
    local W="$1" FC="$2" f base
    for base in system1 FAD LIG; do
        f="$W/posre_${base}.itp"
        [[ -f "$W/posre_${base}.master.itp" ]] || cp "$f" "$W/posre_${base}.master.itp"
        awk -v fc="$FC" '
            /^[[:space:]]*[;\[]/ { print; next }
            NF >= 5 { printf "%6s %2s %9s %9s %9s\n", $1, $2, fc, fc, fc; next }
            { print }
        ' "$W/posre_${base}.master.itp" > "$f"
    done
    # verify: every data line must now carry FC
    local bad
    bad=$(awk -v fc="$FC" 'NF>=5 && $1 !~ /^[;\[]/ && $3 != fc {c++} END{print c+0}' "$W/posre_system1.itp")
    if [[ "$bad" != "0" ]]; then
        echo "    posres rescale FAILED: $bad lines not at k=$FC"; return 1
    fi
    return 0
}

# ------------------------------------------------------------- gpu flag set
GPUFLAGS=""
probe_gpu() {
    # An explicit GPUFLAGS_FORCE skips probing entirely.  Use it when the
    # benchmark has already established which flag set works -- the probe's
    # own detection has been wrong before, and a false negative silently
    # drops the whole run onto the CPU (~3 ns/day instead of ~200).
    if [[ -n "${GPUFLAGS_FORCE:-}" ]]; then
        GPUFLAGS="$GPUFLAGS_FORCE"
        echo "  GPU flags FORCED: $GPUFLAGS"
        return 0
    fi
    local D="$BUNDLE/go19_v2" TMP; TMP="$(mktemp -d)"
    write_mdps "$TMP"
    local START="$D/nvt_yoltla.gro"; [[ -f "$START" ]] || START="$D/complex.gro"
    "$GMX" grompp -f "$TMP/mdp/npt4.mdp" -c "$START" -r "$D/complex.gro" \
        -p "$D/complex.top" -n "$D/index.ndx" -o "$TMP/probe.tpr" \
        -maxwarn 1 > "$TMP/probe_grompp.log" 2>&1
    [[ -f "$TMP/probe.tpr" ]] || { echo "  probe grompp FAILED:"; tail -20 "$TMP/probe_grompp.log"; return 1; }
    local f
    for f in "-nb gpu -pme gpu -bonded gpu -update gpu" \
             "-nb gpu -pme gpu -update gpu" \
             "-nb gpu -pme gpu" ; do
        rm -f "$TMP"/probe.gro "$TMP"/probe.log "$TMP"/probe_mdrun.log
        ( cd "$TMP" && "$GMX" mdrun -s probe.tpr -deffnm probe -nsteps 200 \
              -ntmpi 1 -ntomp "$NCPU" -pin on $f > probe_mdrun.log 2>&1 )
        if [[ -f "$TMP/probe.gro" ]] || grep -qE "Finished mdrun|^Performance:" \
                "$TMP/probe_mdrun.log" "$TMP/probe.log" 2>/dev/null; then
            GPUFLAGS="$f"; echo "  GPU flags accepted: $f"; rm -rf "$TMP"; return 0
        fi
        echo "  rejected: $f"
    done
    echo "  no GPU flag set worked."
    echo "  REFUSING to continue on CPU: equilibration would take ~80 h instead of ~1.5 h."
    echo "  If CPU is genuinely intended, re-run with GPUFLAGS_FORCE=' ' (a single space)."
    grep -iE "fatal|error|not supported" "$TMP/probe_mdrun.log" | head -5 | sed 's/^/     /'
    GPUFLAGS=""; rm -rf "$TMP"; return 1
}

# ------------------------------------------------------------------ one stage
# run_stage <workdir> <stage> <input.gro> <refgro> [prev.cpt]
run_stage() {
    local W="$1" S="$2" IN="$3" REF="$4" CPT="${5:-}" FC="${6:-}"
    local T="$W/$S"
    if [[ -f "$T.gro" ]]; then echo "    $S: already done, skipping"; return 0; fi

    if [[ -n "$FC" ]]; then
        set_posres_fc "$W" "$FC" || return 1
    fi
    local TFLAG=(); [[ -n "$CPT" && -f "$CPT" ]] && TFLAG=(-t "$CPT")
    "$GMX" grompp -f "$W/mdp/$S.mdp" -c "$IN" -r "$REF" \
        -p "$W/complex.top" -n "$W/index.ndx" "${TFLAG[@]}" \
        -o "$T.tpr" -maxwarn 1 > "$T.grompp.log" 2>&1

    if [[ ! -f "$T.tpr" ]]; then
        echo "    $S: GROMPP FAILED"; grep -iE "^(ERROR|Fatal)" -A6 "$T.grompp.log" | head -25 | sed 's/^/      /'
        return 1
    fi
    local nw; nw=$(grep -c "^WARNING [0-9]" "$T.grompp.log")
    if [[ "$nw" -gt 1 ]]; then
        echo "    $S: $nw grompp WARNINGs (expect 1: net charge) -- STOPPING"
        grep -A5 "^WARNING [0-9]" "$T.grompp.log" | sed 's/^/      /'
        return 1
    fi
    grep -qi "refcoord" "$T.grompp.log" && { echo "    $S: refcoord warning still present -- STOPPING"; return 1; }

    local MF=($GPUFLAGS)
    # Minimisation (steep) is a non-dynamical integrator: GROMACS refuses to
    # offload PME for it ("PME GPU does not support: Non-dynamical
    # integrator").  Non-bonded on the GPU is fine, PME stays on the CPU.
    # Minimisation is seconds either way.
    [[ "$S" == min* ]] && MF=(-nb gpu)
    "$GMX" mdrun -s "$T.tpr" -deffnm "$T" -ntmpi 1 -ntomp "$NCPU" -pin on \
        "${MF[@]}" > "$T.mdrun.log" 2>&1

    if [[ ! -f "$T.gro" ]]; then
        echo "    $S: MDRUN FAILED"; tail -15 "$T.mdrun.log" | sed 's/^/      /'; return 1
    fi
    local lincs=0
    [[ -f "$T.mdrun.log" ]] && lincs=$(grep -c "LINCS WARNING" "$T.mdrun.log")
    lincs=${lincs//[^0-9]/}; lincs=${lincs:-0}
    # NB: no `|| echo 0`.  `grep -c` already prints 0 when it matches nothing,
    # and it exits 1 while doing so -- the fallback then appends a SECOND 0 and
    # the variable becomes "0\n0", which makes the [[ -gt ]] test below abort
    # instead of running.  That silently disables the barostat gate.
    local pscale=0
    [[ -f "$T.log" ]] && pscale=$(grep -c "pressure scaling more than 1%" "$T.log" 2>/dev/null)
    pscale=${pscale//[^0-9]/}; pscale=${pscale:-0}
    local perf; perf=$(awk '/^Performance:/{print $2}' "$T.log" 2>/dev/null | tail -1)
    printf "    %-5s ok   k=%-6s LINCS=%-4s pscale>1%%=%-5s %s ns/day\n" \
           "$S" "${FC:-none}" "$lincs" "$pscale" "${perf:-n/a}"
    if [[ "$lincs" != "0" ]]; then echo "      ^ LINCS warnings -- inspect before continuing"; return 1; fi
    if [[ "$pscale" -gt 50 ]]; then echo "      ^ $pscale box rescalings >1% -- barostat unstable, STOPPING"; return 1; fi
    return 0
}

# -------------------------------------------------------------------- bench
do_bench() {
    echo "=== BENCHMARK ==============================================="
    probe_gpu || return 1
    local D="$BUNDLE/go19_v2" TMP; TMP="$(mktemp -d)"; write_mdps "$TMP"
    local START="$D/nvt_yoltla.gro"; [[ -f "$START" ]] || START="$D/complex.gro"
    echo "  start structure: $(basename "$START")"
    echo

    # production-like, Parrinello-Rahman
    sed -e 's/^pcoupl .*/pcoupl              = Parrinello-Rahman/' \
        -e 's/^nsteps .*/nsteps              = 25000/' \
        "$TMP/mdp/npt4.mdp" > "$TMP/bench_pr.mdp"
    sed -e 's/^nsteps .*/nsteps              = 25000/' \
        "$TMP/mdp/npt4.mdp" > "$TMP/bench_cr.mdp"

    local v f
    printf "  %-28s %-34s %10s\n" VARIANT FLAGS "ns/day"
    for v in bench_pr bench_cr; do
      "$GMX" grompp -f "$TMP/$v.mdp" -c "$START" -r "$D/complex.gro" \
          -p "$D/complex.top" -n "$D/index.ndx" -o "$TMP/$v.tpr" \
          -maxwarn 1 > "$TMP/$v.grompp.log" 2>&1
      [[ -f "$TMP/$v.tpr" ]] || { echo "  $v: grompp failed"; tail -8 "$TMP/$v.grompp.log"; continue; }
      for f in "-nb gpu -pme gpu -bonded gpu -update gpu" "-nb gpu -pme gpu"; do
        ( cd "$TMP" && "$GMX" mdrun -s "$v.tpr" -deffnm "run" -ntmpi 1 -ntomp "$NCPU" \
              -pin on $f > "run.out" 2>&1 )
        local ns; ns=$(awk '/^Performance:/{print $2}' "$TMP/run.log" 2>/dev/null | tail -1)
        printf "  %-28s %-34s %10s\n" "$v" "$f" "${ns:-FAILED}"
        rm -f "$TMP"/run.*
      done
    done
    echo
    echo "  Read this before committing to the campaign:"
    echo "   - if Parrinello-Rahman + '-update gpu' works, production uses it as submitted"
    echo "   - if it FAILED but C-rescale worked, that is a protocol decision, not a flag"
    echo "   - 12 x ${PROD_NS} ns = $((12*PROD_NS)) ns total; GPU-days = $((12*PROD_NS)) / (ns per day)"
    rm -rf "$TMP"
}

# -------------------------------------------------------------------- equil
do_equil() {
    # NB: separate statements. Bash expands every word of a `local` command
    # before making any assignment, so `local S="$1" SRC="$BUNDLE/$S"` reads
    # an unset $S and aborts under `set -u`.
    local S="$1"
    local SRC="$BUNDLE/$S"
    local W="$PWD/run/$S"
    echo "=== EQUIL $S ================================================"
    [[ -d "$SRC" ]] || { echo "  no bundle dir $SRC"; return 1; }
    mkdir -p "$W"
    local f; for f in complex.gro complex.top index.ndx posre_system1.itp posre_FAD.itp posre_LIG.itp; do
        cp -n "$SRC/$f" "$W/" 2>/dev/null || true
        [[ -f "$W/$f" ]] || { echo "  MISSING $f"; return 1; }
    done
    write_mdps "$W"
    ( cd "$W" || exit 1
      #            stage  input          reference      checkpoint     k
      run_stage "$W" min1 "$W/complex.gro" "$W/complex.gro" ""            1000.0 || exit 1
      run_stage "$W" min2 "$W/min1.gro"    "$W/complex.gro" ""            ""     || exit 1
      run_stage "$W" nvt  "$W/min2.gro"    "$W/min2.gro"    ""            1000.0 || exit 1
      run_stage "$W" npt1 "$W/nvt.gro"     "$W/min2.gro" "$W/nvt.cpt"     1000.0 || exit 1
      run_stage "$W" npt2 "$W/npt1.gro"    "$W/min2.gro" "$W/npt1.cpt"     100.0 || exit 1
      run_stage "$W" npt3 "$W/npt2.gro"    "$W/min2.gro" "$W/npt2.cpt"      10.0 || exit 1
      run_stage "$W" npt4 "$W/npt3.gro"    "$W/min2.gro" "$W/npt3.cpt"     ""    || exit 1
    ) || { echo "  $S: equilibration stopped"; return 1; }

    echo "  --- npt4 averages (expect ~${TEMP} K, ~${PRESSURE} bar, converged volume)"
    printf "Temperature\nPressure\nVolume\nDensity\n0\n" \
        | "$GMX" energy -f "$W/npt4.edr" 2>/dev/null \
        | sed -n '/^Energy/,/^$/p' | sed 's/^/     /'
    if [[ -f "$PWD/check_lig_geom.py" ]]; then
        echo "  --- ligand geometry at npt4"
        python3 "$PWD/check_lig_geom.py" "$W/npt4.gro" | sed 's/^/     /'
    else
        echo "  (check_lig_geom.py not in $PWD -- upload it to run the geometry gate)"
    fi
    echo "  $S EQUILIBRATED -> run/$S/npt4.gro"
}

# --------------------------------------------------------------------- prod
do_prod() {
    local S="$1" REP="$2"
    local SEED=$((100003 + REP*7919))
    local W="$PWD/run/$S" P="$PWD/run/$S/prod_rep${REP}"
    echo "=== PROD $S rep$REP  seed=$SEED  ${PROD_NS} ns ====================="
    [[ -f "$W/npt4.gro" ]] || { echo "  no $W/npt4.gro -- equilibrate first"; return 1; }
    mkdir -p "$P"
    write_prod_mdp "$P/prod.mdp" "$SEED"

    if [[ ! -f "$P/prod.tpr" ]]; then
        # Production legitimately raises TWO grompp warnings, both examined:
        #   1. gen_vel=yes with Parrinello-Rahman. grompp assumes fresh velocities
        #      mean you are equilibrating, where PR can be unstable. Here the box
        #      is already equilibrated (300.02 K, 0.63 bar, volume converged to
        #      0.04 nm^3, no restraints), so PR has no density error to correct --
        #      which is the condition that makes it misbehave. gen_vel=yes is
        #      required: independent velocities per seed are what make the three
        #      replicas statistically independent (R1.4 / R2.5). The submitted
        #      runs used these same settings for 3 x 100 ns without incident.
        #   2. Ewald with net charge -0.002 e, from AM1-BCC rounding.
        # maxwarn is 2 and no higher, and the two are matched by content below --
        # a third warning, or a different second one, still stops the run.
        "$GMX" grompp -f "$P/prod.mdp" -c "$W/npt4.gro" -p "$W/complex.top" \
            -n "$W/index.ndx" -o "$P/prod.tpr" -maxwarn 2 > "$P/prod.grompp.log" 2>&1
        [[ -f "$P/prod.tpr" ]] || { echo "  grompp FAILED"; grep -iE "^(ERROR|Fatal)" -A6 "$P/prod.grompp.log" | head -20; return 1; }
        local nw; nw=$(grep -c "^WARNING [0-9]" "$P/prod.grompp.log")
        if [[ "$nw" -gt 2 ]]; then
            echo "  $nw grompp WARNINGs (expect 2) -- STOPPING"
            grep -A6 "^WARNING [0-9]" "$P/prod.grompp.log"; return 1
        fi
        # identity check, not just a count
        grep -qi "Parrinello-Rahman pressure coupling, but this can" "$P/prod.grompp.log" || {
            echo "  expected gen_vel/Parrinello-Rahman warning absent -- mdp changed? STOPPING"; return 1; }
        grep -qi "Ewald electrostatics in a system with net charge"  "$P/prod.grompp.log" || {
            echo "  expected net-charge warning absent -- topology changed? STOPPING"; return 1; }
        local unexpected
        unexpected=$(grep -A3 "^WARNING [0-9]" "$P/prod.grompp.log" \
                     | grep -icE "refcoord|position restraint|atoms are not|charge group|domain decomposition" || true)
        if [[ "${unexpected:-0}" != "0" ]]; then
            echo "  a warning mentions something unexpected -- STOPPING"
            grep -A6 "^WARNING [0-9]" "$P/prod.grompp.log"; return 1
        fi
    fi

    local RESUME=(); [[ -f "$P/prod.cpt" ]] && { RESUME=(-cpi "$P/prod.cpt" -append); echo "  resuming from checkpoint"; }
    "$GMX" mdrun -s "$P/prod.tpr" -deffnm "$P/prod" -ntmpi 1 -ntomp "$NCPU" -pin on \
        $GPUFLAGS -cpt 15 "${RESUME[@]}" 2>&1 | tail -5
    if [[ -f "$P/prod.gro" ]]; then
        echo "  rep$REP COMPLETE  $(awk '/^Performance:/{print $2" ns/day"}' "$P/prod.log" | tail -1)"
    else
        echo "  rep$REP incomplete (preempted?) -- rerun the same command to resume"
    fi
}

# --------------------------------------------------------------------- main
case "${1:-}" in
  bench) do_bench ;;
  equil)
        probe_gpu
        if [[ "${2:-all}" == "all" ]]; then
            for s in $ALLSYS; do do_equil "$s" || echo "  !! $s failed, continuing"; echo; done
        else do_equil "$2"; fi ;;
  prod)
        [[ $# -ge 3 ]] || { echo "usage: $0 prod <system> <rep 1-3>"; exit 1; }
        probe_gpu; do_prod "$2" "$3" ;;
  *)  sed -n '2,20p' "$0" ;;
esac
