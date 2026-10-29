#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
#  07_equilibrate.sh  —  Equilibration: NVT (heating) -> NPT (density)
#                         GROMACS/GPU version  [PATCHED]
#
#  Fixes applied vs original:
#    [1] Hardcoded WORKDIR cd (config.sh had WORKDIR="")
#    [2] srun -n 1 launcher for gmx_mpi (was -ntmpi 1 — invalid for MPI build)
#    [3] Helper strips any pre-existing 'define' line in the MDP before writing
#        a new one. Fixes the bug where npt4 inherited '-DPOSRES' from a cp
#        of npt1.mdp, so Stage 5 was NOT actually free of restraints.
#    [4] -r reference structure is min/min2.gro (was the pre-minimization
#        $GMX_GRO, which would pull atoms toward un-minimized geometry).
#    [5] -t <prev>.cpt is passed to grompp between stages for clean
#        thermostat/barostat/box state continuation.
#
#  Stage 1  NVT  500 ps  |  V-rescale 300K  |  posres k=1000
#  Stage 2  NPT  500 ps  |  C-rescale 1 bar |  posres k=1000
#  Stage 3  NPT  500 ps  |                  |  posres k=100
#  Stage 4  NPT  500 ps  |                  |  posres k=10
#  Stage 5  NPT  500 ps  |                  |  no restraints  (FREE)
#
#  Total equilibration: 2.5 ns
#
#  Submit with: sbatch 07_equilibrate.sh
# ══════════════════════════════════════════════════════════════════════════════
#SBATCH --job-name=equil
#SBATCH --partition=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=equil_%j.log
#SBATCH --error=equil_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jccarrillo@uach.mx

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./config.sh
load_gromacs
module load python/3.11.9

# Log fijo adicional
exec > >(tee -a 07_equilibrate.log) 2>&1

echo ""
echo "══════════════════════════════════════════════════"
echo "  Step 07 — Equilibration (NVT -> NPT)"
echo "  Job ID : $SLURM_JOB_ID"
echo "  Node   : $SLURMD_NODENAME"
echo "══════════════════════════════════════════════════"

[[ ! -f "min/min2.gro" ]] && {
    echo "[ERROR] min/min2.gro not found — run Step 06 first."; exit 1; }

mkdir -p equil
NCPU="${SLURM_CPUS_PER_TASK:-10}"

# ──────────────────────────────────────────────────────────────────────────────
# Helper: run one equilibration stage
# Usage: run_stage STAGE_NAME INPUT_GRO MDP_FILE [POSRES_DEFINE]
# ──────────────────────────────────────────────────────────────────────────────
run_stage() {
    local STAGE="$1"
    local INPUT="$2"
    local MDP="$3"
    local DEFINE="${4:-}"

    echo ""
    echo "  --- $STAGE ---"

    # [3] Strip any pre-existing 'define' line so cp'd MDPs don't carry stale
    # POSRES directives. Then write the requested define (if any).
    sed -i '/^define[[:space:]]*=/d' "$MDP"
    if [[ -n "$DEFINE" ]]; then
        echo "define = $DEFINE" >> "$MDP"
    fi

    # [5] Pass previous stage's checkpoint for proper continuation.
    # Empty TFLAG for nvt (fresh start, gen_vel=yes).
    local TFLAG=""
    case "$STAGE" in
        npt1) [[ -f equil/nvt.cpt  ]] && TFLAG="-t equil/nvt.cpt"  ;;
        npt2) [[ -f equil/npt1.cpt ]] && TFLAG="-t equil/npt1.cpt" ;;
        npt3) [[ -f equil/npt2.cpt ]] && TFLAG="-t equil/npt2.cpt" ;;
        npt4) [[ -f equil/npt3.cpt ]] && TFLAG="-t equil/npt3.cpt" ;;
    esac

    # grompp -- preprocess (no MPI launcher needed, single process)
    # [4] -r points to min/min2.gro (the actual minimized state)
    $GMX grompp \
        -f "$MDP" \
        -c "$INPUT" \
        -r min/min2.gro \
        -p "$GMX_TOP" \
        -n "$GMX_NDX" \
        $TFLAG \
        -o "equil/${STAGE}.tpr" \
        -maxwarn 5 \
        2>&1 | tee "equil/${STAGE}_grompp.log"

    # [2] mdrun under srun -n 1 (gmx_mpi is an MPI build; -ntmpi is invalid)
    srun -n 1 $GMX mdrun \
        -v \
        -deffnm "equil/${STAGE}" \
        -s "equil/${STAGE}.tpr" \
        -ntomp "$NCPU" \
        2>&1 | tee "equil/${STAGE}_mdrun.log"

    echo "  ✓ $STAGE complete -> equil/${STAGE}.gro"
}

# ══════════════════════════════════════════════════════════════════════════════
#  Stage 1 — NVT: 500 ps, V-rescale 300K, strong restraints (k=1000)
# ══════════════════════════════════════════════════════════════════════════════
cat > equil/nvt.mdp << EOF
; NVT equilibration — 500 ps, ${TEMP} K, position restraints k=1000

integrator          = md
nsteps              = 250000        ; 500 ps at dt=0.002
dt                  = 0.002

; Output
nstxout-compressed  = 500           ; .xtc every 1 ps
nstenergy           = 500
nstlog              = 500

; Bonds
continuation        = no            ; fresh start (generate velocities)
constraint_algorithm= lincs
constraints         = h-bonds
lincs_iter          = 1
lincs_order         = 4

; Neighbor searching
cutoff-scheme       = Verlet
nstlist             = 10
ns_type             = grid

; Electrostatics
coulombtype         = PME
rcoulomb            = 1.0
pme_order           = 4
fourierspacing      = 0.16

; Van der Waals
rvdw                = 1.0

; Temperature — V-rescale (velocity rescaling with stochastic term)
tcoupl              = V-rescale
tc-grps             = System
tau_t               = 0.1
ref_t               = ${TEMP}

; Pressure — NVT (no pressure coupling)
pcoupl              = no

; Periodic boundary
pbc                 = xyz
DispCorr            = EnerPres

; Generate velocities at target temperature
gen_vel             = yes
gen_temp            = ${TEMP}
gen_seed            = -1
EOF

# Set posres k=1000 (already set by step 05, rerun here for safety)
python3 posres_setup.py 1000
run_stage "nvt" "min/min2.gro" "equil/nvt.mdp" "-DPOSRES"

# ══════════════════════════════════════════════════════════════════════════════
#  Stage 2 — NPT: 500 ps, C-rescale, full restraints (k=1000)
# ══════════════════════════════════════════════════════════════════════════════
cat > equil/npt1.mdp << EOF
; NPT stage 1 — 500 ps, density relaxation, posres k=1000

integrator          = md
nsteps              = 250000
dt                  = 0.002

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

; Temperature
tcoupl              = V-rescale
tc-grps             = System
tau_t               = 0.1
ref_t               = ${TEMP}

; Pressure — C-rescale (modern equilibration barostat)
pcoupl              = C-rescale
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = ${PRESSURE}
compressibility     = 4.5e-5

pbc                 = xyz
DispCorr            = EnerPres
gen_vel             = no
EOF

run_stage "npt1" "equil/nvt.gro" "equil/npt1.mdp" "-DPOSRES"

# ══════════════════════════════════════════════════════════════════════════════
#  Stage 3 — NPT: 500 ps, release restraints to k=100
# ══════════════════════════════════════════════════════════════════════════════
cp equil/npt1.mdp equil/npt2.mdp
python3 posres_setup.py 100
run_stage "npt2" "equil/npt1.gro" "equil/npt2.mdp" "-DPOSRES"

# ══════════════════════════════════════════════════════════════════════════════
#  Stage 4 — NPT: 500 ps, release restraints to k=10
# ══════════════════════════════════════════════════════════════════════════════
cp equil/npt1.mdp equil/npt3.mdp
python3 posres_setup.py 10
run_stage "npt3" "equil/npt2.gro" "equil/npt3.mdp" "-DPOSRES"

# ══════════════════════════════════════════════════════════════════════════════
#  Stage 5 — NPT: 500 ps, NO restraints (free equilibration)
#  Note: helper will strip the inherited 'define = -DPOSRES' from npt1.mdp.
# ══════════════════════════════════════════════════════════════════════════════
cp equil/npt1.mdp equil/npt4.mdp
run_stage "npt4" "equil/npt3.gro" "equil/npt4.mdp"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════════"
echo "  Step 07 complete."
echo "  Final structure : equil/npt4.gro"
echo "  Final checkpoint: equil/npt4.cpt"
echo "  Total equilibration: 2.5 ns"
echo ""
echo "  Check equilibration quality:"
echo "    echo Temperature | $GMX energy -f equil/npt4.edr -o equil/temperature.xvg"
echo "    echo Pressure    | $GMX energy -f equil/npt4.edr -o equil/pressure.xvg"
echo "    echo Density     | $GMX energy -f equil/npt4.edr -o equil/density.xvg"
echo "    - Temperature should be stable at ${TEMP} K"
echo "    - Pressure should fluctuate around ${PRESSURE} bar (large fluctuations are normal)"
echo "    - Density should converge (~1.0 g/cm^3 for OPC water)"
echo ""
echo "  NEXT: sbatch 08_production.sh"
echo "══════════════════════════════════════════════════"
