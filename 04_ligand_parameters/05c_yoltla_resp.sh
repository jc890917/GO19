#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════
#  05c_yoltla_resp.sh  — replaces 05b
#
#  Root cause of the 05/05b failure:
#      antechamber: error while loading shared libraries: libimf.so
#  libimf.so is Intel's math library. AmberTools19 on Yoltla was built with
#  the Intel compilers, so its binaries need that runtime on LD_LIBRARY_PATH.
#  I had exported PATH by hand instead of loading the modules. Your own
#  config.sh load_amber() gets this right, and this script now uses the same
#  module set:
#      amber/22, amber/ambertools19, intel-oneapi/compilers-2022.2.1,
#      anaconda/2021.04, gromacs/2024.3-gpu
#
#  It also now VERIFIES each binary actually executes, instead of only
#  checking that the file exists and is marked executable — which is why the
#  earlier probe wrongly reported these tools as working.
#
#  Inputs in the current directory: LIG_all.esp  LIG_optA.sdf  nconf.txt
#  Run:  bash 05c_yoltla_resp.sh 2>&1 | tee cp3_report.txt
# ══════════════════════════════════════════════════════════════════════════
set -uo pipefail
HR="──────────────────────────────────────────────────────────────────"
die () { echo; echo "[STOP] $*"; echo "Send me cp3_report.txt and the *.log files."; exit 1; }

# ── modules, exactly as the working pipeline loads them ──────────────────
echo "$HR"; echo " 0. modules"; echo "$HR"
if ! type module >/dev/null 2>&1; then
    for INIT in /etc/profile.d/modules.sh /usr/share/Modules/init/bash \
                /opt/ohpc/admin/lmod/lmod/init/bash /etc/profile.d/lmod.sh; do
        [[ -r "$INIT" ]] && { . "$INIT"; break; }
    done
fi
type module >/dev/null 2>&1 || die "the 'module' command is not available in this shell"
module purge 2>/dev/null
for M in amber/22 amber/ambertools19 intel-oneapi/compilers-2022.2.1 \
         anaconda/2021.04 gromacs/2024.3-gpu; do
    if module load "$M" 2>/dev/null; then echo "  loaded  $M"
    else echo "  [!] could not load $M (continuing)"; fi
done

# ── verify the tools RUN, not merely exist ───────────────────────────────
echo; echo "$HR"; echo " 1. verifying AmberTools actually executes"; echo "$HR"
runs () {   # runs <binary> -> 0 if it executes without a loader error
    local out; out=$("$1" -h 2>&1 </dev/null); local rc=$?
    if grep -qi "error while loading shared libraries\|cannot open shared object" <<<"$out"; then
        echo "$out" | head -2 | sed 's/^/      /'; return 2
    fi
    return 0
}
for b in antechamber respgen espgen resp parmchk2; do
    P=$(command -v "$b" 2>/dev/null) || die "$b not on PATH after module load"
    if runs "$b"; then printf "  ✓ %-12s %s\n" "$b" "$P"
    else
        echo "  ✗ $b still cannot load its libraries."
        echo "    LD_LIBRARY_PATH is:"; tr ':' '\n' <<<"${LD_LIBRARY_PATH:-}" | sed 's/^/      /' | head -12
        echo "    Try in an interactive shell:  module load intel-oneapi/compilers-2022.2.1"
        echo "    then:  ldd \$(command -v $b) | grep -i 'not found'"
        die "$b is unusable"
    fi
done

# ── inputs ───────────────────────────────────────────────────────────────
echo; echo "$HR"; echo " 2. inputs"; echo "$HR"
for f in LIG_all.esp LIG_optA.sdf nconf.txt; do
    [[ -s "$f" ]] || die "missing or empty: $f"
done
NCONF=$(tr -dc '0-9' < nconf.txt)
NATM=$(awk 'NR==1{print $1+0; exit}' LIG_all.esp)
NPTS=$(awk 'NR==1{print $2+0; exit}' LIG_all.esp)
echo "  LIG_all.esp : $NCONF conformer(s), $NATM atoms, $NPTS ESP points"
[[ "$NATM" == "16" ]] || die "LIG_all.esp header says $NATM atoms, expected 16"

# ── 3. atom types (no charge method; respgen needs only types + bonds) ────
echo; echo "$HR"; echo " 3. GAFF2 atom types"; echo "$HR"
BUILT=""
for FMT in mdl sdf; do
    rm -f LIG.ac
    echo "  trying -fi $FMT ..."
    antechamber -i LIG_optA.sdf -fi "$FMT" -o LIG.ac -fo ac \
                -at gaff2 -nc -2 -rn LIG > "ac_$FMT.log" 2>&1
    if [[ -s LIG.ac ]]; then echo "  ✓ LIG.ac written with -fi $FMT"; BUILT="$FMT"; break; fi
    echo "    failed. Tail:"; tail -12 "ac_$FMT.log" | sed 's/^/      /'
done
if [[ -z "$BUILT" ]]; then
    echo; echo "  Formats this antechamber accepts:"; antechamber -L 2>&1 | head -30 | sed 's/^/    /'
    die "could not build LIG.ac"
fi
echo; echo "  atom types:"
awk '/^ATOM/{printf "    %2d  %-5s %-5s\n", $2, $3, $NF}' LIG.ac
NAC=$(grep -c '^ATOM' LIG.ac)
[[ "$NAC" == "16" ]] || die "LIG.ac has $NAC atoms, expected 16"

# ── 4. two-stage RESP, AmberTools reference implementation ───────────────
echo; echo "$HR"; echo " 4. two-stage RESP"; echo "$HR"
respgen -i LIG.ac -o LIG.respin1 -f resp1 -n "$NCONF" > respgen1.log 2>&1
[[ -s LIG.respin1 ]] || { tail -15 respgen1.log | sed 's/^/    /'; die "respgen resp1 failed"; }
respgen -i LIG.ac -o LIG.respin2 -f resp2 -n "$NCONF" > respgen2.log 2>&1
[[ -s LIG.respin2 ]] || { tail -15 respgen2.log | sed 's/^/    /'; die "respgen resp2 failed"; }
echo "  ✓ respin1 / respin2"

resp -O -i LIG.respin1 -o LIG.respout1 -p LIG.esout1 -t LIG.qout1 -e LIG_all.esp > resp1.log 2>&1
[[ -s LIG.qout1 ]] || { tail -25 LIG.respout1 2>/dev/null | sed 's/^/    /'
                        tail -10 resp1.log | sed 's/^/    /'; die "RESP stage 1 failed"; }
resp -O -i LIG.respin2 -o LIG.respout2 -p LIG.esout2 -q LIG.qout1 -t LIG.qout2 -e LIG_all.esp > resp2.log 2>&1
[[ -s LIG.qout2 ]] || { tail -25 LIG.respout2 2>/dev/null | sed 's/^/    /'; die "RESP stage 2 failed"; }
echo "  ✓ qout1 / qout2"
echo
echo "  fit quality (also the .esp format check — a value near 1 means the"
echo "  ESP file was misread; tell me if you see that):"
grep -iE "rms" LIG.respout1 | head -3 | sed 's/^/    stage1: /'
grep -iE "rms" LIG.respout2 | head -3 | sed 's/^/    stage2: /'
echo
echo "  stage-2 charges:"
tr -s ' \t' '\n' < LIG.qout2 | grep -v '^$' \
  | awk '{printf "    %2d  %+9.5f\n", NR, $1; s+=$1} END{printf "    SUM %+9.5f   (must be -2.00000)\n", s}'

antechamber -i LIG.ac -fi ac -o LIG_resp.mol2 -fo mol2 \
            -c rc -cf LIG.qout2 -at gaff2 -rn LIG > ac_resp.log 2>&1
[[ -s LIG_resp.mol2 ]] || { tail -15 ac_resp.log | sed 's/^/    /'; die "LIG_resp.mol2 not written"; }
echo "  ✓ LIG_resp.mol2"

# ── 5. AM1-BCC on the SAME geometry ─────────────────────────────────────
# maxcyc=0 keeps sqm at a single point. Without it sqm minimises in vacuum
# and transfers the ammonium proton to the carboxylate — the same class of
# error that caused this whole rerun — and the two arms would no longer be
# comparable.
echo; echo "$HR"; echo " 5. AM1-BCC arm, identical geometry"; echo "$HR"
antechamber -i LIG_optA.sdf -fi "$BUILT" -o LIG_bcc.mol2 -fo mol2 \
            -c bcc -nc -2 -at gaff2 -rn LIG -s 2 \
            -ek "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0, ndiis_attempts=700" \
            > ac_bcc.log 2>&1
if [[ ! -s LIG_bcc.mol2 ]]; then
    echo "  failed. Tail of ac_bcc.log:"; tail -25 ac_bcc.log | sed 's/^/    /'
    echo "  (if sqm objects to maxcyc=0, tell me — I will handle the geometry differently)"
    die "AM1-BCC arm failed"
fi
echo "  ✓ LIG_bcc.mol2"
awk '/@<TRIPOS>ATOM/{f=1;next}/@<TRIPOS>BOND/{f=0}f{printf "    %2d %-5s %-5s %+9.5f\n",$1,$2,$6,$9; s+=$9}
     END{printf "    SUM %+9.5f\n",s}' LIG_bcc.mol2

# ── 6. bonded parameters must match between arms ─────────────────────────
echo; echo "$HR"; echo " 6. frcmod check"; echo "$HR"
parmchk2 -i LIG_bcc.mol2  -f mol2 -o LIG_bcc.frcmod  -s gaff2 > pc_bcc.log  2>&1
parmchk2 -i LIG_resp.mol2 -f mol2 -o LIG_resp.frcmod -s gaff2 > pc_resp.log 2>&1
if diff -q <(grep -v '^Remark' LIG_bcc.frcmod 2>/dev/null) \
           <(grep -v '^Remark' LIG_resp.frcmod 2>/dev/null) >/dev/null 2>&1; then
    echo "  ✓ identical — only the charges differ between the two arms"
else
    echo "  [!] frcmod files differ:"
    diff <(grep -v '^Remark' LIG_bcc.frcmod) <(grep -v '^Remark' LIG_resp.frcmod) 2>/dev/null \
        | head -20 | sed 's/^/    /'
fi
for f in LIG_bcc.frcmod LIG_resp.frcmod; do
    grep -q ATTN "$f" 2>/dev/null && { echo "  [!] $f has ATTN parameters:"; grep ATTN "$f" | sed 's/^/    /'; }
done

echo; echo "$HR"
echo " DONE. Send me:"
echo "   cp3_report.txt  LIG_bcc.mol2  LIG_resp.mol2  LIG.qout2"
echo "   LIG.respout1  LIG.respout2  LIG.ac"
echo "$HR"
