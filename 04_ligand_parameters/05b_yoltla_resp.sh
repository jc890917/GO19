#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════
#  05b_yoltla_resp.sh  — replaces 05_yoltla_resp.sh
#
#  Fixes in this version:
#    * NO 'set -e'. The previous script died silently at an unguarded
#      antechamber call, before its own error handler could report anything.
#      Every external command is now checked explicitly and prints why it
#      failed.
#    * '-fi mdl' instead of '-fi sdf'. AmberTools 19 names the MDL/SDF
#      reader 'mdl'; 'sdf' is not a recognised format there.
#    * dropped the '-c wc' attempt, which cannot work on an input carrying
#      no charges. respgen needs only connectivity and atom types, so the
#      .ac is built with no charge method at all.
#
#  Inputs in the current directory: LIG_all.esp  LIG_optA.sdf  nconf.txt
#  Run:  bash 05b_yoltla_resp.sh 2>&1 | tee cp3_report.txt
# ══════════════════════════════════════════════════════════════════════════
set -uo pipefail
AMBBIN=/LUSTRE/yoltla/programas/amber/ambertools19_INTEL-parallel_studio_xe_2017_update4_IMPI/amber18/bin
export PATH="$AMBBIN:$PATH"
export AMBERHOME="$(dirname "$AMBBIN")"
HR="──────────────────────────────────────────────────────────────────"
die () { echo; echo "[STOP] $*"; echo "Send me cp3_report.txt and every *.log in this directory."; exit 1; }

echo "$HR"; echo " environment"; echo "$HR"
for b in antechamber respgen espgen resp parmchk2; do
    command -v "$b" >/dev/null || die "$b not found on PATH"
    printf "  %-12s %s\n" "$b" "$(command -v $b)"
done
echo "  AMBERHOME = $AMBERHOME"
for f in LIG_all.esp LIG_optA.sdf nconf.txt; do
    [[ -s "$f" ]] || die "missing or empty: $f"
done
NCONF=$(tr -dc '0-9' < nconf.txt)
NATM=$(awk 'NR==1{print $1+0; exit}' LIG_all.esp)
NPTS=$(awk 'NR==1{print $2+0; exit}' LIG_all.esp)
echo "  LIG_all.esp : $NCONF conformer(s), $NATM atoms, $NPTS ESP points"
[[ "$NATM" == "16" ]] || die "expected 16 atoms in LIG_all.esp, header says $NATM"

# ── 1. topology + GAFF2 atom types (no charge method) ────────────────────
echo; echo "$HR"; echo " 1. atom types via antechamber"; echo "$HR"
BUILT=""
for FMT in mdl sdf; do
    rm -f LIG.ac
    echo "  trying -fi $FMT ..."
    antechamber -i LIG_optA.sdf -fi "$FMT" -o LIG.ac -fo ac \
                -at gaff2 -nc -2 -rn LIG > "antechamber_ac_$FMT.log" 2>&1
    RC=$?
    if [[ -s LIG.ac ]]; then
        echo "  ✓ LIG.ac written with -fi $FMT"; BUILT="$FMT"; break
    fi
    echo "    failed (exit $RC). Last lines:"
    tail -12 "antechamber_ac_$FMT.log" | sed 's/^/      /'
done
if [[ -z "$BUILT" ]]; then
    echo; echo "  Neither format worked. Available antechamber formats:"
    antechamber -L 2>&1 | head -30 | sed 's/^/    /'
    die "could not build LIG.ac from LIG_optA.sdf"
fi
echo
echo "  GAFF2 atom types:"
awk '/^ATOM/{printf "    %2d  %-5s %-5s\n", $2, $3, $NF}' LIG.ac
NAC=$(grep -c '^ATOM' LIG.ac)
[[ "$NAC" == "16" ]] || die "LIG.ac has $NAC atoms, expected 16"

# ── 2. two-stage RESP ────────────────────────────────────────────────────
echo; echo "$HR"; echo " 2. two-stage RESP (AmberTools reference implementation)"; echo "$HR"
respgen -i LIG.ac -o LIG.respin1 -f resp1 -n "$NCONF" > respgen1.log 2>&1
[[ -s LIG.respin1 ]] || { tail -15 respgen1.log | sed 's/^/    /'; die "respgen stage 1 failed"; }
respgen -i LIG.ac -o LIG.respin2 -f resp2 -n "$NCONF" > respgen2.log 2>&1
[[ -s LIG.respin2 ]] || { tail -15 respgen2.log | sed 's/^/    /'; die "respgen stage 2 failed"; }
echo "  ✓ respin1 / respin2 written"

resp -O -i LIG.respin1 -o LIG.respout1 -p LIG.esout1 -t LIG.qout1 -e LIG_all.esp > resp1.log 2>&1
[[ -s LIG.qout1 ]] || { echo "  resp stage 1 output tail:"; tail -25 LIG.respout1 2>/dev/null | sed 's/^/    /'
                        tail -10 resp1.log | sed 's/^/    /'; die "RESP stage 1 failed"; }
resp -O -i LIG.respin2 -o LIG.respout2 -p LIG.esout2 -q LIG.qout1 -t LIG.qout2 -e LIG_all.esp > resp2.log 2>&1
[[ -s LIG.qout2 ]] || { tail -25 LIG.respout2 2>/dev/null | sed 's/^/    /'; die "RESP stage 2 failed"; }
echo "  ✓ qout1 / qout2 written"

echo
echo "  fit quality — this is also the .esp format check:"
grep -iE "rms" LIG.respout1 | head -3 | sed 's/^/    stage1: /'
grep -iE "rms" LIG.respout2 | head -3 | sed 's/^/    stage2: /'
echo "    A sane RESP fit gives a small relative RMS. A value near or above 1"
echo "    would mean the ESP file was misread — tell me if you see that."
echo
echo "  stage-2 charges (LIG.qout2):"
cat LIG.qout2 | tr -s ' ' '\n' | grep -v '^$' | awk '{printf "    %2d  %+9.5f\n", NR, $1; s+=$1} END{printf "    SUM %+9.5f  (must be -2.00000)\n", s}'

antechamber -i LIG.ac -fi ac -o LIG_resp.mol2 -fo mol2 \
            -c rc -cf LIG.qout2 -at gaff2 -rn LIG > antechamber_resp.log 2>&1
[[ -s LIG_resp.mol2 ]] || { tail -15 antechamber_resp.log | sed 's/^/    /'; die "LIG_resp.mol2 not written"; }
echo "  ✓ LIG_resp.mol2"

# ── 3. AM1-BCC on the SAME geometry ──────────────────────────────────────
# maxcyc=0 keeps sqm at a single point. Without it sqm minimises in vacuum
# and transfers the ammonium proton, so the two arms would no longer be
# comparable — the same class of error that caused this whole rerun.
echo; echo "$HR"; echo " 3. AM1-BCC arm, identical geometry, no internal optimisation"; echo "$HR"
antechamber -i LIG_optA.sdf -fi "$BUILT" -o LIG_bcc.mol2 -fo mol2 \
            -c bcc -nc -2 -at gaff2 -rn LIG -s 2 \
            -ek "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0, ndiis_attempts=700" \
            > antechamber_bcc.log 2>&1
if [[ ! -s LIG_bcc.mol2 ]]; then
    echo "  failed. Last 25 lines of antechamber_bcc.log:"
    tail -25 antechamber_bcc.log | sed 's/^/    /'
    echo "  (if sqm complains about maxcyc=0, tell me — I will pre-optimise instead)"
    die "AM1-BCC arm failed"
fi
echo "  ✓ LIG_bcc.mol2"
awk '/@<TRIPOS>ATOM/{f=1;next}/@<TRIPOS>BOND/{f=0}f{printf "    %2d %-5s %-5s %+9.5f\n",$1,$2,$6,$9; s+=$9}END{printf "    SUM %+9.5f\n",s}' LIG_bcc.mol2

# ── 4. bonded parameters must be identical between the two arms ───────────
echo; echo "$HR"; echo " 4. frcmod check"; echo "$HR"
parmchk2 -i LIG_bcc.mol2  -f mol2 -o LIG_bcc.frcmod  -s gaff2 > parmchk2_bcc.log  2>&1
parmchk2 -i LIG_resp.mol2 -f mol2 -o LIG_resp.frcmod -s gaff2 > parmchk2_resp.log 2>&1
if diff -q <(grep -v '^Remark' LIG_bcc.frcmod) <(grep -v '^Remark' LIG_resp.frcmod) >/dev/null 2>&1; then
    echo "  ✓ identical — only charges differ between the two arms"
else
    echo "  [!] frcmod files differ:"
    diff <(grep -v '^Remark' LIG_bcc.frcmod) <(grep -v '^Remark' LIG_resp.frcmod) | head -20 | sed 's/^/    /'
fi
for f in LIG_bcc.frcmod LIG_resp.frcmod; do
    grep -q ATTN "$f" 2>/dev/null && { echo "  [!] $f has ATTN parameters:"; grep ATTN "$f" | sed 's/^/    /'; }
done

echo; echo "$HR"
echo " DONE. Send me:"
echo "   cp3_report.txt  LIG_bcc.mol2  LIG_resp.mol2  LIG.qout2"
echo "   LIG.respout1  LIG.respout2  LIG.ac"
echo "$HR"
