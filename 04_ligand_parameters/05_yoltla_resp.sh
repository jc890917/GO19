#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════
#  05_yoltla_resp.sh — the official RESP fit + the AM1-BCC arm, on Yoltla.
#
#  Gaussian is NOT licensed to this account (both g09 install directories
#  are Permission denied for groups l.100073 / lancad / tl.100073), so the
#  quantum chemistry was done elsewhere with PySCF at the same level of
#  theory: HF/6-31G* IEF-PCM(water) optimisation, then a vacuum HF/6-31G*
#  Merz-Kollman ESP. This script does the parts that need AmberTools, whose
#  `resp` is the reference implementation of the fit.
#
#  Inputs it expects in the same directory:
#      LIG_all.esp     multi-conformer ESP (Amber format, coords in Bohr)
#      LIG_optA.sdf    the OPTIMISED conformer A geometry
#      nconf.txt       how many conformers are in LIG_all.esp
#  Run:  bash 05_yoltla_resp.sh
# ══════════════════════════════════════════════════════════════════════════
set -euo pipefail
AMBBIN=/LUSTRE/yoltla/programas/amber/ambertools19_INTEL-parallel_studio_xe_2017_update4_IMPI/amber18/bin
export PATH="$AMBBIN:$PATH"; export AMBERHOME="$(dirname "$AMBBIN")"
for b in antechamber respgen espgen resp parmchk2; do
    command -v "$b" >/dev/null || { echo "[ERROR] $b not found"; exit 1; }
done
for f in LIG_all.esp LIG_optA.sdf nconf.txt; do
    [[ -f "$f" ]] || { echo "[ERROR] missing $f"; exit 1; }
done
NCONF=$(tr -dc '0-9' < nconf.txt)
echo "✓ AmberTools found; LIG_all.esp carries $NCONF conformer(s)"

# ── topology + GAFF2 atom types, from the optimised geometry ──────────────
antechamber -i LIG_optA.sdf -fi sdf -o LIG.ac -fo ac -at gaff2 -nc -2 -rn LIG \
            -c wc -cf LIG_zero.dat > antechamber_ac.log 2>&1 || true
if [[ ! -f LIG.ac ]]; then
    # no-charge route: let antechamber assign types only
    antechamber -i LIG_optA.sdf -fi sdf -o LIG.ac -fo ac -at gaff2 -nc -2 -rn LIG \
                > antechamber_ac.log 2>&1
fi
[[ -f LIG.ac ]] || { echo "[ERROR] could not build LIG.ac"; tail -30 antechamber_ac.log; exit 1; }
echo ""
echo "GAFF2 atom types assigned:"
awk '/^ATOM/{printf "  %2d  %-5s %-5s\n", $2, $3, $NF}' LIG.ac

# ── two-stage RESP with AmberTools ───────────────────────────────────────
respgen -i LIG.ac -o LIG.respin1 -f resp1 -n "$NCONF"
respgen -i LIG.ac -o LIG.respin2 -f resp2 -n "$NCONF"

resp -O -i LIG.respin1 -o LIG.respout1 -p LIG.esout1 -t LIG.qout1 -e LIG_all.esp
resp -O -i LIG.respin2 -o LIG.respout2 -p LIG.esout2 -q LIG.qout1 -t LIG.qout2 -e LIG_all.esp
[[ -f LIG.qout2 ]] || { echo "[ERROR] RESP stage 2 failed"; tail -40 LIG.respout2; exit 1; }

# RRMS is the format check: a mis-read .esp gives a hopeless fit.
echo ""
echo "── fit quality (this is also the .esp format check) ──"
grep -iE "relative rms|rms" LIG.respout1 | head -4 | sed 's/^/  stage1: /'
grep -iE "relative rms|rms" LIG.respout2 | head -4 | sed 's/^/  stage2: /'
echo "  A good RESP fit has relative RMS around 0.1-0.3. Anything near or"
echo "  above 1.0 means the ESP file was misread — stop and tell me."

antechamber -i LIG.ac -fi ac -o LIG_resp.mol2 -fo mol2 \
            -c rc -cf LIG.qout2 -at gaff2 -rn LIG > antechamber_resp.log 2>&1
[[ -f LIG_resp.mol2 ]] || { echo "[ERROR] LIG_resp.mol2 not written"; tail -30 antechamber_resp.log; exit 1; }
echo "  ✓ LIG_resp.mol2"

# ── AM1-BCC arm, SAME geometry, no internal optimisation ─────────────────
# maxcyc=0 keeps sqm at a single point; otherwise it minimises in vacuum and
# transfers the ammonium proton, and the two arms stop being comparable.
echo ""
echo "── AM1-BCC arm ──"
antechamber -i LIG_optA.sdf -fi sdf -o LIG_bcc.mol2 -fo mol2 \
            -c bcc -nc -2 -at gaff2 -rn LIG -s 2 \
            -ek "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0, ndiis_attempts=700" \
            > antechamber_bcc.log 2>&1
[[ -f LIG_bcc.mol2 ]] || { echo "[ERROR] AM1-BCC arm failed"; tail -40 antechamber_bcc.log; exit 1; }
echo "  ✓ LIG_bcc.mol2"

parmchk2 -i LIG_bcc.mol2  -f mol2 -o LIG_bcc.frcmod  -s gaff2
parmchk2 -i LIG_resp.mol2 -f mol2 -o LIG_resp.frcmod -s gaff2
echo ""
if diff -q <(grep -v '^Remark' LIG_bcc.frcmod) <(grep -v '^Remark' LIG_resp.frcmod) >/dev/null; then
    echo "  ✓ frcmod files identical — only charges differ"
else
    echo "  [!] frcmod files differ:"; diff <(grep -v '^Remark' LIG_bcc.frcmod) <(grep -v '^Remark' LIG_resp.frcmod) | head -20
fi
for f in LIG_bcc.frcmod LIG_resp.frcmod; do
    grep -q ATTN "$f" && { echo "  [!] $f has ATTN parameters:"; grep ATTN "$f"; }
done

echo ""
echo "══════════════════════════════════════════════════════════════════════"
echo " Done. Send me back:  LIG_bcc.mol2  LIG_resp.mol2  LIG.qout2"
echo "                      LIG.respout1  LIG.respout2  LIG.ac"
echo " I will run the comparison and the numpy cross-check against them."
echo "══════════════════════════════════════════════════════════════════════"
