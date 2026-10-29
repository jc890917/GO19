#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════
#  cp4_build_wt.sh — build the corrected wild-type system
#
#  Run on Yoltla:   bash cp4_build_wt.sh 2>&1 | tee cp4_build.log
#
#  What changes relative to the published build:
#    * LIG.mol2 is the CORRECTED glyphosate dianion (NH2+, COO-, PO3(2-)),
#      AM1-BCC charges, coordinates from the seed404 rank-1 docked pose —
#      the centroid of the catalytic ensemble across 10 docking runs.
#      Verified: C(reacting)-FAD N5 = 3.64 A, P-N5 = 6.57 A, no clash < 2.8 A.
#    * a Solute index group is created and the production mdp writes only
#      solute coordinates, cutting each 100 ns trajectory from ~5 GB to ~45 MB.
#  Unchanged: protein_clean.pdb (protein only — the spurious FAD N1-H never
#  entered the MD), FAD_bcc.mol2, FAD.frcmod, ff19SB/GAFF2/OPC, box, ions.
#
#  NOTE on the earlier failures: AmberTools binaries here need the Intel
#  runtime (libimf.so). config.sh's load_amber() loads it; exporting PATH by
#  hand does not. This script uses load_amber() and then VERIFIES the tools
#  actually execute before doing anything.
# ══════════════════════════════════════════════════════════════════════════
set -uo pipefail
HR="──────────────────────────────────────────────────────────────────"
die(){ echo; echo "[STOP] $*"; echo "Send me cp4_build.log and any *.log in the build folder."; exit 1; }

WT_OLD="/LUSTRE/home/lancad/2026/l.100073/alumnos/ian/go19"
DEST="/LUSTRE/home/lancad/2026/l.100073/alumnos/ian/go19_v2"
HERE="$(cd "$(dirname "$0")" && pwd)"

echo "$HR"; echo " 0. inputs"; echo "$HR"
[[ -d "$WT_OLD" ]] || die "cannot find the original build at $WT_OLD"
[[ -f "$HERE/LIG.mol2" ]] || die "LIG.mol2 must sit next to this script"
NEED_OLD=(config.sh protein_clean.pdb FAD_bcc.mol2 FAD.frcmod
          05_build_system.sh 06_minimize.sh 07_equilibrate.sh 08_production_100ns.sh)
for f in "${NEED_OLD[@]}"; do [[ -f "$WT_OLD/$f" ]] || die "missing $WT_OLD/$f"; done
# LIG.frcmod comes from the CP3 run (parmchk2 on LIG_bcc.mol2)
FRC=""
for c in "$HERE/LIG.frcmod" "$HERE/LIG_bcc.frcmod" "$HOME/cp3_charges/LIG_bcc.frcmod"; do
    [[ -f "$c" ]] && { FRC="$c"; break; }
done
[[ -n "$FRC" ]] || die "cannot find LIG_bcc.frcmod (or LIG.frcmod). It was produced by 05c_yoltla_resp.sh in ~/cp3_charges."
echo "  original build : $WT_OLD"
echo "  new ligand     : $HERE/LIG.mol2"
echo "  ligand frcmod  : $FRC"
echo "  destination    : $DEST"

echo; echo "$HR"; echo " 1. modules, and proof the tools run"; echo "$HR"
cd "$WT_OLD" && source ./config.sh && load_amber
for b in tleap parmchk2 antechamber; do
    P=$(command -v "$b" 2>/dev/null) || die "$b not on PATH after load_amber"
    OUT=$("$b" -h 2>&1 </dev/null | head -3)
    grep -qi "error while loading shared libraries\|cannot open shared object" <<<"$OUT" \
        && { echo "$OUT" | sed 's/^/    /'; die "$b cannot load its libraries (libimf.so?)"; }
    printf "  ✓ %-12s %s\n" "$b" "$P"
done
python3 -c "import parmed" 2>/dev/null && echo "  ✓ parmed importable" || die "parmed not importable"

echo; echo "$HR"; echo " 2. verify the new ligand BEFORE building"; echo "$HR"
NAT=$(awk '/@<TRIPOS>ATOM/{f=1;next}/@<TRIPOS>BOND/{f=0}f&&NF{n++}END{print n+0}' "$HERE/LIG.mol2")
QSUM=$(awk '/@<TRIPOS>ATOM/{f=1;next}/@<TRIPOS>BOND/{f=0}f&&NF{s+=$9}END{printf "%.4f", s}' "$HERE/LIG.mol2")
echo "  atoms: $NAT   charge sum: $QSUM"
[[ "$NAT" == "16" ]] || die "LIG.mol2 has $NAT atoms, expected 16"
awk -v q="$QSUM" 'BEGIN{exit !(q < -1.98 && q > -2.02)}' || die "charge sum $QSUM is not -2"
# protonation: exactly 2 hydrogens bonded to N, none to O
python3 - "$HERE/LIG.mol2" << 'PYCHK'
import sys, math
p=sys.argv[1]; sec=None; at=[]; bd=[]
for ln in open(p):
    if ln.startswith("@<TRIPOS>"): sec=ln.split(">")[1].strip(); continue
    t=ln.split()
    if sec=="ATOM" and len(t)>=9: at.append((t[1],t[5],float(t[2]),float(t[3]),float(t[4])))
    elif sec=="BOND" and len(t)>=4: bd.append((int(t[1])-1,int(t[2])-1))
el=lambda s: "P" if s.split(".")[0].upper().startswith("P") else s[0].upper()
nH=sum(1 for a,b in bd if {el(at[a][1]),el(at[b][1])}=={"N","H"})
oH=sum(1 for a,b in bd if {el(at[a][1]),el(at[b][1])}=={"O","H"})
print(f"  H bonded to N: {nH} (need 2)   H bonded to O: {oH} (need 0)")
if nH!=2 or oH!=0:
    print("  *** WRONG PROTONATION STATE — this is the error being corrected. ABORT. ***"); sys.exit(1)
print("  >>> protonation correct: ammonium, carboxylate, phosphonate dianion")
PYCHK
[[ $? -eq 0 ]] || die "ligand protonation check failed"

echo; echo "$HR"; echo " 3. assemble the build folder"; echo "$HR"
if [[ -d "$DEST/gromacs" ]]; then
    echo "  $DEST/gromacs already exists — refusing to overwrite."
    echo "  Move it aside first if you really want to rebuild."; exit 0
fi
mkdir -p "$DEST"
for f in "${NEED_OLD[@]}"; do cp "$WT_OLD/$f" "$DEST/"; done
cp "$HERE/LIG.mol2" "$DEST/LIG.mol2"
cp "$FRC" "$DEST/LIG.frcmod"
# record what the old ligand was, for the data deposit
cp "$WT_OLD/LIG.mol2" "$DEST/LIG_original_submitted.mol2" 2>/dev/null || true
sed -i "s#^WORKDIR=.*#WORKDIR=\"$DEST\"#" "$DEST/config.sh"
sed -i "s#^cd /LUSTRE.*#cd \"\$SLURM_SUBMIT_DIR\"#" "$DEST/06_minimize.sh" "$DEST/07_equilibrate.sh"
grep -q 'cd "\$SLURM_SUBMIT_DIR"' "$DEST/06_minimize.sh"    || die "could not patch cd in 06_minimize.sh"
grep -q 'cd "\$SLURM_SUBMIT_DIR"' "$DEST/07_equilibrate.sh" || die "could not patch cd in 07_equilibrate.sh"
echo "  ✓ folder assembled, config.sh and 06/07 paths patched"

echo; echo "$HR"; echo " 4. tleap + parmed (05_build_system.sh)"; echo "$HR"
( cd "$DEST" && bash 05_build_system.sh > build.log 2>&1 )
RC=$?
if [[ ! -f "$DEST/gromacs/complex.gro" ]]; then
    echo "  build failed (exit $RC). Tail of build.log:"; tail -40 "$DEST/build.log" | sed 's/^/    /'
    echo; echo "  tleap errors, if any:"; grep -i "error\|fatal\|could not" "$DEST"/*leap*.log 2>/dev/null | head -15 | sed 's/^/    /'
    die "05_build_system.sh did not produce gromacs/complex.gro"
fi
NATOMS=$(sed -n '2p' "$DEST/gromacs/complex.gro" | tr -d ' ')
CHG=$(grep -i "Total unperturbed charge" "$DEST"/leap.log "$DEST"/tleap.log 2>/dev/null | tail -1)
echo "  ✓ system built: $NATOMS atoms"
echo "    final $CHG"
awk -v s="$(sed 's/.*: *//' <<<"$CHG")" 'BEGIN{exit !(s<0.01 && s>-0.01)}' 2>/dev/null \
    || echo "    [!] confirm this is ~0 before equilibrating"
grep -iE "WARNING|ERROR" "$DEST"/leap.log 2>/dev/null | head -10 | sed 's/^/    leap: /' || true

echo; echo "$HR"; echo " 5. Solute index group + trajectory-size fix"; echo "$HR"
cd "$DEST/gromacs"
load_gromacs 2>/dev/null || true
GMX="${GMX:-gmx_mpi}"; command -v "$GMX" >/dev/null || GMX=gmx
command -v "$GMX" >/dev/null || die "no gmx on PATH after load_gromacs"
echo "  residue names present:"
awk 'NR>2 && NF>3 {print substr($0,6,5)}' complex.gro | sed 's/ //g' | sort -u | tr '\n' ' ' | fold -w 100 | sed 's/^/    /'
echo
printf 'q\n' | "$GMX" make_ndx -f complex.gro -o index_base.ndx > mkndx.log 2>&1
"$GMX" select -s complex.gro -on solute.ndx \
       -select '"Solute" not resname WAT SOL HOH TIP3 OPC NA CL "Na+" "Cl-" K "K+"' \
       > select.log 2>&1
if [[ ! -s solute.ndx ]]; then
    echo "  gmx select failed:"; tail -15 select.log | sed 's/^/    /'
    die "could not build the Solute group — check the residue names listed above"
fi
NSOL=$(grep -v '^\[' solute.ndx | tr -s ' \n' '\n' | grep -c '[0-9]')
echo "  Solute group: $NSOL atoms  (expected ~5,800 = protein + FAD + LIG)"
if (( NSOL < 4500 || NSOL > 7500 )); then
    die "Solute group size $NSOL is implausible — the resname exclusion is wrong"
fi
cat index_base.ndx solute.ndx > index.ndx
echo "  ✓ index.ndx written with a Solute group"
cd "$DEST"
if grep -q '^compressed-x-grps' 08_production_100ns.sh; then
    echo "  08_production_100ns.sh already sets compressed-x-grps"
else
    sed -i 's|^\(nstxout-compressed  = 5000.*\)$|\1\ncompressed-x-grps   = Solute   ; solute only: ~45 MB per 100 ns instead of ~5 GB|' \
        08_production_100ns.sh
    grep -q '^compressed-x-grps' 08_production_100ns.sh \
        && echo "  ✓ compressed-x-grps = Solute added to the production mdp" \
        || echo "  [!] could not patch the mdp automatically — add 'compressed-x-grps = Solute' by hand"
fi
grep -q '\-n .*index.ndx\|GMX_NDX' 08_production_100ns.sh \
    || echo "  [!] check that 08_production_100ns.sh passes '-n gromacs/index.ndx' to grompp"

echo; echo "$HR"; echo " DONE — nothing has been submitted."; echo "$HR"
echo " Built in: $DEST"
echo " Before equilibrating, send me:"
echo "   cp4_build.log   $DEST/build.log   $DEST/leap.log"
echo " I will check the charge, atom count, ion counts and the ligand's"
echo " position in the built system against the docked pose."
