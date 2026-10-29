#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════
#  cp5_build_mutants.sh — build R53A, R272L and R53A/R272L from go19_v2
#
#  Run on the Yoltla LOGIN NODE (no sbatch, no GPU, minutes per system):
#      bash cp5_build_mutants.sh 2>&1 | tee cp5_build.log
#
#  Everything is inherited from the corrected wild type in go19_v2: the
#  glyphosate dianion with AM1-BCC charges at the seed404 docked pose, FAD,
#  ff19SB/GAFF2/OPC, the 12 A box, 0.15 M NaCl, and the patched run scripts.
#  The ONLY difference between each mutant and the WT is the point mutation.
#
#  Mutations are applied to protein_clean.pdb by deleting the side-chain
#  atoms the target residue does not have; tleap rebuilds the rest from its
#  own residue template.
#      ARG -> ALA : keep N CA C O CB
#      ARG -> LEU : keep N CA C O CB CG, rename CD -> CD1  (tleap adds CD2)
#
#  Gates, learned the hard way during CP3/CP4:
#    * module load via config.sh load_amber(), then PROVE tleap executes —
#      these binaries need the Intel runtime (libimf.so) and a bare PATH
#      export is not enough
#    * no 'set -e': every step reports its own failure with a log tail
#    * abort unless residue 53/272 really is ARG before touching it
#    * the mutant protein charge must come out an EXACT integer (-5 or -6);
#      a fractional value means the mutation left a broken residue
#    * the Solute index group is built by MERGING make_ndx groups, never by
#      'gmx select ... resname "Na+"' — GROMACS treats quoted strings as
#      regular expressions, so "Na+" means "N followed by one or more a"
#      and silently fails to match the residue name Na+
# ══════════════════════════════════════════════════════════════════════════
set -uo pipefail
HR="──────────────────────────────────────────────────────────────────"
die(){ echo; echo "[STOP] $*"; echo "Send me cp5_build.log and the build.log of the system that failed."; exit 1; }

HERE="$(cd "$(dirname "$0")" && pwd)"
IAN="/LUSTRE/home/lancad/2026/l.100073/alumnos/ian"
WT="$IAN/go19_v2"
MUTANTS=("R53A" "R272L" "R53A_R272L")
PARTITION="q1d-20p"          # the gpus partition has 16 of 24 nodes down
NCPU=20                      # q*-20p nodes have 20 cores

echo "$HR"; echo " 0. inputs"; echo "$HR"
[[ -d "$WT" ]] || die "corrected wild type not found at $WT — run CP4 first"
NEED=(config.sh protein_clean.pdb LIG.mol2 LIG.frcmod FAD_bcc.mol2 FAD.frcmod
      05_build_system.sh 06_minimize.sh 07_equilibrate.sh 08_production_100ns.sh)
for f in "${NEED[@]}"; do [[ -f "$WT/$f" ]] || die "missing $WT/$f"; done
[[ -f "$HERE/prmtop_charge.py" ]] || die "prmtop_charge.py must sit next to this script"
echo "  source wild type : $WT"
echo "  mutants          : ${MUTANTS[*]}"
echo "  parent directory : $IAN"

echo; echo "$HR"; echo " 1. modules, and proof the tools run"; echo "$HR"
cd "$WT" && source ./config.sh && load_amber
for b in tleap parmchk2; do
    P=$(command -v "$b" 2>/dev/null) || die "$b not on PATH after load_amber"
    OUT=$("$b" -h 2>&1 </dev/null | head -3)
    grep -qi "error while loading shared libraries\|cannot open shared object" <<<"$OUT" \
        && { echo "$OUT" | sed 's/^/    /'; die "$b cannot load its libraries (libimf.so?)"; }
    printf "  ✓ %-10s %s\n" "$b" "$P"
done
python3 -c "import parmed" 2>/dev/null && echo "  ✓ parmed importable" || die "parmed not importable"

echo; echo "$HR"; echo " 2. the inherited ligand is the corrected one"; echo "$HR"
NAT=$(awk '/@<TRIPOS>ATOM/{f=1;next}/@<TRIPOS>BOND/{f=0}f&&NF{n++}END{print n+0}' "$WT/LIG.mol2")
QSUM=$(awk '/@<TRIPOS>ATOM/{f=1;next}/@<TRIPOS>BOND/{f=0}f&&NF{s+=$9}END{printf "%.4f",s}' "$WT/LIG.mol2")
echo "  LIG.mol2: $NAT atoms, charge $QSUM"
[[ "$NAT" == "16" ]] || die "inherited LIG.mol2 has $NAT atoms, expected 16"
python3 - "$WT/LIG.mol2" << 'PYCHK'
import sys
p=sys.argv[1]; sec=None; at=[]; bd=[]
for ln in open(p):
    if ln.startswith("@<TRIPOS>"): sec=ln.split(">")[1].strip(); continue
    t=ln.split()
    if sec=="ATOM" and len(t)>=9: at.append(t[5])
    elif sec=="BOND" and len(t)>=4: bd.append((int(t[1])-1,int(t[2])-1))
el=lambda s:"P" if s.split(".")[0].upper().startswith("P") else s[0].upper()
nH=sum(1 for a,b in bd if {el(at[a]),el(at[b])}=={"N","H"})
oH=sum(1 for a,b in bd if {el(at[a]),el(at[b])}=={"O","H"})
print(f"  H on N: {nH} (need 2)   H on O: {oH} (need 0)")
sys.exit(0 if (nH==2 and oH==0) else 1)
PYCHK
[[ $? -eq 0 ]] || die "the inherited ligand has the WRONG protonation state"
echo "  >>> ammonium / carboxylate / phosphonate dianion — correct"

# ── mutation helper ──────────────────────────────────────────────────────
mutate () {   # mutate IN OUT RESNUM NEWRES
    local IN="$1" OUT="$2" RES="$3" NEW="$4" OLD
    OLD=$(awk -v r="$RES" '($1=="ATOM"||$1=="HETATM") && substr($0,23,4)+0==r {print substr($0,18,3); exit}' "$IN")
    [[ "$OLD" == "ARG" ]] || { echo "    residue $RES is '$OLD', expected ARG"; return 1; }
    awk -v r="$RES" -v new="$NEW" '
      function trim(s){ gsub(/ /,"",s); return s }
      {
        if (($1=="ATOM"||$1=="HETATM") && substr($0,23,4)+0==r) {
          name = trim(substr($0,13,4)); keep = 0
          if (new=="ALA" && (name=="N"||name=="CA"||name=="C"||name=="O"||name=="CB")) keep=1
          if (new=="LEU" && (name=="N"||name=="CA"||name=="C"||name=="O"||name=="CB"||name=="CG"||name=="CD")) keep=1
          if (!keep) next
          line = substr($0,1,17) new substr($0,21)
          if (new=="LEU" && name=="CD") line = substr(line,1,12) " CD1" substr(line,17)
          print line
        } else print
      }' "$IN" > "$OUT"
    local NKEPT
    NKEPT=$(awk -v r="$RES" '($1=="ATOM") && substr($0,23,4)+0==r' "$OUT" | wc -l)
    echo "    ARG$RES -> $NEW   ($NKEPT atoms kept, tleap rebuilds the rest)"
    return 0
}

declare -A EXPECTED_CHARGE=( [R53A]=-5 [R272L]=-5 [R53A_R272L]=-6 )
SUMMARY=""

for M in "${MUTANTS[@]}"; do
    D="$IAN/go19_v2_$M"
    echo; echo "$HR"; echo " ══ $M ══  -> $D"; echo "$HR"
    ALREADY=0
    if [[ -f "$D/gromacs/complex.gro" ]]; then
        echo "  already built — skipping tleap, but re-running every check and"
        echo "  completing the index group and mdp patch."
        ALREADY=1
    fi
    mkdir -p "$D"
  if [[ "$ALREADY" -eq 0 ]]; then
    for f in "${NEED[@]}"; do cp "$WT/$f" "$D/"; done
    cp "$WT/LIG_original_submitted.mol2" "$D/" 2>/dev/null || true
    cp "$D/protein_clean.pdb" "$D/protein_clean_WT.pdb"

    echo "  applying the mutation:"
    case "$M" in
        R53A)       mutate "$D/protein_clean_WT.pdb" "$D/protein_clean.pdb" 53  ALA || die "$M mutation failed" ;;
        R272L)      mutate "$D/protein_clean_WT.pdb" "$D/protein_clean.pdb" 272 LEU || die "$M mutation failed" ;;
        R53A_R272L) mutate "$D/protein_clean_WT.pdb" "$D/tmp.pdb"            53  ALA || die "$M mutation 1 failed"
                    mutate "$D/tmp.pdb"              "$D/protein_clean.pdb" 272 LEU || die "$M mutation 2 failed"
                    rm -f "$D/tmp.pdb" ;;
    esac
    echo "  residue identity after mutation:"
    awk '$1=="ATOM" && substr($0,13,4)~/ CA /{r=substr($0,23,4)+0; if(r==53||r==272) printf "    %-3d %s\n", r, substr($0,18,3)}' \
        "$D/protein_clean.pdb"

    # point config/scripts at this folder; move off the dead gpus partition
    sed -i "s#^WORKDIR=.*#WORKDIR=\"$D\"#" "$D/config.sh"
    sed -i "s#^cd /LUSTRE.*#cd \"\$SLURM_SUBMIT_DIR\"#" "$D/06_minimize.sh" "$D/07_equilibrate.sh"
    for s in 06_minimize.sh 07_equilibrate.sh; do
        sed -i -e "s/^#SBATCH --partition=.*/#SBATCH --partition=$PARTITION/" \
               -e "/^#SBATCH --gres=gpu/d" \
               -e "s/^#SBATCH --cpus-per-task=.*/#SBATCH --cpus-per-task=$NCPU/" "$D/$s"
        grep -q "partition=$PARTITION" "$D/$s" || die "could not set the partition in $D/$s"
    done
    grep -q 'cd "\$SLURM_SUBMIT_DIR"' "$D/06_minimize.sh" || die "could not patch cd in $D/06_minimize.sh"
    echo "  ✓ config.sh, cd paths and SLURM headers patched ($PARTITION, $NCPU cores)"

    echo "  building (tleap + parmed) ..."
    ( cd "$D" && bash 05_build_system.sh > build.log 2>&1 )
    if [[ ! -f "$D/gromacs/complex.gro" ]]; then
        echo "  build FAILED. Tail of build.log:"; tail -40 "$D/build.log" | sed 's/^/    /'
        grep -i "error\|fatal" "$D"/*leap*.log 2>/dev/null | head -10 | sed 's/^/    /'
        die "$M did not produce gromacs/complex.gro"
    fi
  else
    echo "  residue identity (existing build):"
    awk '$1=="ATOM" && substr($0,13,4)~/ CA /{r=substr($0,23,4)+0; if(r==53||r==272) printf "    %-3d %s\n", r, substr($0,18,3)}' \
        "$D/protein_clean.pdb"
  fi

    # ── the mutation-specific gate: protein charge must be an exact integer ──
    # ── verify the TOPOLOGY, not tleap's log text ────────────────────────
    # Two log-parsing attempts failed here (an awk field index that hit the
    # word "is", then a latch that read past PROT into COMPLEX). Read the
    # prmtop instead: it is what actually gets simulated.
    EXP="${EXPECTED_CHARGE[$M]}"
    echo "  topology check (complex.prmtop):"
    PRMOUT=$(python3 "$HERE/prmtop_charge.py" "$D/complex.prmtop") || die "could not read $D/complex.prmtop"
    sed '$d' <<< "$PRMOUT"
    RES=$(tail -1 <<< "$PRMOUT")
    python3 - "$RES" "$EXP" << 'PYQ'
import sys, re
res, exp = sys.argv[1], float(sys.argv[2])
kv = dict(re.findall(r"(\w+)=(-?[\d.]+)", res))
if not kv: print("    [!] could not read the prmtop summary"); sys.exit(1)
ok = True
pq, pa = float(kv["protein"]), float(kv["protein_atoms"])
lq, la = float(kv["lig"]), float(kv["lig_atoms"])
fq, tq = float(kv["fad"]), float(kv["total"])
def chk(label, got, want, tol, unit="e"):
    global ok
    good = abs(got - want) <= tol
    print(f"    {'OK ' if good else '*** '}{label:<26}{got:>12.6f}  (expect {want:g} +/- {tol:g} {unit})")
    ok = ok and good
chk("protein charge", pq, exp, 1e-5)
chk("ligand charge", lq, -1.999001, 1e-4)
chk("FAD charge", fq, -2.002999, 1e-4)
chk("system total charge", tq, 0.0, 0.01)
chk("ligand atom count", la, 16, 0, "atoms")
if not ok:
    print("    the mutation or the parameters are wrong — do not proceed"); sys.exit(1)
print(f"    protein has {int(pa)} atoms")
PYQ
    [[ $? -eq 0 ]] || die "$M failed the topology check"

    NATOMS=$(sed -n '2p' "$D/gromacs/complex.gro" | tr -d ' ')
    NWAT=$(grep -oP '\d+(?= solvent molecules will remain)' "$D/build.log" 2>/dev/null | tail -1)
    PCHG=$(sed -n 's/.*protein=\(-\?[0-9.]*\).*/\1/p' <<< "$RES")
    TCHG=$(sed -n 's/.*total=\(-\?[0-9.]*\).*/\1/p' <<< "$RES")
    NSOLP=$(sed -n 's/.*protein_atoms=\([0-9]*\).*/\1/p' <<< "$RES")
    echo "  system: $NATOMS atoms, ${NWAT:-?} waters"

    # ── Solute index group by MERGING make_ndx groups ────────────────────
    ( cd "$D/gromacs"
      load_gromacs 2>/dev/null || true
      G="${GMX:-gmx_mpi}"; command -v "$G" >/dev/null || G=gmx
      printf 'q\n' | "$G" make_ndx -f complex.gro -o index_base.ndx > mk1.log 2>&1
      NG=$(grep -c '^\[' index_base.ndx)
      # group numbers for Protein / FAD / LIG, read from make_ndx's own listing
      PN=$(awk '/^ *[0-9]+ Protein +:/{print $1; exit}' mk1.log)
      FN=$(awk '/^ *[0-9]+ FAD +:/{print $1; exit}' mk1.log)
      LN=$(awk '/^ *[0-9]+ LIG +:/{print $1; exit}' mk1.log)
      echo "    merging groups: Protein=$PN FAD=$FN LIG=$LN  (index_base has $NG groups)"
      printf "%s | %s | %s\nname %s Solute\nq\n" "$PN" "$FN" "$LN" "$NG" \
        | "$G" make_ndx -f complex.gro -n index_base.ndx -o index.ndx > mk2.log 2>&1
    )
    NSOL=$(awk '/^\[ Solute/{f=1;next} /^\[/{f=0} f' "$D/gromacs/index.ndx" | tr -s ' \n' '\n' | grep -c '[0-9]')
    MN=$(awk '/^\[ Solute/{f=1;next} /^\[/{f=0} f' "$D/gromacs/index.ndx" | tr -s ' \n' '\n' | grep '[0-9]' | sort -n | head -1)
    MX=$(awk '/^\[ Solute/{f=1;next} /^\[/{f=0} f' "$D/gromacs/index.ndx" | tr -s ' \n' '\n' | grep '[0-9]' | sort -n | tail -1)
    echo "  Solute group: $NSOL atoms  (min $MN, max $MX)"
    if [[ -z "$NSOL" || "$NSOL" -lt 1 ]]; then die "$M Solute group is empty"; fi
    if (( MX - MN + 1 != NSOL )); then
        die "$M Solute group is NOT contiguous — it has picked up ions or water"
    fi
    echo "    ✓ contiguous: exactly the first $NSOL atoms (protein + FAD + LIG)"

    # ── trajectory size: write only the solute ───────────────────────────
    if grep -q '^compressed-x-grps' "$D/08_production_100ns.sh"; then
        echo "  compressed-x-grps already present"
    else
        sed -i 's|^\(nstxout-compressed  = 5000.*\)$|\1\ncompressed-x-grps   = Solute   ; solute only|' \
            "$D/08_production_100ns.sh"
        grep -q '^compressed-x-grps' "$D/08_production_100ns.sh" \
            && echo "  ✓ compressed-x-grps = Solute added" \
            || echo "  [!] could not patch the mdp — add 'compressed-x-grps = Solute' by hand"
    fi

    # ── is the ligand still in the site? ─────────────────────────────────
    python3 - "$D/complex_solvated.pdb" << 'PYGEO'
import sys, math
lig={}; fad={}
for ln in open(sys.argv[1]):
    if not ln.startswith(("ATOM","HETATM")): continue
    rn=ln[17:20].strip(); an=ln[12:16].strip()
    try: c=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
    except ValueError: continue
    if rn=="LIG": lig[an]=c
    elif rn=="FAD": fad[an]=c
d=lambda a,b: math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))
print(f"  ligand in the built system: {len(lig)} atoms")
print(f"    reacting C -> FAD N5 : {d(lig['C1'],fad['N5']):.2f} A   (WT built: 3.64)")
print(f"    P          -> FAD N5 : {d(lig['P1'],fad['N5']):.2f} A   (WT built: 6.57)")
PYGEO

    SUMMARY+=$'\n'"  $M  $NATOMS atoms | total ${TCHG:-?} | protein ${PCHG:-?} (${NSOLP:-?} atoms) | Solute ${NSOL:-?}"
done

echo; echo "$HR"; echo " SUMMARY"; echo "$HR"
echo "  WT (go19_v2): 72826 atoms, charge -0.002000, protein -4.000000, Solute 5812"
echo "$SUMMARY"
echo
echo "  Expected: R53A and R272L protein charge -5 (one Arg removed),"
echo "  R53A_R272L -6. Atom counts: ARG(24 atoms) -> ALA(10) loses 14,"
echo "  ARG -> LEU(19) loses 5, so roughly 5698 / 5707 / 5693 protein atoms."
echo "  Water counts differ slightly between systems because the solute"
echo "  bounding box changes — that is expected and correct."
echo
echo "$HR"
echo " Nothing submitted. Send me cp5_build.log before minimising."
echo "$HR"
