#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
#  setup_mutants.sh — Build R53A, R272L and R53A/R272L systems for GO19
#
#  Run ONCE on the Yoltla login node (from any folder):
#      bash setup_mutants.sh
#  Creates  .../ian/go19_R53A, go19_R272L, go19_R53A_R272L  next to go19/
#
#  Takes the files already prepared for the WT (same docked pose, same FAD,
#  same AM1-BCC/GAFF2 parameters, same ff19SB/OPC setup), applies each point
#  mutation to protein_clean.pdb, and builds each system with the WT's own
#  05_build_system.sh. The only difference from the WT is the mutation itself.
# ══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

# ── EDIT IF NEEDED ────────────────────────────────────────────────────────────
WT_DIR="/LUSTRE/home/lancad/2026/l.100073/alumnos/ian/go19"   # WT system folder
MUTANTS=("R53A" "R272L" "R53A_R272L")
# ──────────────────────────────────────────────────────────────────────────────

HERE="$(cd "$(dirname "$0")" && pwd)"
PARENT="$(dirname "$WT_DIR")"   # .../alumnos/ian — mutants go next to go19/ so md_scripts_reps works unchanged

# Files that have to come from the WT folder
NEEDED=(config.sh protein_clean.pdb LIG.mol2 LIG.frcmod FAD_bcc.mol2 FAD.frcmod
        05_build_system.sh 06_minimize.sh 07_equilibrate.sh 08_production_100ns.sh)
for f in "${NEEDED[@]}"; do
    [[ -f "$WT_DIR/$f" ]] || { echo "[ERROR] Missing $WT_DIR/$f"; exit 1; }
done
echo "✓ All WT input files found in $WT_DIR"

# ── Mutation function (pure awk, no extra dependencies) ──────────────────────
# mutate IN OUT RESNUM NEWRES
#   ARG→ALA : keep N CA C O CB
#   ARG→LEU : keep N CA C O CB CG, rename CD→CD1 (tleap rebuilds CD2 and all H)
mutate() {
    local IN="$1" OUT="$2" RES="$3" NEW="$4"
    local OLD
    OLD=$(awk -v r="$RES" '($1=="ATOM"||$1=="HETATM") && substr($0,23,4)+0==r {print substr($0,18,3); exit}' "$IN")
    if [[ "$OLD" != "ARG" ]]; then
        echo "[ERROR] Residue $RES in $IN is '$OLD', expected ARG. Numbering changed? Aborting."
        exit 1
    fi
    awk -v r="$RES" -v new="$NEW" '
    function trim(s){ gsub(/ /,"",s); return s }
    {
      if (($1=="ATOM"||$1=="HETATM") && substr($0,23,4)+0==r) {
        name = trim(substr($0,13,4))
        keep = 0
        if (new=="ALA" && (name=="N"||name=="CA"||name=="C"||name=="O"||name=="CB")) keep=1
        if (new=="LEU" && (name=="N"||name=="CA"||name=="C"||name=="O"||name=="CB"||name=="CG"||name=="CD")) keep=1
        if (!keep) next
        line = substr($0,1,17) new substr($0,21)
        if (new=="LEU" && name=="CD") line = substr(line,1,12) " CD1" substr(line,17)
        print line
      } else print
    }' "$IN" > "$OUT"
    echo "   ARG$RES → $NEW applied"
}

for M in "${MUTANTS[@]}"; do
    echo ""
    echo "══════════ $M ══════════"
    D="$PARENT/go19_$M"
    if [[ -d "$D/gromacs" ]]; then
        echo "   $D/gromacs already exists — skipping (delete the folder to rebuild)."
        continue
    fi
    mkdir -p "$D"
    for f in "${NEEDED[@]}"; do cp "$WT_DIR/$f" "$D/"; done

    # Mutate the protein
    cp "$D/protein_clean.pdb" "$D/protein_clean_WT.pdb"
    case "$M" in
        R53A)       mutate "$D/protein_clean_WT.pdb" "$D/protein_clean.pdb" 53  ALA ;;
        R272L)      mutate "$D/protein_clean_WT.pdb" "$D/protein_clean.pdb" 272 LEU ;;
        R53A_R272L) mutate "$D/protein_clean_WT.pdb" "$D/tmp.pdb"            53  ALA
                    mutate "$D/tmp.pdb"              "$D/protein_clean.pdb" 272 LEU
                    rm -f "$D/tmp.pdb" ;;
    esac

    # Point config and the 06/07 scripts at THIS folder (WT scripts had a
    # hard-coded 'cd' to another protein folder)
    sed -i "s#^WORKDIR=.*#WORKDIR=\"$D\"#" "$D/config.sh"
    sed -i "s#^cd /LUSTRE.*#cd \"\$SLURM_SUBMIT_DIR\"#" "$D/06_minimize.sh" "$D/07_equilibrate.sh"
    grep -q 'cd "\$SLURM_SUBMIT_DIR"' "$D/06_minimize.sh"    || { echo "[ERROR] could not patch cd in 06"; exit 1; }
    grep -q 'cd "\$SLURM_SUBMIT_DIR"' "$D/07_equilibrate.sh" || { echo "[ERROR] could not patch cd in 07"; exit 1; }


    # Build the solvated system (tleap + parmed + posres + index), as for WT
    ( cd "$D" && bash 05_build_system.sh > build.log 2>&1 ) \
        || { echo "[ERROR] 05_build_system.sh failed for $M — see $D/build.log"; exit 1; }

    NAT=$(sed -n '2p' "$D/gromacs/complex.gro" | tr -d ' ')
    CHG=$(grep -i "Total unperturbed charge" "$D/tleap.log" | tail -1 || true)
    echo "   ✓ System built: $NAT atoms  | final $CHG   (must be 0.000000)"
    echo "   Residue check:"
    grep -E " (53|272) " <(awk '$1=="ATOM" && $3=="CA" {print $4, $6, ""}' "$D/protein_clean.pdb") || true
done

echo ""
echo "══════════════════════════════════════════════════"
echo " Setup done. Review each */build.log (tleap warnings), then:"
echo "     bash launch_all.sh"
echo "══════════════════════════════════════════════════"
