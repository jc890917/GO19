#!/usr/bin/env bash
# Run on the Ubuntu analysis PC.
#   bash collect_from_analysis_pc.sh            -> ~/GO19_release/
# Produces:
#   ~/GO19_release/scripts/              -> copy into 06_md_analysis/scripts/ (GitHub)
#   ~/GO19_release/mmgbsa_topologies/    -> copy into 05_md_setup/ (GitHub; check sizes)
#   ~/GO19_release/GO19_MD_trajectories.tar.gz  -> upload to Zenodo (NOT GitHub)
set -uo pipefail
SRC="${SRC:-$HOME/GO19_MD}"
OUT="${OUT:-$HOME/GO19_release}"
mkdir -p "$OUT/scripts" "$OUT/mmgbsa_topologies" "$OUT/traj"

echo "== analysis scripts"
find "$SRC" -maxdepth 3 -type f \( -name "*.py" -o -name "*.sh" -o -name "*.in" -o -name "*.mdin" \) \
     -not -path "*/results/*" -not -path "*/mdpocket/*/*" -print -exec cp -n {} "$OUT/scripts/" \;

echo "== stripped topologies used by MMPBSA.py / ProLIF"
cp -v "$SRC"/topol/*.prmtop "$OUT/mmgbsa_topologies/" 2>/dev/null

echo "== trajectories (fit.xtc + frame0.pdb per system/replicate)"
find "$SRC/all" -name fit.xtc | while read -r x; do
  d=$(dirname "$x"); rel=${d#"$SRC/all/"}; mkdir -p "$OUT/traj/$rel"
  cp "$x" "$OUT/traj/$rel/"; [ -f "$d/frame0.pdb" ] && cp "$d/frame0.pdb" "$OUT/traj/$rel/"
  echo "   $rel"
done
cat > "$OUT/traj/README.txt" <<'TXT'
GO19 MD trajectories: wild type (go19_v2), R53A, R272L, R53A/R272L; 3 x 100 ns each.
Solute only (protein + FAD + glyphosate), PBC-corrected (-pbc mol -center, -pbc nojump),
fitted on the protein backbone. frame0.pdb is the matching topology/first frame.
Simulation inputs and analysis scripts: see the GitHub/Zenodo software record.
TXT
tar czf "$OUT/GO19_MD_trajectories.tar.gz" -C "$OUT" traj
du -sh "$OUT"/* ; echo "Check the list above: delete any scripts that were superseded before copying them to GitHub."
