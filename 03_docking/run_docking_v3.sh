#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════════
#  run_docking_v2.sh — GO19 redocking with correctly protonated glyphosate (CP2)
#
#  Changes vs run_docking.sh (original, used for the manuscript):
#   - Ligand: glifosato_m2.sdf  (net -2: NH2+, COO-, PO3 2-; 6 H) instead of the
#             neutral glifosato.pdbqt. SDF is read directly by GNINA, so no
#             OpenBabel step can re-protonate it.
#   - Receptor: go19_receptor_fixed.pdb = go19.pdb with the spurious FAD N1-H
#               removed (oxidised FAD; N3-H kept). Nothing else changed.
#   - Flexible residues: A:314, A:272, A:257 (R272 added: it occludes the pocket
#               in the AlphaFold model).
#   - --out_flex saves the flexible side-chain conformations of every pose.
#   - 5 independent seeds to test convergence of the top pose.
#  Unchanged: box centre on FAD N5, 20 x 20 x 20 A, CNN rescoring, GNINA 1.3.2.
#
#  Run from this folder:   bash run_docking_v2.sh
#  (If there is no usable GPU, add  --no_gpu  to the gnina command.)
# ══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

RECEPTOR="go19_receptor_fixed.pdb"
LIGANDO="glifosato_m2.sdf"
CENTER_X=-2.857; CENTER_Y=1.645; CENTER_Z=4.71      # FAD N5
SIZE_X=20; SIZE_Y=20; SIZE_Z=20
FLEXRES="A:257"
EXH=32
NMODES=20
SEEDS=(43 102 203 304 405)

gnina --version | head -1 | tee gnina_version.txt

for S in "${SEEDS[@]}"; do
  echo "── seed $S ──"
  gnina \
    --receptor "$RECEPTOR" \
    --ligand "$LIGANDO" \
    --center_x $CENTER_X --center_y $CENTER_Y --center_z $CENTER_Z \
    --size_x $SIZE_X --size_y $SIZE_Y --size_z $SIZE_Z \
    --flexres "$FLEXRES" \
    --out "poses_seed${S}.sdf" \
    --out_flex "flex_seed${S}.pdb" \
    --log "log_seed${S}.txt" \
    --exhaustiveness $EXH \
    --num_modes $NMODES \
    --seed $S \
    --cnn_scoring rescore
done
echo "Done. Send me the folder back (poses_seed*.sdf, flex_seed*.pdb, log_seed*.txt)."
