#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
#  launch_all.sh — Submit the full mutant campaign in one go
#
#  Run on the Yoltla login node AFTER setup_mutants.sh:
#      bash launch_all.sh
#
#  Per mutant:  minimization → equilibration → 3 × 100 ns production
#  chained with SLURM dependencies, so every job enters the queue now and
#  starts as soon as its predecessor finishes. Uses the WT's own 06/07/08
#  scripts unchanged (same MDP settings and the same per-replica seeds).
# ══════════════════════════════════════════════════════════════════════════════
set -euo pipefail
MUTANTS=("R53A" "R272L" "R53A_R272L")
HERE="$(cd "$(dirname "$0")" && pwd)"
PARENT="/LUSTRE/home/lancad/2026/l.100073/alumnos/ian"

for M in "${MUTANTS[@]}"; do
    D="$PARENT/go19_$M"
    [[ -f "$D/gromacs/complex.gro" ]] || { echo "[ERROR] $M not built — run setup_mutants.sh first"; exit 1; }
    cd "$D"
    J_MIN=$(sbatch --parsable --job-name="${M}_min" 06_minimize.sh)
    J_EQ=$(sbatch --parsable --job-name="${M}_eq" --dependency=afterok:${J_MIN} 07_equilibrate.sh)
    for REP in 1 2 3; do
        J_P=$(sbatch --parsable --job-name="${M}_r${REP}" \
              --output="prod_rep${REP}_%j.out" --error="prod_rep${REP}_%j.err" \
              --dependency=afterok:${J_EQ} 08_production_100ns.sh "$REP")
        echo "$M rep$REP → job $J_P"
    done
    echo "$M: min=$J_MIN  equil=$J_EQ"
    echo "$M min=$J_MIN eq=$J_EQ" >> "$HERE/submitted_jobs.txt"
done
echo ""
echo "All submitted. Check with:  squeue -u \$USER"
