#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
#  06_minimize.sh  —  Energy minimization (2-stage, GROMACS/GPU)
#
#  Stage 1: Restrained steepest descent  — solute fixed, solvent relaxes
#  Stage 2: Unrestrained steepest descent — full system relaxes
#
#  Submit with: sbatch 06_minimize.sh
# ══════════════════════════════════════════════════════════════════════════════
#SBATCH --job-name=min
#SBATCH --partition=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=min_%j.log
#SBATCH --error=min_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jccarrillo@uach.mx

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./config.sh
load_gromacs

echo ""
echo "══════════════════════════════════════════════════"
echo "  Step 06 — Energy minimization (GROMACS)"
echo "  Job ID : $SLURM_JOB_ID"
echo "  Node   : $SLURMD_NODENAME"
echo "══════════════════════════════════════════════════"

# ── Check input files ─────────────────────────────────────────────────────────
for f in "$GMX_GRO" "$GMX_TOP" "$GMX_NDX"; do
    [[ ! -f "$f" ]] && { echo "[ERROR] Missing: $f — run Step 05 first."; exit 1; }
done

mkdir -p min

# ══════════════════════════════════════════════════════════════════════════════
#  Stage 1 — Restrained minimization (solute held, solvent/ions relax)
# ══════════════════════════════════════════════════════════════════════════════
cat > min/min1.mdp << 'EOF'
; Restrained steepest descent -- relax solvent while holding solute
integrator      = steep
emtol           = 1000.0        ; stop when max force < 1000 kJ/mol/nm
emstep          = 0.01
nsteps          = 50000

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
ns_type         = grid

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0
pme_order       = 4
fourierspacing  = 0.16

; Van der Waals
rvdw            = 1.0

; No bonds constraints during minimization
constraints     = none

; Periodic boundary
pbc             = xyz

; Position restraints on solute (activate #ifdef POSRES in topology)
define          = -DPOSRES
EOF

echo "  Preprocessing Stage 1 (restrained)..."
$GMX grompp \
    -f min/min1.mdp \
    -c "$GMX_GRO" \
    -r "$GMX_GRO" \
    -p "$GMX_TOP" \
    -n "$GMX_NDX" \
    -o min/min1.tpr \
    -maxwarn 5 \
    2>&1 | tee min/min1_grompp.log

echo "  Running Stage 1..."
$GMX mdrun \
    -v \
    -deffnm min/min1 \
    -s min/min1.tpr \
    -ntomp "${SLURM_CPUS_PER_TASK:-10}" \
    2>&1 | tee min/min1_mdrun.log

echo "  ✓ Stage 1 complete."

# ══════════════════════════════════════════════════════════════════════════════
#  Stage 2 — Unrestrained minimization (everything relaxes)
# ══════════════════════════════════════════════════════════════════════════════
cat > min/min2.mdp << 'EOF'
; Unrestrained steepest descent -- full system relaxes
integrator      = steep
emtol           = 500.0
emstep          = 0.01
nsteps          = 50000

cutoff-scheme   = Verlet
nstlist         = 10
ns_type         = grid

coulombtype     = PME
rcoulomb        = 1.0
pme_order       = 4
fourierspacing  = 0.16

rvdw            = 1.0
constraints     = none
pbc             = xyz

; No restraints
EOF

echo "  Preprocessing Stage 2 (unrestrained)..."
$GMX grompp \
    -f min/min2.mdp \
    -c min/min1.gro \
    -p "$GMX_TOP" \
    -n "$GMX_NDX" \
    -o min/min2.tpr \
    -maxwarn 5 \
    2>&1 | tee min/min2_grompp.log

echo "  Running Stage 2..."
$GMX mdrun \
    -v \
    -deffnm min/min2 \
    -s min/min2.tpr \
    -ntomp "${SLURM_CPUS_PER_TASK:-10}" \
    2>&1 | tee min/min2_mdrun.log

echo "  ✓ Stage 2 complete."

# ── Sanity check ─────────────────────────────────────────────────────────────
echo ""
echo "  Potential energy at end of minimization:"
$GMX energy -f min/min2.edr -o min/min2_energy.xvg << 'EOF' 2>/dev/null | grep -A2 "Potential"
Potential
EOF
echo ""
echo "  Output restart : min/min2.gro"
echo "  Review min/min1.log and min/min2.log — energy should decrease smoothly."
echo "  NEXT: sbatch 07_equilibrate.sh"
