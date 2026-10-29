#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
#  config.sh  —  Master configuration for the MD preparation pipeline
#  All other scripts source this file. Edit ONLY this file to change settings.
# ══════════════════════════════════════════════════════════════════════════════

# ── Module system initialization (required for SLURM jobs) ────────────────────
source /etc/profile.d/modules.sh

# ── Working directory (Yotla) ─────────────────────────────────────────────────
WORKDIR="/LUSTRE/home/lancad/2026/l.100073/alumnos/ian/go19_v2"

# ── Input files ───────────────────────────────────────────────────────────────
PROTEIN_PDB="go19.pdb"      # Original PDB containing protein + FAD
LIGAND_SDF="selected_pose.sdf"     # Extracted best pose from GNINA (Step 01)

# ── Residue / molecule names ──────────────────────────────────────────────────
LIGAND_NAME="LIG"                  # 3-letter residue name for the docked ligand
FAD_NAME="FAD"                     # Residue name for FAD in the PDB

# ── Charges  ──────────────────────────────────────────────────────────────────
# Glyphosate at pH 7.4:
#   Phosphonic –OH ×2  → deprotonated → -2
#   Carboxylic –COOH   → deprotonated → -1
#   Amine –NH          → protonated   → +1
#   Net charge = -2
LIGAND_CHARGE=-2                   # Glyphosate at physiological pH 7.4
FAD_CHARGE=-2                      # FAD (oxidized, FADH0 → FAD²⁻) net charge = -2
                                   # Change to -1 if your FAD is semiquinone (FADH•)
                                   # Change to  0 if your FAD is fully reduced (FADH₂)

# ── Force field choices ───────────────────────────────────────────────────────
PROTEIN_FF="ff19SB"                # Protein force field
LIGAND_FF="gaff2"                  # Small molecule force field (ligand + FAD)
CHARGE_METHOD="bcc"                # AM1-BCC charges
WATER_MODEL="opc"                  # OPC water (recommended with ff19SB)
WATER_BOX="OPCBOX"                 # tleap box name for OPC

# ── Solvation ─────────────────────────────────────────────────────────────────
BOX_BUFFER=12.0                    # Water box buffer in Angstroms
SALT_CONC=0.15                     # NaCl concentration in mol/L

# ── Production MD ─────────────────────────────────────────────────────────────
PRODUCTION_NS=50                  # !! SET THIS: production length in nanoseconds
TEMP=300                           # Temperature in Kelvin
PRESSURE=1.0                       # Pressure in bar

# ── SLURM settings (Yotla) ────────────────────────────────────────────────────
SLURM_PARTITION="gpus"             # GPU partition on Yotla (nodes: ngk1-24)
SLURM_GPUS=1                       # Number of GPUs per job
SLURM_CPUS=10                       # Half of 20 CPUs per ngk node
SLURM_MEM="32G"                    # Safe headroom on 64 GB nodes
SLURM_EMAIL="jccarrillo@uach.mx"       # !! SET THIS: your email for job notifications

# ── Derived values (do not edit) ─────────────────────────────────────────────
PRODUCTION_STEPS=$(echo "$PRODUCTION_NS * 500000" | bc)   # steps at dt=0.002 ps
PROTEIN_CLEAN="protein_clean.pdb"
FAD_CLEAN="FAD_clean.pdb"
LIGAND_MOL2="${LIGAND_NAME}.mol2"
LIGAND_FRCMOD="${LIGAND_NAME}.frcmod"
FAD_MOL2="FAD_bcc.mol2"
FAD_FRCMOD="${FAD_NAME}.frcmod"
COMPLEX_PRMTOP="complex.prmtop"
COMPLEX_INPCRD="complex.inpcrd"
 
# ── GROMACS file paths ────────────────────────────────────────────────────────
GMX_GRO="gromacs/complex.gro"      # GROMACS coordinate file
GMX_TOP="gromacs/complex.top"      # GROMACS topology file
GMX_NDX="gromacs/index.ndx"        # GROMACS index file
GMX="gmx_mpi"                      # GROMACS executable name on Yotla

# ── Module loads ──────────────────────────────────────────────────────────────
load_amber() {
    module load amber/22
    module load amber/ambertools19
    module load intel-oneapi/compilers-2022.2.1
    module load anaconda/2021.04
    module load gromacs/2024.3-gpu
}

load_openbabel() {
    module load openbabel/2.2.3
}

load_anaconda() {
    module load anaconda/2021.04
}

load_gromacs() {
    module load gromacs/2024.3-gpu
}
# ── Tool paths (Yotla) ────────────────────────────────────────────────────────
REDUCE="/LUSTRE/yoltla/programas/amber/ambertools19_INTEL-parallel_studio_xe_2017_update4_IMPI/amber18/bin/reduce"
