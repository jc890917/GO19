#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
#  05_build_system.sh  —  Build the full solvated system with tleap,
#                          then convert to GROMACS format
#
#  Combines: protein (ff19SB) + FAD (GAFF2) + ligand (GAFF2)
#  Solvates with OPC water, neutralizes, adds 0.15 M NaCl
#  Converts AMBER topology to GROMACS using parmed
#  Sets up position restraint files for equilibration
#
#  Input  : protein_clean.pdb, FAD.mol2, FAD.frcmod, LIG.mol2, LIG.frcmod
#  Output : complex.prmtop, complex.inpcrd   (AMBER -- kept as checkpoint)
#           gromacs/complex.gro              (GROMACS coordinates)
#           gromacs/complex.top              (GROMACS topology, patched for posres)
#           gromacs/posre_*.itp              (position restraint files, k=1000)
#           posres_setup.py                  (reusable posres generator)
# ══════════════════════════════════════════════════════════════════════════════

set -euo pipefail
source ./config.sh
load_amber

echo ""
echo "══════════════════════════════════════════════════"
echo "  Step 05 — Build solvated system with tleap"
echo "══════════════════════════════════════════════════"

# ── 1. Check all required files ───────────────────────────────────────────────
echo "  Checking required input files..."
MISSING=0
for f in "$PROTEIN_CLEAN" "$FAD_MOL2" "$FAD_FRCMOD" "$LIGAND_MOL2" "$LIGAND_FRCMOD"; do
    if [[ ! -f "$f" ]]; then
        echo "  [ERROR] Missing: $f"
        MISSING=1
    else
        echo "  ✓ Found: $f"
    fi
done
[[ $MISSING -eq 1 ]] && { echo "  Run Steps 02-04 first."; exit 1; }

# ── 2. Preliminary tleap to count water molecules ────────────────────────────
#  addIonsRand with count=0 only neutralizes -- it does NOT add salt.
#  We need the actual water count to calculate how many NaCl pairs give 0.15 M.
#  Formula: N_pairs = round(0.15 * N_waters / 55.5)
echo ""
echo "  Pass 1 -- counting water molecules for ion calculation..."

cat > tleap_pre.in << EOF
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.opc
loadamberparams ${LIGAND_FRCMOD}
${LIGAND_NAME} = loadmol2 ${LIGAND_MOL2}
loadamberparams ${FAD_FRCMOD}
${FAD_NAME} = loadmol2 ${FAD_MOL2}
PROT = loadpdb ${PROTEIN_CLEAN}
COMPLEX = combine { PROT ${FAD_NAME} ${LIGAND_NAME} }
solvateBox COMPLEX ${WATER_BOX} ${BOX_BUFFER}
quit
EOF

tleap -f tleap_pre.in 2>&1 | tee tleap_pre.log

# Extract water count from tleap output
NWAT=$(grep -i "added.*residues\|solvent molecules\|mols" tleap_pre.log \
       | grep -oE '[0-9]+' | tail -1 || true)

if [[ -z "$NWAT" ]] || [[ "$NWAT" -eq 0 ]]; then
    NWAT=$(grep -c "^WAT\b\|^WATP\b" tleap_pre.log 2>/dev/null || true)
fi

if [[ -z "$NWAT" ]] || [[ "$NWAT" -eq 0 ]]; then
    echo ""
    echo "  [ERROR] Could not auto-detect water count from tleap_pre.log."
    echo "  Check tleap_pre.log for the number of water molecules added,"
    echo "  then set NIONS manually and rerun."
    exit 1
else
    echo "  Waters in box : $NWAT"
    NIONS=$(python3 -c "import math; print(round(0.15 * ${NWAT} / 55.5))")
    echo "  NaCl pairs for 0.15 M : $NIONS"
fi

# ── 3. Write final tleap input script ────────────────────────────────────────
echo ""
echo "  Pass 2 -- writing tleap.in with correct ion counts..."

cat > tleap.in << EOF
# ============================================================
#  tleap.in -- System builder
#  Protein: ff19SB | Small molecules: GAFF2 | Water: OPC
# ============================================================

source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.opc

loadamberparams ${LIGAND_FRCMOD}
${LIGAND_NAME} = loadmol2 ${LIGAND_MOL2}
check ${LIGAND_NAME}

loadamberparams ${FAD_FRCMOD}
${FAD_NAME} = loadmol2 ${FAD_MOL2}
check ${FAD_NAME}

PROT = loadpdb ${PROTEIN_CLEAN}
check PROT

COMPLEX = combine { PROT ${FAD_NAME} ${LIGAND_NAME} }
check COMPLEX
charge COMPLEX

solvateBox COMPLEX ${WATER_BOX} ${BOX_BUFFER}

# Step 1: neutralize (Na+ only, count=0 = auto)
# Step 2: add ${NIONS} NaCl pairs for 0.15 M
#         = round(0.15 * ${NWAT} waters / 55.5)
addIonsRand COMPLEX Na+ 0
addIonsRand COMPLEX Na+ ${NIONS} Cl- ${NIONS}

check COMPLEX
charge COMPLEX
saveamberparm COMPLEX ${COMPLEX_PRMTOP} ${COMPLEX_INPCRD}
savepdb COMPLEX complex_solvated.pdb

quit
EOF

echo "  ✓ tleap.in written."

# ── 4. Run tleap ─────────────────────────────────────────────────────────────
echo ""
echo "  Running tleap..."
tleap -f tleap.in 2>&1 | tee tleap.log

# ── 5. Check outputs ─────────────────────────────────────────────────────────
echo ""
echo "  Checking tleap output..."

if [[ ! -f "$COMPLEX_PRMTOP" ]] || [[ ! -f "$COMPLEX_INPCRD" ]]; then
    echo ""
    echo "  [ERROR] tleap did not produce topology/coordinate files."
    echo "  Review tleap.log for errors."
    exit 1
fi

echo "  ✓ $COMPLEX_PRMTOP"
echo "  ✓ $COMPLEX_INPCRD"
echo "  ✓ complex_solvated.pdb"

# ── 6. Print system summary ───────────────────────────────────────────────────
echo ""
echo "  System summary from tleap.log:"
grep -E "Total residues|Total atoms|Total charge|System charge|solute|solvent|Na\+|Cl\-" tleap.log | head -20 || true

# ── 7. Verify no serious warnings ────────────────────────────────────────────
echo ""
if grep -qi "WARNING\|ERROR\|Could not" tleap.log; then
    echo "  [!] tleap reported warnings/errors -- review tleap.log before proceeding."
    grep -i "WARNING\|ERROR\|Could not" tleap.log | head -20
else
    echo "  ✓ No critical warnings in tleap.log."
fi

# ── 8. Convert AMBER topology to GROMACS format ──────────────────────────────
echo ""
echo "══════════════════════════════════════════════════"
echo "  Converting AMBER topology -> GROMACS (parmed)"
echo "══════════════════════════════════════════════════"

mkdir -p gromacs

python3 << 'PYEOF'
import sys
try:
    import parmed as pmd
except ImportError:
    print("  [ERROR] parmed not found. It should be available with AmberTools.")
    sys.exit(1)

print("  Loading AMBER topology + coordinates...")
amber = pmd.load_file('complex.prmtop', 'complex.inpcrd')
print("  Writing gromacs/complex.top ...")
amber.save('gromacs/complex.top', overwrite=True)
print("  Writing gromacs/complex.gro ...")
amber.save('gromacs/complex.gro', overwrite=True)
print("  Done.")
PYEOF

if [[ ! -f "gromacs/complex.top" ]] || [[ ! -f "gromacs/complex.gro" ]]; then
    echo "  [ERROR] GROMACS conversion failed."
    exit 1
fi

echo "  ✓ gromacs/complex.gro"
echo "  ✓ gromacs/complex.top"

# ── 9. Write posres_setup.py to project directory ────────────────────────────
#  This script reads the GRO + TOP files, generates per-molecule position
#  restraint ITP files, and patches the topology (first run only).
#  Rerun with different force constants for each equilibration stage.
echo ""
echo "  Writing posres_setup.py ..."

cat > posres_setup.py << 'POSRES_PY'
#!/usr/bin/env python3
"""
posres_setup.py -- Generate GROMACS position restraint files + patch topology.

Usage:
    python3 posres_setup.py <force_constant>

    force_constant : kJ/mol/nm^2  (e.g. 1000, 100, 10)

Reads  : gromacs/complex.gro, gromacs/complex.top
Writes : gromacs/posre_<MOL>.itp for each non-solvent molecule type
         Patches gromacs/complex.top with #ifdef POSRES blocks (first run only)
         Subsequent runs only update the posre force constants.

Stages:
    k=1000  NVT + NPT stage 1
    k=100   NPT stage 2
    k=10    NPT stage 3
    (no call needed for unrestrained runs)
"""
import re, sys, os

GRO = 'gromacs/complex.gro'
TOP = 'gromacs/complex.top'
FC  = float(sys.argv[1]) if len(sys.argv) > 1 else 1000.0

SOLVENT = {
    'WAT','SOL','HOH','TIP3','TIP4','SPC',
    'NA','CL','Na+','Cl-','SOD','CLA','ION'
}

# ── Parse topology ────────────────────────────────────────────────────────────
with open(TOP) as f:
    top_text = f.read()

def mol_atom_counts(text):
    counts = {}
    for block in re.split(r'(?=\[\s*moleculetype\s*\])', text):
        if not re.match(r'\s*\[\s*moleculetype\s*\]', block):
            continue
        mm = re.search(r'\[\s*moleculetype\s*\]\s*\n(?:\s*;[^\n]*\n)*\s*(\S+)', block)
        if not mm:
            continue
        name = mm.group(1)
        am = re.search(r'\[\s*atoms\s*\](.*?)(?=\[|\Z)', block, re.DOTALL)
        if am:
            n = sum(1 for ln in am.group(1).split('\n')
                    if ln.strip() and not ln.strip().startswith(';'))
            counts[name] = n
    return counts

def mol_list(text):
    ms = re.search(r'\[\s*molecules\s*\](.*?)$', text, re.DOTALL)
    if not ms:
        return []
    out = []
    for ln in ms.group(1).split('\n'):
        ln = ln.strip()
        if not ln or ln.startswith(';'):
            continue
        p = ln.split()
        if len(p) >= 2:
            try:
                out.append((p[0], int(p[1])))
            except ValueError:
                pass
    return out

atom_counts = mol_atom_counts(top_text)
molecules   = mol_list(top_text)

print(f"Molecule types     : {list(atom_counts.keys())}")
print(f"System composition : {molecules}")

# ── Read GRO atoms ────────────────────────────────────────────────────────────
with open(GRO) as f:
    f.readline()
    n_total = int(f.readline().strip())
    gro_atoms = []
    for _ in range(n_total):
        ln = f.readline()
        gro_atoms.append((ln[5:10].strip(), ln[10:15].strip()))

# ── Generate posre files ──────────────────────────────────────────────────────
cursor = 0
generated = []

for mol_name, n_copies in molecules:
    n_per = atom_counts.get(mol_name, 0)
    if n_per == 0:
        print(f"  Warning: no atom count for '{mol_name}'")
        continue

    if mol_name in SOLVENT:
        cursor += n_per * n_copies
        continue

    mol_atoms = gro_atoms[cursor : cursor + n_per]
    out = f'gromacs/posre_{mol_name}.itp'
    n_heavy = 0

    with open(out, 'w') as f:
        f.write(f'; Position restraints for {mol_name}\n')
        f.write(f'; FC = {FC:.1f} kJ/mol/nm^2\n')
        f.write('[ position_restraints ]\n')
        f.write('; atom  type    fx        fy        fz\n')
        for idx, (rname, aname) in enumerate(mol_atoms, 1):
            is_H = (aname.startswith('H') or aname.startswith('h') or
                    (aname and aname[0].isdigit() and 'H' in aname.upper()))
            if not is_H:
                f.write(f'{idx:5d}  1  {FC:8.1f}  {FC:8.1f}  {FC:8.1f}\n')
                n_heavy += 1

    print(f"  Generated {out}  ({n_heavy} heavy atoms, FC={FC:.1f})")
    generated.append(mol_name)
    cursor += n_per * n_copies

# ── Patch topology (first run only) ──────────────────────────────────────────
if '#ifdef POSRES' in top_text:
    print('\n  Topology already patched. Posre files updated with new FC.')
    sys.exit(0)

with open(TOP) as f:
    lines = f.readlines()

new_lines    = []
in_mol       = False
current_mol  = None

for line in lines:
    s = line.strip()

    if re.match(r'^\[\s*moleculetype\s*\]', s):
        in_mol = True; current_mol = None

    elif in_mol and current_mol is None and s and not s.startswith(';') and not s.startswith('['):
        p = s.split()
        if p:
            current_mol = p[0]

    if in_mol and current_mol and current_mol in generated:
        if re.match(r'^\[\s*(moleculetype|system)\s*\]', s):
            pf = f'posre_{current_mol}.itp'
            new_lines += ['\n',
                          f'; Position restraints for {current_mol}\n',
                          '#ifdef POSRES\n',
                          f'#include "{pf}"\n',
                          '#endif\n', '\n']
            in_mol = False; current_mol = None
    elif in_mol and current_mol and re.match(r'^\[\s*(moleculetype|system)\s*\]', s):
        in_mol = False; current_mol = None

    new_lines.append(line)

with open(TOP, 'w') as f:
    f.writelines(new_lines)

print(f'\n  Topology patched for: {generated}')
print("  Use 'define = -DPOSRES' in MDP files to activate restraints.")
print(f'\nDone. FC = {FC:.1f} kJ/mol/nm^2')
POSRES_PY

chmod +x posres_setup.py
echo "  ✓ posres_setup.py written."

# ── 10. Run posres_setup.py (k=1000, initial setup) ──────────────────────────
echo ""
echo "  Setting up position restraints (k=1000 kJ/mol/nm^2)..."
python3 posres_setup.py 1000

# ── 11. Generate GROMACS index file ──────────────────────────────────────────
echo ""
echo "  Generating GROMACS index file..."
load_gromacs
echo "q" | $GMX make_ndx -f gromacs/complex.gro -o gromacs/index.ndx 2>&1 | tee gromacs/make_ndx.log
echo "  ✓ gromacs/index.ndx"

# ── 12. Summary ──────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════════"
echo "  Step 05 complete."
echo ""
echo "  AMBER files (checkpoint):"
echo "    $COMPLEX_PRMTOP"
echo "    $COMPLEX_INPCRD"
echo "    complex_solvated.pdb"
echo ""
echo "  GROMACS files (for simulation):"
echo "    gromacs/complex.gro    <- coordinates"
echo "    gromacs/complex.top    <- topology (patched for POSRES)"
echo "    gromacs/index.ndx      <- index file"
echo "    gromacs/posre_*.itp    <- restraints (k=1000 kJ/mol/nm^2)"
echo "    posres_setup.py        <- rerun with new k to update restraints"
echo ""
echo "  NEXT: sbatch 06_minimize.sh"
echo "══════════════════════════════════════════════════"