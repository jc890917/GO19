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
