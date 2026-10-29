#!/usr/bin/env python3
"""Verify that position restraints were ACTUALLY active in a minimisation stage.

Usage: python3 check_posres.py <before.gro> <after.gro>

Compares solute heavy-atom positions between two .gro files.  Solute = every
residue that is not water or a monatomic ion, i.e. protein + FAD + LIG.

IMPORTANT — this script does NOT decide whether restraints were active.
At the k = 1000 kJ/mol/nm^2 used here, moving an atom 1 A costs only ~5
kJ/mol, so charged side-chain tips (Glu/Asp carboxylates) routinely travel
1-2 A even when restraints ARE active.  There is no absolute displacement
threshold that separates the two cases.

Use it for one of two comparisons instead:
  (a) same stage, patched vs unpatched topology -> steepest descent is
      deterministic, so any systematic reduction is due to the restraints;
  (b) across systems that should behave alike (WT vs mutants).

The authoritative check on whether restraints reached the run input is
  gmx dump -s min/min1.tpr | grep -icE "posres|Position Rest"
which reads what grompp actually compiled.  `define = -DPOSRES` with no
matching `#ifdef POSRES` in the topology yields zero, silently.
"""
import sys, math

SOLVENT = {"WAT", "HOH", "SOL", "TIP3", "OPC", "Na+", "Cl-", "NA", "CL", "K+"}

def read_gro(p):
    L = open(p).read().splitlines()
    n = int(L[1])
    out = []
    for ln in L[2:2+n]:
        rn = ln[5:10].strip()
        an = ln[10:15].strip()
        if rn in SOLVENT or an.startswith("H"):
            continue
        out.append((int(ln[0:5]), rn, an,
                    (float(ln[20:28])*10, float(ln[28:36])*10, float(ln[36:44])*10)))
    return out

a, b = read_gro(sys.argv[1]), read_gro(sys.argv[2])
if len(a) != len(b):
    sys.exit(f"atom-count mismatch: {len(a)} vs {len(b)} — not the same system")

d = lambda u, v: math.sqrt(sum((u[i]-v[i])**2 for i in range(3)))
disp = [(d(x[3], y[3]), f"{x[1]}{x[0]}:{x[2]}") for x, y in zip(a, b)]
disp.sort(reverse=True)

mx = disp[0][0]
mean = sum(v for v, _ in disp)/len(disp)
rmsd = math.sqrt(sum(v*v for v, _ in disp)/len(disp))

print(f"  {sys.argv[1]}  ->  {sys.argv[2]}")
print(f"  solute heavy atoms compared : {len(disp)}")
print(f"  mean displacement           : {mean:6.3f} A")
print(f"  RMSD                        : {rmsd:6.3f} A")
print(f"  max  displacement           : {mx:6.3f} A   ({disp[0][1]})")
print("  largest 8 movers            : " + ", ".join(f"{w} {v:.2f}" for v, w in disp[:8]))
tips = sum(1 for v, w in disp[:20]
           if w.split(":")[1] in ("OE1","OE2","OD1","OD2","NZ","NH1","NH2","OH"))
print(f"  of the 20 largest movers, {tips} are charged/polar sidechain tips")
print("  NO VERDICT: see the header — displacement alone cannot settle this.")
print("  Authoritative: gmx dump -s <stage>.tpr | grep -icE \"posres|Position Rest\"")
