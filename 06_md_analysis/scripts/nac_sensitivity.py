#!/usr/bin/env python3
"""
R1 comment 1 - sensitivity of the near-attack-conformation (NAC) fraction to the
distance/angle thresholds.  Recomputes geometry directly from the fitted trajectories.

Usage (analysis PC, env mdana):
    python nac_sensitivity.py ~/GO19_MD/all
It searches <root> recursively for fit.xtc with a frame0.pdb in the same folder,
infers system/replica from the two parent folder names, and uses 20-100 ns.

Atoms: glyphosate C1 (reacting CH2) with H3/H4; FAD N5.
Per frame: d(C1-N5); for each H: d(H-N5) and angle C1-H...N5.
"closest H" = the H with the shorter H...N5 distance in that frame.
NAC (closest)  : d(C1-N5) <= D  and  angle(closest H) >= A
NAC (either H) : d(C1-N5) <= D  and  max(angle H3, angle H4) >= A
"""
import sys, glob, os, itertools
import numpy as np, MDAnalysis as mda

root = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '~/GO19_MD/all')
T0, T1 = 20000.0, 100000.0          # ps
DS = [3.3, 3.5, 3.7, 4.0]
AS = [120, 130, 140, 150]

def ang(a, b, c):                    # angle at b (degrees)
    v1, v2 = a - b, c - b
    cos = (v1 * v2).sum(-1) / (np.linalg.norm(v1, axis=-1) * np.linalg.norm(v2, axis=-1))
    return np.degrees(np.arccos(np.clip(cos, -1, 1)))

runs = sorted(glob.glob(os.path.join(root, '**', 'fit.xtc'), recursive=True))
if not runs: sys.exit(f'no fit.xtc under {root}')
rows = []
for xtc in runs:
    d = os.path.dirname(xtc); top = os.path.join(d, 'frame0.pdb')
    if not os.path.exists(top): print('skip (no frame0.pdb):', d); continue
    rep = os.path.basename(d); sysname = os.path.basename(os.path.dirname(d))
    u = mda.Universe(top, xtc)
    sel = lambda s: u.select_atoms(s)
    C1, H3, H4 = sel('resname LIG and name C1'), sel('resname LIG and name H3'), sel('resname LIG and name H4')
    N5 = sel('resname FAD and name N5')
    if min(len(C1), len(H3), len(H4), len(N5)) != 1:
        print('atom selection problem in', d, len(C1), len(H3), len(H4), len(N5)); continue
    dc, a3, a4, dh3, dh4 = [], [], [], [], []
    for ts in u.trajectory:
        if ts.time < T0 or ts.time > T1: continue
        c, h3, h4, n = C1.positions[0], H3.positions[0], H4.positions[0], N5.positions[0]
        dc.append(np.linalg.norm(c - n)); dh3.append(np.linalg.norm(h3 - n)); dh4.append(np.linalg.norm(h4 - n))
        a3.append(ang(c, h3, n)); a4.append(ang(c, h4, n))
    dc, a3, a4, dh3, dh4 = map(np.array, (dc, a3, a4, dh3, dh4))
    a_close = np.where(dh3 <= dh4, a3, a4); a_any = np.maximum(a3, a4)
    r = dict(system=sysname, rep=rep, frames=len(dc), dC1N5=dc.mean(), ang_close=a_close.mean())
    for D, A in itertools.product(DS, AS):
        r[f'close_{D}_{A}'] = 100 * np.mean((dc <= D) & (a_close >= A))
        r[f'any_{D}_{A}'] = 100 * np.mean((dc <= D) & (a_any >= A))
    rows.append(r)
    print(f'{sysname:>12s} {rep:>6s} frames={len(dc):5d}  C1-N5={dc.mean():.2f} A  angle(closest H)={a_close.mean():.1f} deg  '
          f'NAC(3.5,130,closest)={r["close_3.5_130"]:.1f}%')

import csv
keys = list(rows[0].keys())
with open('nac_sensitivity_per_rep.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

systems = sorted({r['system'] for r in rows})
for mode in ['close', 'any']:
    print(f'\n=== % NAC frames, mean +/- SD (n-1) over replicas  [{mode} H]  rows: C1-N5 <= D (A); cols: C-H...N5 >= A (deg)')
    for s in systems:
        R = [r for r in rows if r['system'] == s]
        print(f'\n{s}  (n={len(R)})\n  D\\A ' + ''.join(f'{A:>14d}' for A in AS))
        for D in DS:
            cells = []
            for A in AS:
                v = np.array([r[f'{mode}_{D}_{A}'] for r in R])
                cells.append(f'{v.mean():6.1f} +/- {v.std(ddof=1) if len(v) > 1 else 0:4.1f}')
            print(f'  {D:<4} ' + ''.join(f'{c:>14s}' for c in cells))
print('\nper-replica values written to nac_sensitivity_per_rep.csv')
