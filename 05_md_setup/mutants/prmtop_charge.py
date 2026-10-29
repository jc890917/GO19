#!/usr/bin/env python3
"""Report per-component charges straight from an AMBER prmtop.

Reads the topology that will actually be simulated, instead of parsing tleap's
log text. Two attempts at log parsing failed on this project -- an awk field
index that landed on the word "is", then a latch that never reset and read past
'Checking PROT' into 'Checking COMPLEX' -- so check the authoritative object.

Usage:  python3 prmtop_charge.py complex.prmtop
Last line is machine readable:
        RESULT protein=<q> protein_atoms=<n> lig=<q> lig_atoms=<n> fad=<q> total=<q>
"""
import sys

def sections(path):
    out, flag, fmt, buf = {}, None, None, []
    for ln in open(path):
        if ln.startswith("%FLAG"):
            if flag: out[flag] = (fmt, buf)
            flag = ln.split()[1]; fmt = None; buf = []
        elif ln.startswith("%FORMAT"):
            fmt = ln.strip()
        elif flag:
            buf.append(ln.rstrip("\n"))
    if flag: out[flag] = (fmt, buf)
    return out

def fixed(lines, width):
    vals = []
    for ln in lines:
        for i in range(0, len(ln), width):
            tok = ln[i:i+width].strip()
            if tok: vals.append(tok)
    return vals

p = sys.argv[1]
S = sections(p)
for need in ("CHARGE", "RESIDUE_LABEL", "RESIDUE_POINTER", "POINTERS"):
    if need not in S: sys.exit(f"{p}: missing %FLAG {need}")

ptr  = [int(x) for x in fixed(S["POINTERS"][1], 8)]
natom, nres = ptr[0], ptr[11]
chg  = [float(x)/18.2223 for x in fixed(S["CHARGE"][1], 16)]
lab  = [x for x in fixed(S["RESIDUE_LABEL"][1], 4)]
rp   = [int(x) for x in fixed(S["RESIDUE_POINTER"][1], 8)]
if len(chg) != natom: sys.exit(f"{p}: CHARGE has {len(chg)} entries, NATOM says {natom}")
if len(lab) != nres:  sys.exit(f"{p}: {len(lab)} residue labels, NRES says {nres}")

NONPROT = {"FAD","LIG","WAT","HOH","Na+","Cl-","NA","CL","K+"}
groups = {}
for i, name in enumerate(lab):
    a0 = rp[i] - 1
    a1 = (rp[i+1] - 1) if i + 1 < nres else natom
    key = name if name in NONPROT else "PROTEIN"
    g = groups.setdefault(key, [0.0, 0, 0])
    g[0] += sum(chg[a0:a1]); g[1] += a1 - a0; g[2] += 1

print(f"  {'component':<12}{'charge':>12}{'atoms':>9}{'residues':>10}")
order = ["PROTEIN","FAD","LIG","Na+","Cl-","WAT"]
for k in order + [k for k in groups if k not in order]:
    if k not in groups: continue
    q, na, nr = groups[k]
    print(f"  {k:<12}{q:>12.6f}{na:>9}{nr:>10}")
tot = sum(g[0] for g in groups.values())
print(f"  {'TOTAL':<12}{tot:>12.6f}{natom:>9}{nres:>10}")
g = lambda k, i, d=0: groups[k][i] if k in groups else d
print(f"RESULT protein={g('PROTEIN',0):.6f} protein_atoms={g('PROTEIN',1)} "
      f"lig={g('LIG',0):.6f} lig_atoms={g('LIG',1)} fad={g('FAD',0):.6f} total={tot:.6f}")
