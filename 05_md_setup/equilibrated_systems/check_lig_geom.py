#!/usr/bin/env python3
"""Catalytic geometry + phosphonate salt-bridge check for a .gro or .pdb.

Usage: python3 check_lig_geom.py <file.gro|file.pdb> [resid1 resid2 ...]

Part 1  hydride-transfer geometry, against the docked seed404 reference:
          C1-N5 3.64 A, C1-C4X 3.91 A, P1-N5 6.57 A
Part 2  the manuscript's phosphonate claim (Arg53 + Arg272).
        Default probe residues 53 and 272; override on the command line.
        The residue NAME found at each number is printed, so a numbering
        offset between the paper and the topology is visible, not silent.
Part 3  unbiased scan: every Arg/Lys/His/Tyr/Ser/Thr donor within 5.0 A of
        any phosphonate oxygen, whatever its number.  If the paper's Arg53
        is a different residue number here, it shows up in this list.

LIG atom names (from CP4_build/LIG.mol2):
  N1 ammonium   C1 reacting CH2   C2/O1/O2 carboxylate
  C3 P1         O3/O4/O5 phosphonate
"""
import sys, math, collections

p = sys.argv[1]
probes = [int(x) for x in sys.argv[2:]] or [53, 272]

# atoms: list of (resid, resname, atomname, xyz)
atoms = []
if p.endswith(".gro"):
    L = open(p).read().splitlines()
    n = int(L[1])
    for ln in L[2:2+n]:
        ri = int(ln[0:5]); rn = ln[5:10].strip(); an = ln[10:15].strip()
        xyz = (float(ln[20:28])*10, float(ln[28:36])*10, float(ln[36:44])*10)
        atoms.append((ri, rn, an, xyz))
else:
    for ln in open(p):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        try:
            ri = int(ln[22:26])
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        atoms.append((ri, ln[17:20].strip(), ln[12:16].strip(), xyz))

lig = {a[2]: a[3] for a in atoms if a[1] == "LIG"}
fad = {a[2]: a[3] for a in atoms if a[1] == "FAD"}
d = lambda a, b: math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))
if not lig or "N5" not in fad:
    sys.exit(f"{p}: could not find LIG / FAD N5")

# ── Part 1 ────────────────────────────────────────────────────────────────
ref = {"C1-N5": 3.64, "C1-C4X": 3.91, "P1-N5": 6.57}
got = {"C1-N5":  d(lig["C1"], fad["N5"]),
       "C1-C4X": d(lig["C1"], fad["C4X"]),
       "P1-N5":  d(lig["P1"], fad["N5"])}
print(f"  {p}   ({len(lig)} LIG atoms, {len(atoms)} total)")
print("  -- hydride-transfer geometry --")
for k in ref:
    delta = got[k] - ref[k]
    flag = "" if abs(delta) < 0.5 else ("   <<< moved" if abs(delta) < 1.5 else "   <<< MOVED A LOT")
    print(f"    {k:<8}{got[k]:>7.2f} A   (docked {ref[k]:.2f}, change {delta:+.2f}){flag}")

# ── Part 2 ────────────────────────────────────────────────────────────────
PHOS = [o for o in ("O3", "O4", "O5") if o in lig]
CARB = [o for o in ("O1", "O2") if o in lig]
DONOR = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",),
         "HIS": ("ND1", "NE2"), "HIE": ("ND1", "NE2"),
         "HID": ("ND1", "NE2"), "HIP": ("ND1", "NE2"),
         "TYR": ("OH",), "SER": ("OG",), "THR": ("OG1",),
         "ASN": ("ND2",), "GLN": ("NE2",), "TRP": ("NE1",)}

byres = collections.defaultdict(list)
resname = {}
for ri, rn, an, xyz in atoms:
    byres[ri].append((an, xyz))
    resname[ri] = rn

print("  -- phosphonate anchoring at the probed residue numbers --")
print("     (docked seed404: Arg53 NH1 5.85 A; NO phosphonate partner <= 3.5 A)")
for ri in probes:
    if ri not in byres:
        print(f"    residue {ri}: NOT PRESENT in this file")
        continue
    rn = resname[ri]
    names = DONOR.get(rn)
    print(f"    residue {ri} is {rn}" + ("" if names else "   <<< no sidechain donor; mutated or numbering offset?"))
    if not names:
        continue
    crd = dict(byres[ri])
    for nm in names:
        if nm not in crd:
            print(f"      {nm:<4} absent")
            continue
        dp = min((d(crd[nm], lig[o]), o) for o in PHOS) if PHOS else (float("nan"), "-")
        dc = min((d(crd[nm], lig[o]), o) for o in CARB) if CARB else (float("nan"), "-")
        mark = "  SALT BRIDGE" if dp[0] <= 3.5 else ("  close" if dp[0] <= 4.5 else "")
        print(f"      {nm:<4} -> phosphonate {dp[1]} {dp[0]:6.2f} A{mark:<13}"
              f"   carboxylate {dc[1]} {dc[0]:6.2f} A")

# ── Part 3 ────────────────────────────────────────────────────────────
# Unbiased scan of each ligand charge group against every sidechain /
# backbone-carbonyl partner.  Docked reference (CP2, 6 catalytic seeds):
#   carboxylate : Arg314 NE 2.79-3.02, Arg314 NH2 2.84-3.01,
#                 Arg272 NH2 2.88-3.08, Tyr257 OH 3.03-3.08,
#                 Pro341 O 3.24-3.43
#   ammonium    : Pro341 O 2.80-3.21 (5 of 6 seeds)
#   phosphonate : nothing <= 3.5 A; Arg53 NH1 5.0-7.4 A
GROUPS = [("phosphonate", PHOS), ("carboxylate", CARB),
          ("ammonium", [a for a in ("N1",) if a in lig])]

# acceptors too, so backbone C=O partners (Pro341 O) are not missed
BACKBONE_O = ("O", "OXT")

print("  -- unbiased scan: partners within 5.0 A of each ligand charge group --")
for label, gatoms in GROUPS:
    if not gatoms:
        continue
    hits = []
    for ri, lst in byres.items():
        rn = resname[ri]
        names = DONOR.get(rn, ())
        for an, xyz in lst:
            if rn in ("LIG", "WAT", "HOH", "SOL", "Na+", "Cl-"):
                continue
            if an not in names and an not in BACKBONE_O:
                continue
            for o in gatoms:
                dd = d(xyz, lig[o])
                if dd <= 5.0:
                    hits.append((dd, f"{rn}{ri}:{an}", o))
    print(f"    {label} ({'/'.join(gatoms)}):")
    if not hits:
        print("      none <= 5.0 A")
        continue
    seen = set()
    for dd, who, o in sorted(hits):
        mark = "  SALT BRIDGE / H-BOND" if dd <= 3.5 else ""
        print(f"      {who:<14} -- {o}  {dd:6.2f} A{mark}")
        seen.add(who.split(":")[0])
