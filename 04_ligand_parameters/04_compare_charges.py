#!/usr/bin/env python3
"""
04_compare_charges.py — the deliverable for Reviewer 2, point 3.

Compares, on the SAME corrected geometry:
    AM1-BCC   vs   3-conformer two-stage RESP
and, for scale, against the ORIGINAL (mis-protonated) LIG.mol2 that the
published simulations used.

Usage:  python3 04_compare_charges.py [path/to/original/LIG.mol2]
"""
import sys, os, math

def read_mol2(path):
    atoms = []   # (name, type, charge, element)
    with open(path) as fh:
        sec = None
        for ln in fh:
            if ln.startswith("@<TRIPOS>"):
                sec = ln.strip().split(">")[1]; continue
            if sec != "ATOM":
                continue
            t = ln.split()
            if len(t) < 9:
                continue
            name, atype, q = t[1], t[5], float(t[8])
            el = atype.split(".")[0].upper()
            el = "P" if el.startswith("P") and el != "PD" else el[0]
            atoms.append((name, atype, q, el))
    return atoms

def fmt(v, w=8, p=4):
    return f"{v:{w}.{p}f}"

BCC, RESP = "LIG_bcc.mol2", "LIG_resp.mol2"
for f in (BCC, RESP):
    if not os.path.exists(f):
        sys.exit(f"missing {f} — run 03_derive_charges.sh first")

b, r = read_mol2(BCC), read_mol2(RESP)
if len(b) != len(r):
    sys.exit(f"atom count mismatch: {len(b)} vs {len(r)}")

# Canonical atom order, set in gen_conformers.py and shared by every file:
#  1 N+   2 C(react)  3 C(coo)  4 O  5 O-  6 C(P)  7 P  8 O  9 O-  10 O-
# 11,12 H-N   13,14 H-C(react)   15,16 H-C(P)
GROUPS = {
    "phosphonate  PO3": [6, 7, 8, 9],
    "  P only":         [6],
    "  the 3 O":        [7, 8, 9],
    "carboxylate COO":  [2, 3, 4],
    "  C only":         [2],
    "  the 2 O":        [3, 4],
    "ammonium NH2+":    [0, 10, 11],
    "CH2 (reacting)":   [1, 12, 13],
    "CH2 (to P)":       [5, 14, 15],
}

print("=" * 78)
print("PER-ATOM CHARGES — AM1-BCC vs RESP, identical geometry (bound conformer)")
print("=" * 78)
print(f"{'#':>3} {'name':<6}{'type':<7}{'AM1-BCC':>10}{'RESP':>10}{'RESP-BCC':>11}")
dmax, dsum2 = 0.0, 0.0
for i, (ab, ar) in enumerate(zip(b, r), 1):
    d = ar[2] - ab[2]
    dmax = max(dmax, abs(d)); dsum2 += d * d
    flag = "  <<" if abs(d) >= 0.10 else ""
    print(f"{i:>3} {ab[0]:<6}{ab[1]:<7}{fmt(ab[2],10)}{fmt(ar[2],10)}{fmt(d,11)}{flag}")
print(f"{'':>3} {'TOTAL':<13}{fmt(sum(a[2] for a in b),10)}{fmt(sum(a[2] for a in r),10)}")
print(f"\n  max |Δq| = {dmax:.4f} e      RMS Δq = {math.sqrt(dsum2/len(b)):.4f} e")

print()
print("=" * 78)
print("GROUP CHARGES  (what actually drives a salt bridge)")
print("=" * 78)
print(f"{'group':<20}{'AM1-BCC':>10}{'RESP':>10}{'Δ':>10}")
for g, idx in GROUPS.items():
    qb = sum(b[i][2] for i in idx); qr = sum(r[i][2] for i in idx)
    print(f"{g:<20}{fmt(qb,10)}{fmt(qr,10)}{fmt(qr-qb,10)}")

# equivalence check: chemically identical oxygens must carry identical charges
print()
print("EQUIVALENCE CHECK (resonance-equivalent atoms must match)")
for label, idx in (("phosphonate O (3)", [7, 8, 9]), ("carboxylate O (2)", [3, 4])):
    for nm, s in (("AM1-BCC", b), ("RESP", r)):
        qs = [s[i][2] for i in idx]
        spread = max(qs) - min(qs)
        verdict = "ok" if spread < 0.01 else ("MINOR" if spread < 0.05 else "NOT EQUIVALENT")
        print(f"  {label:<20}{nm:<9} spread = {spread:.4f} e   [{verdict}]")

# scale of the protonation fix, for context
if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
    old = read_mol2(sys.argv[1])
    print()
    print("=" * 78)
    print(f"FOR SCALE — the original mis-protonated set ({os.path.basename(sys.argv[1])})")
    print("=" * 78)
    print(f"  {len(old)} atoms vs {len(b)} — different structure, so only group sums compare.")
    # group by element+neighbourhood is unreliable across files; sum by element
    def by_el(s):
        out = {}
        for nm, ty, q, el in s:
            out[el] = out.get(el, 0.0) + q
        return out
    eo, eb, er = by_el(old), by_el(b), by_el(r)
    print(f"\n  {'element':<10}{'ORIGINAL':>12}{'AM1-BCC':>12}{'RESP':>12}")
    for el in sorted(set(eo) | set(eb)):
        print(f"  {el:<10}{fmt(eo.get(el,0),12)}{fmt(eb.get(el,0),12)}{fmt(er.get(el,0),12)}")
    print(f"  {'TOTAL':<10}{fmt(sum(eo.values()),12)}{fmt(sum(eb.values()),12)}{fmt(sum(er.values()),12)}")
    print("\n  The ORIGINAL column is AM1-BCC computed on NEUTRAL glyphosate with the")
    print("  sum forced to -2 by '-nc -2'. Compare its O total against either corrected")
    print("  column: that gap is the error the rerun removes, and it is the number to")
    print("  weigh the AM1-BCC-vs-RESP gap against.")

# what the difference is worth, in energy
print()
print("=" * 78)
print("WHAT THE DIFFERENCE IS WORTH")
print("=" * 78)
dO = max(abs(r[i][2] - b[i][2]) for i in (7, 8, 9))
q_arg_H = 0.35          # typical ff19SB charge on an Arg guanidinium H
rsep = 2.9              # a tight salt-bridge H...O distance, A
for eps in (1.0, 4.0, 20.0):
    dE = 332.06 * dO * q_arg_H / (eps * rsep)
    print(f"  largest phosphonate-O Δq = {dO:.4f} e  ->  ΔE ≈ {dE:5.2f} kcal/mol "
          f"per H-bond at r={rsep} Å, ε={eps:g}")
print("  ε=1 is the unscreened upper bound; a solvent-exposed salt bridge in GB")
print("  behaves closer to the ε=4-20 rows. Compare against the replica-to-replica")
print("  spread of the published MM-GBSA, ±6.35 kcal/mol.")

print()
print("=" * 78)
print("DECISION RULE (set in advance, so the answer is not chosen after the fact)")
print("=" * 78)
if dmax <= 0.05:
    print("  max |Δq| ≤ 0.05 e  ->  AM1-BCC is ADEQUATE for this ligand.")
    print("  Keep AM1-BCC, and answer Reviewer 2 with this table. The concern is")
    print("  real in general but does not apply measurably here.")
elif dmax <= 0.10:
    print("  0.05 < max |Δq| ≤ 0.10 e  ->  BORDERLINE.")
    print("  RESP is nearly free for a 16-atom ligand and the systems are being")
    print("  rebuilt anyway, so use RESP and report both. No reason to defend the")
    print("  weaker input when the stronger one costs an hour.")
else:
    print("  max |Δq| > 0.10 e  ->  USE RESP.")
    print("  The divergence is large enough to change salt-bridge strengths, which")
    print("  is exactly where the paper's phosphonate-anchoring claim lives.")
print()
print("  Either way: this table goes in the response letter. It converts an")
print("  assertion into a measurement, which is the only kind of rebuttal that")
print("  closes a methods objection.")
