#!/usr/bin/env python3
"""Verify the PCM optimisation did NOT transfer the ammonium proton.
Run after the Gaussian jobs:  python3 02_check_geometry.py confA.log confB.log confC.log
A zwitterion that collapsed to the neutral form invalidates the whole charge set —
this is the same class of error that forced the rerun, so it is checked explicitly.
"""
import sys, re, math

def last_geometry(path):
    """final 'Standard orientation' / 'Input orientation' block"""
    txt = open(path, errors="replace").read()
    blocks = re.findall(
        r"(?:Standard|Input) orientation:.*?-{20,}\n(.*?)\n\s*-{20,}",
        txt, re.S)
    if not blocks:
        return None
    rows = []
    for ln in blocks[-1].splitlines():
        t = ln.split()
        if len(t) >= 6 and t[0].isdigit():
            rows.append((int(t[1]), float(t[-3]), float(t[-2]), float(t[-1])))
    return rows

Z = {1: "H", 6: "C", 7: "N", 8: "O", 15: "P"}
def d(a, b):
    return math.sqrt(sum((a[i+1]-b[i+1])**2 for i in range(3)))

ok_all = True
for path in sys.argv[1:]:
    g = last_geometry(path)
    print(f"\n=== {path} ===")
    if not g:
        print("  !! no geometry block found — job probably failed"); ok_all = False; continue
    N = [a for a in g if a[0] == 7]
    O = [a for a in g if a[0] == 8]
    H = [a for a in g if a[0] == 1]
    if len(N) != 1 or len(O) != 5 or len(H) != 6:
        print(f"  !! unexpected composition N={len(N)} O={len(O)} H={len(H)}"); ok_all = False
    nH = [h for h in H if d(h, N[0]) < 1.25]
    oH = [(h, min(O, key=lambda o: d(h, o))) for h in H]
    oH = [(h, o) for h, o in oH if d(h, o) < 1.25]
    print(f"  H on N : {len(nH)}  (need 2)   distances: "
          + ", ".join(f"{d(h,N[0]):.3f}" for h in nH))
    print(f"  H on O : {len(oH)}  (need 0)   "
          + (", ".join(f"{d(h,o):.3f}" for h, o in oH) if oH else "-"))
    if len(nH) == 2 and not oH:
        print("  >>> zwitterion INTACT — geometry usable")
    else:
        print("  >>> PROTON TRANSFER DETECTED — DO NOT USE THIS CONFORMER")
        print("      (re-run with the N-H distances frozen, or drop this conformer)")
        ok_all = False
print("\nALL CONFORMERS OK" if ok_all else "\n*** AT LEAST ONE CONFORMER IS BAD ***")
sys.exit(0 if ok_all else 1)
