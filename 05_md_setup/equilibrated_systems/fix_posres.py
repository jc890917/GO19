#!/usr/bin/env python3
"""Insert the position-restraint include blocks that posres_setup.py never did.

posres_setup.py's insertion test was unreachable: it required current_mol to be
set on a '[ moleculetype ]' line, but the branch above it cleared current_mol on
exactly those lines. It printed success and changed nothing, in every build.

This walks the moleculetype blocks properly and inserts, at the END of each
target molecule's block:
      #ifdef POSRES
      #include "posre_<name>.itp"
      #endif

Usage:  python3 fix_posres.py gromacs/complex.top [--apply]
Without --apply it only reports what it would do.
"""
import re, sys, shutil, os

TOP = sys.argv[1]
APPLY = "--apply" in sys.argv
TARGETS = {"system1", "FAD", "LIG"}

lines = open(TOP).readlines()
if any('#ifdef POSRES' in l for l in lines):
    print("  already contains '#ifdef POSRES' — nothing to do"); sys.exit(0)

mt  = [n for n, l in enumerate(lines) if re.match(r'^\s*\[\s*moleculetype\s*\]', l)]
sysl= [n for n, l in enumerate(lines) if re.match(r'^\s*\[\s*system\s*\]', l)]
if not mt: sys.exit("  no [ moleculetype ] sections found")
ends = mt[1:] + [sysl[0] if sysl else len(lines)]

blocks = []
for s, e in zip(mt, ends):
    name = None
    for n in range(s + 1, e):
        t = lines[n].strip()
        if not t or t.startswith(';'): continue
        if t.startswith('['): break
        name = t.split()[0]; break
    blocks.append((s, name, e))

print(f"  {TOP}: {len(blocks)} moleculetype blocks")
todo = []
for s, name, e in blocks:
    mark = ""
    if name in TARGETS:
        itp = os.path.join(os.path.dirname(TOP) or ".", f"posre_{name}.itp")
        if os.path.exists(itp):
            todo.append((s, name, e)); mark = "  <- insert posres here"
        else:
            mark = f"  <- TARGET but {itp} missing!"
    print(f"    lines {s+1:>7}-{e:<7} {str(name):<10}{mark}")

if not todo: sys.exit("  nothing to insert")
if not APPLY:
    print(f"\n  would insert {len(todo)} block(s). Re-run with --apply"); sys.exit(0)

shutil.copy(TOP, TOP + ".bak")
for s, name, e in reversed(todo):          # back to front keeps indices valid
    lines[e:e] = ['\n', f'; Position restraints for {name}\n',
                  '#ifdef POSRES\n', f'#include "posre_{name}.itp"\n',
                  '#endif\n', '\n']
open(TOP, 'w').writelines(lines)
print(f"\n  inserted {len(todo)} block(s); original saved as {os.path.basename(TOP)}.bak")
