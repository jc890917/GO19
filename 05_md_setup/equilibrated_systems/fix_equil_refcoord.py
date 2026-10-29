#!/usr/bin/env python3
"""Patch 07_equilibrate.sh for the posres + pressure-coupling bug.

Usage: python3 fix_equil_refcoord.py 07_equilibrate.sh [--apply]

Three changes:

1. refcoord-scaling.  With position restraints active, the barostat rescales
   atom coordinates but NOT the posres reference coordinates unless
   refcoord-scaling is set.  The restrained atoms are then displaced relative
   to their references on every box scaling, which spikes the virial, which
   makes the barostat overcorrect -- a divergent oscillation.  Observed in the
   WT NPT: box volume flapping 592 <-> 688 nm^3, pressure +-25000 bar,
   T = 395 K against ref_t = 300, conserved energy running from -8e5 to +2.8e8,
   LINCS failure and SIGSEGV at 376 ps.

   grompp warns about this explicitly:
     "You are using pressure coupling with absolute position restraints,
      this will give artifacts. Use the refcoord_scaling option."

   The line is added in run_stage() rather than in each heredoc, so it applies
   to every pressure-coupled stage, including any added later.  'com' scales
   the reference set's centre of mass with the box without distorting its
   internal geometry.

2. -maxwarn 5 -> 1.  The warning above was suppressed by -maxwarn 5 and a
   six-hour run was lost.  1 permits the known benign net-charge warning
   (-0.002 e from AM1-BCC rounding) and stops the job on anything new.

3. Stage skip.  run_stage returns early if equil/<stage>.gro already exists,
   so the script can be resubmitted to redo a failed NPT without repeating a
   good NVT (~4 h per system).
"""
import sys, re, shutil

path = sys.argv[1]
apply_ = "--apply" in sys.argv
src = open(path).read()
orig = src
changes = []

# ---- 1. stage skip guard -------------------------------------------------
anchor = '    echo ""\n    echo "  --- $STAGE ---"\n'
guard = anchor + '''
    # [6] Skip a stage that already produced its .gro, so a failed later stage
    # can be redone without repeating a good earlier one.
    if [[ -f "equil/${STAGE}.gro" ]]; then
        echo "  = $STAGE already complete (equil/${STAGE}.gro exists) -- skipping"
        return 0
    fi
'''
if 'already complete (equil/${STAGE}.gro' in src:
    changes.append("stage-skip guard: already present")
elif anchor in src:
    src = src.replace(anchor, guard, 1)
    changes.append("stage-skip guard: INSERTED")
else:
    changes.append("stage-skip guard: ANCHOR NOT FOUND -- check manually")

# ---- 2. refcoord-scaling -------------------------------------------------
anchor2 = '''    if [[ -n "$DEFINE" ]]; then
        echo "define = $DEFINE" >> "$MDP"
    fi
'''
add2 = anchor2 + '''
    # [7] Position restraints + pressure coupling REQUIRE refcoord-scaling.
    # Without it the barostat rescales coordinates but not the restraint
    # reference, driving a divergent box oscillation. grompp warns; -maxwarn
    # used to hide it.
    sed -i '/^refcoord[-_]scaling[[:space:]]*=/d' "$MDP"
    if [[ -n "$DEFINE" ]] && grep -qiE '^[[:space:]]*pcoupl[[:space:]]*=[[:space:]]*(berendsen|c-rescale|parrinello-rahman|mttk)' "$MDP"; then
        echo "refcoord-scaling = com" >> "$MDP"
        echo "  [refcoord-scaling = com added -- posres with pressure coupling]"
    fi
'''
if 'refcoord-scaling = com' in src:
    changes.append("refcoord-scaling: already present")
elif anchor2 in src:
    src = src.replace(anchor2, add2, 1)
    changes.append("refcoord-scaling: INSERTED")
else:
    changes.append("refcoord-scaling: ANCHOR NOT FOUND -- check manually")

# ---- 3. maxwarn ----------------------------------------------------------
n = len(re.findall(r'-maxwarn\s+5', src))
if n:
    src = re.sub(r'-maxwarn\s+5', '-maxwarn 1', src)
    changes.append(f"-maxwarn 5 -> 1: {n} occurrence(s)")
elif '-maxwarn 1' in src:
    changes.append("-maxwarn: already 1")
else:
    changes.append("-maxwarn: no '-maxwarn 5' found -- check manually")

print(f"  {path}")
for c in changes:
    print(f"    {c}")

if src == orig:
    print("\n  nothing to change.")
elif apply_:
    shutil.copy(path, path + ".prerefcoord.bak")
    open(path, "w").write(src)
    print(f"\n  APPLIED; original saved as {path}.prerefcoord.bak")
else:
    print("\n  dry run. Re-run with --apply")
