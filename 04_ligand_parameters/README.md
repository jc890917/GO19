# CP3 — AM1-BCC vs RESP for glyphosate (the cheap comparison)

Goal: find out whether the charge method measurably changes this ligand, before
spending any MD time on the answer. 16 atoms, so the whole thing is minutes of
compute. The output is a table for the response letter either way.

## Run order

```bash
# 1. login node — answers the two open questions. Send me probe_report.txt.
bash 00_probe_env.sh

# 2. put the CPU partition from probe_report.txt into 01_run_gaussian.slurm
#    (replace CHANGE_ME), then:
sbatch 01_run_gaussian.slurm          # array of 3, one per conformer

# 3. when the array finishes — MANDATORY check, do not skip
python3 02_check_geometry.py confA.log confB.log confC.log

# 4. login node — builds both charge sets on the same geometry
bash 03_derive_charges.sh

# 5. the answer
python3 04_compare_charges.py /LUSTRE/.../alumnos/ian/go19/LIG.mol2
```

Step 5's optional argument is the **original mis-protonated** `LIG.mol2`, so the
report also shows how the AM1-BCC-vs-RESP gap compares with the protonation error
we already found. That comparison is the point: it tells us whether this
checkpoint matters at all relative to the fix already in hand.

## Files

| file | what it is |
|---|---|
| `confA.gjf` | bound conformation, taken from the validated seed102 pose 1 |
| `confB.gjf`, `confC.gjf` | two diverse free conformers (1.37 and 0.97 Å heavy-atom RMSD from A) |
| `confA_bound.sdf`, `conf{B,C}_free.sdf` | the same geometries, for inspection |
| `LIG_confA.sdf` | input to the AM1-BCC arm — same geometry, fair comparison |
| `00_probe_env.sh` | licence + partition + AmberTools probe |
| `01_run_gaussian.slurm` | PCM optimisation then vacuum MK ESP, per conformer |
| `02_check_geometry.py` | proton-transfer guard |
| `03_derive_charges.sh` | both arms, plus frcmod diff and ATTN check |
| `04_compare_charges.py` | the comparison table and the decision rule |

All three conformers share one atom ordering, which multi-conformer RESP requires:

```
 1 N+   2 C(reacting)  3 C(carboxyl)  4 O   5 O-
 6 C(to P)   7 P   8 O   9 O-  10 O-
11,12 H-N    13,14 H-C(reacting)   15,16 H-C(to P)
```

## Why the protocol looks like this

Four choices that are easy to get wrong, and each one would have produced a
plausible-looking but wrong charge set:

1. **Optimise under PCM(water), not in vacuum.** Glyphosate at pH 7.4 has an
   NH₂⁺ adjacent to a COO⁻. A vacuum optimisation transfers that proton and
   hands back charges for a neutral molecule — the same failure mode as the
   error that forced this rerun. `02_check_geometry.py` exists to catch it.
2. **Fit the ESP in vacuum, at HF/6-31G\*.** RESP and GAFF were both
   parameterised against vacuum HF/6-31G\* ESP. Fitting in solvent gives
   charges that are arguably better but no longer consistent with the force
   field they are used in.
3. **No diffuse functions.** For a −2 anion, `6-31+G*` or larger lets the extra
   electron delocalise into vacuum and the ESP becomes meaningless. Plain
   6-31G\* is the right basis here precisely because it lacks them.
4. **`maxcyc=0` in the AM1-BCC arm.** Otherwise sqm minimises in vacuum and
   transfers the proton, so the two arms would be compared at different
   geometries — and the comparison would be worthless.

Also checked: that the three resonance-equivalent phosphonate oxygens, and the
two carboxylate oxygens, come out with equal charges. RESP's stage 2 should
enforce this from the GAFF2 atom types; if it does not, the fit needs manual
equivalence groups and the script says so.

## Decision rule, fixed before seeing the result

| max \|Δq\| | action |
|---|---|
| ≤ 0.05 e | keep AM1-BCC; rebut Reviewer 2 with the table |
| 0.05 – 0.10 e | use RESP (it is an hour, and we are rebuilding anyway) |
| > 0.10 e | use RESP; the divergence reaches the phosphonate-anchoring claim |

Stated in advance so the conclusion is not chosen after the numbers are in.

## Not doing

Full RESP on FAD (84 atoms). The flavin half of Reviewer 2's comment is being
answered by argument, not by calculation: FAD never changes redox state in these
simulations, the best-documented flavin model is CHARMM36 and so incompatible
with ff19SB/GAFF2/OPC, no citable AMBER oxidised-FAD set exists, FAD is identical
across all four systems so any error is common-mode, and isoalloxazine geometry
is set by bonded parameters rather than ring charges.
