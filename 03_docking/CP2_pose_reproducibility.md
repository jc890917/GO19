# CP2 — Redocking reproducibility, corrected ligand (CLOSED)

**Receptor** `go19_receptor_fixed.pdb` (spurious FAD N1–H removed)
**Ligand** `glifosato_m2.sdf` — net −2 at pH 7.4 (NH₂⁺, COO⁻, PO₃²⁻), built from SMILES
**Box** 20³ Å centred exactly on FAD N5 (−2.857, 1.645, 4.710)
**GNINA** 1.3.2, exhaustiveness 32, 20 modes, `--cnn_scoring rescore`, `--flexres A:257`
**Seeds** 42, 43, 101, 102, 202, 203, 303, 304, 404, 405 (10 independent runs, 200 poses)
**Reacting carbon** = glycine-moiety CH₂, the hydride donor to N5

---

## Verdict

Your claim holds, and in a stronger form than stated: in **6 of 10 seeds the catalytically
competent pose is rank 1**, and the 4 seeds that miss it are identifiable from the CNN score
alone, without looking at the geometry.

### Rank-1 pose of every seed

| seed | CNN | affinity | C–N5 | C–C4X | H···N5 | C–H···N5 | P–N5 | RMSD vs seed102 | Y257 OH shift | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| 101 | 0.9255 | −4.36 | 3.73 | 3.76 | 3.29 | 105° | 6.54 | 2.41 | 0.12 | **catalytic** |
| 303 | 0.9079 | −4.50 | 3.49 | 3.71 | 2.56 | 144° | 7.06 | 1.60 | 0.13 | **catalytic** |
| 102 | 0.8674 | −4.40 | 3.48 | 3.67 | 2.60 | 138° | 7.15 | 0.00 | 0.05 | **catalytic** |
| 404 | 0.8427 | −4.46 | 3.64 | 3.91 | 3.06 | 114° | 6.57 | 1.56 | 0.14 | **catalytic** |
| 202 | 0.7943 | −4.41 | 3.41 | 3.69 | 2.52 | 138° | 6.39 | 1.86 | 0.09 | **catalytic** |
| 405 | 0.7519 | −4.74 | 3.51 | 3.70 | 2.60 | 140° | 6.15 | 2.26 | 0.04 | **catalytic** |
| 304 | 0.5691 | −4.53 | 9.06 | 8.11 | 8.94 | 93° | 8.80 | 6.17 | 0.03 | decoy |
| 42  | 0.5516 | −4.52 | 9.03 | 8.08 | 8.87 | 95° | 8.68 | 6.23 | 0.01 | decoy |
| 203 | 0.5462 | −5.05 | 8.90 | 7.96 | 8.64 | 100° | 9.08 | 6.59 | 0.89 | decoy |
| 43  | 0.5438 | −5.10 | 8.92 | 7.98 | 8.72 | 97° | 8.95 | 6.62 | 0.02 | decoy |

H···N5 is from tetrahedrally rebuilt hydrogens (GNINA strips non-polar H from its output).
No pose in any run has a heavy-atom clash below 2.71 Å.

### The CNN score is a clean discriminator

```
CNN ≥ 0.75  →  6/10 seeds  →  ALL catalytic   (C–N5 3.41–3.73 Å)
CNN ≤ 0.57  →  4/10 seeds  →  NONE catalytic  (C–N5 8.90–9.06 Å)
gap: 0.7519 vs 0.5691 — bimodal, no overlap
```

Note the inversion: the four decoy runs have the *better* Vina affinities (−4.52 to −5.10)
and the *worse* CNN scores. The CNN rescoring is what recognises the catalytic pose; the
classical score does not. That is a usable, honest selection rule.

### Consensus ensemble (the 6 catalytic rank-1 poses)

| quantity | value |
|---|---|
| C(react)···N5 | **3.54 ± 0.12 Å** (3.41–3.73) |
| C(react)···C4X | 3.74 ± 0.09 Å |
| H···N5 | 2.77 ± 0.32 Å |
| C–H···N5 angle | 130 ± 16° |
| P···N5 | 6.65 ± 0.39 Å |
| pairwise heavy-atom RMSD | 1.97 ± 0.37 Å (1.37–2.64) |
| Tyr257 OH displacement | **0.09 ± 0.04 Å** |

Original MD starting pose: C–N5 3.55 Å. The consensus reproduces it to within 0.01 Å of the
mean — with the corrected −2 ligand and the corrected receptor. **The pose the MD was built
on was right; only the ligand protonation was wrong.** That is the best possible outcome for
the rerun: the binding mode does not have to be re-argued, only re-parameterised.

### Invariant carboxylate anchoring

In all six catalytic poses the carboxylate is held by the same four partners, tightly clustered:

| partner | range across the 6 seeds |
|---|---|
| Arg314 NE | 2.79 – 3.02 Å |
| Arg314 NH2 | 2.84 – 3.01 Å |
| Arg272 NH2 | 2.88 – 3.08 Å |
| Tyr257 OH | 3.03 – 3.08 Å |
| Pro341 O (backbone) | 3.24 – 3.43 Å |

Plus the ammonium to Pro341 O at 2.80–3.21 Å in five of six.

### Tyr257 is not a gate

OH moves 0.04–0.14 Å in every catalytic pose — a rounding error. Y257 is oblique to the
flavin and acts as the fourth carboxylate anchor, not as an obstruction. **No relaxation MD
is needed**, and rigid docking would very likely give the same answer. Worth one control run
(drop `--flexres` entirely) purely to simplify the Methods.

### The phosphonate is not anchored at t = 0

| seed | phosphonate partner ≤3.5 Å | nearest Arg53 N |
|---|---|---|
| 101 | Arg272 NH2 2.89 | NH1 5.75 |
| 303 | Arg272 NH2 3.00 | NH1 5.03 |
| 405 | FAD N3 3.14, Pro341 O 3.20 | NH1 7.41 |
| 102 | none | NH1 5.23 |
| 202 | none | NH1 6.16 |
| 404 | none | NH1 5.85 |

Arg53 is **5.0–7.4 Å away in all six**. Arg272 bridges to the phosphonate in only 2 of 6, and
in those it is the *same* NH2 already bonded to the carboxylate — a shared, not a dedicated,
contact.

The manuscript's R53 + R272 phosphonate-anchoring claim has no support in the docked
complex. It either forms during MD (a specific, testable prediction for the rerun: track
Arg53 NH1/NH2 to phosphonate O over 3 × 100 ns) or it must be removed. Reviewer 1's
point 2 is about exactly this inconsistency, so it needs a definite answer, not a hedge.

### The decoy site

All four failing seeds land in the same place: Ser47 / Pro48 / Asn51 / Thr345, two of them
extending to Thr50 / Gln91 / Ser251. This is a reproducible surface anion site ~9 Å from the
flavin. It can be reported as a probable peripheral/entry site rather than discarded.

---

## Superseded

The earlier three-residue runs (`--flexres A:314,A:272,A:257`, now overwritten; analysis kept
at `analysis/pose_analysis_flexres3_SUPERSEDED.txt`) gave **0 of 100 poses** with C–N5 < 4 Å,
best 4.18 Å at CNN rank 10. The catalytic orientation was still *sampled* there (3 of 5 seeds,
ranks 6–20) but buried under decoys. So the difference between the protocols is ranking, not
sampling — which is the honest way to justify `A:257` in the Methods: fewer flexible degrees
of freedom, better discrimination, verified over 10 seeds.

---

## What this settles for the reviewers

- **Reviewer 1 point 3 (docking validation):** answerable with data now — 10 independent
  seeds, 6/10 top-ranked convergence to a 1.97 Å-RMSD consensus, and a score threshold that
  separates success from failure. Report all ten, including the four failures. Still owed:
  the 1RYI glycolate redocking control.
- **Reviewer 1 point 1 / Reviewer 2 point 6 (hydride geometry):** C–N5 3.54 ± 0.12 Å at t = 0
  relaxing to 4.47 ± 0.49 Å over 3 × 100 ns is a near-attack conformation, not a pre-reactive
  one. The Abstract has to say so.
- **Reviewer 1 point 2 (residue assignment):** the docked complex says carboxylate →
  R272/R314/Y257/P341, phosphonate → unanchored. Any phosphonate assignment must come from
  the MD and be reported with occupancies.

---

## Next: CP3 — charges

The protonation is fixed and the pose is validated. The remaining input-quality objection is
Reviewer 2 point 3 (AM1-BCC on phosphonate and isoalloxazine). Gaussian 09 is available on
Yoltla, so RESP is feasible for both glyphosate and FAD. That is the next checkpoint, and it
should be done before any system is rebuilt — otherwise the 12 × 100 ns get run twice.

**Files** — `CP2_redocking/analysis/`: `analyze_poses.py`, `consensus.py`, `thresholds.py`,
`pose_analysis_flexres257.txt`, `consensus_flexres257.txt`, `thresholds_flexres257.txt`,
`pose_analysis_flexres3_SUPERSEDED.txt`. All pure Python, no dependencies.
