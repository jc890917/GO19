# GO19 glyphosate oxidoreductase — structural and evolutionary analysis

Data, input files and scripts supporting:

> Angel-Lerma LE, Contreras-Montoya E, González-Escobedo R, Moncayo-Luján R, Barraza-Guerrero SI and Carrillo-Campos J (2026). **Structural and evolutionary basis of glyphosate specificity in a naturally occurring glyphosate oxidoreductase.** *Frontiers in Chemical Biology* 5:1996804. doi: 10.3389/fchbi.2026.1996804

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23754903.svg)](https://doi.org/10.5281/zenodo.23754903)

Archived on Zenodo — cite all versions: https://doi.org/10.5281/zenodo.23754903 (v1.0: https://doi.org/10.5281/zenodo.23754904)
Molecular dynamics trajectories (Zenodo, separate record): DOI to be added

GO19 sequences: GenBank PQ065997 (nucleotide) and XER92051 (protein).

## Contents

| Folder | What it contains |
|---|---|
| `01_sequences/` | GO19 protein sequence (XER92051) |
| `02_structure_prediction/` | AlphaFold 3 models of GO19 + FAD and of the reference glycine oxidases / glyphosate oxidoreductases (model 0–4 for GO19, model 0 for references, with confidence summaries); TM-score/RMSD comparison with the 1RYI crystal structure |
| `03_docking/` | GNINA receptor and ligand, docking scripts, poses and logs for the 10 independent seeds, and the pose-reproducibility analysis |
| `04_ligand_parameters/` | Glyphosate (net charge −2) conformers, AM1-BCC and RESP charge derivation and comparison |
| `05_md_setup/` | Wild-type system build (tleap inputs, GAFF2 parameters for glyphosate and FAD, ff19SB/OPC topology, AMBER→GROMACS conversion), mutant build scripts (R53A, R272L, R53A/R272L), equilibrated systems with all `.mdp` files and topologies for the four systems, and the production driver (`production/gpu_run.sh`, which writes the 100-ns production `.mdp`) |
| `06_md_analysis/` | Per-replicate analysis outputs (RMSD, Rg, DSSP, PCA, hydride-transfer geometry, ProLIF fingerprints, contacts), MM-GBSA per-replicate and per-residue results (0–100 ns), NAC sensitivity, MDpocket summary, final tables (Table 1A/1B and supplementary tables), representative structures, and the figure notebook |
| `07_phylogeny_conservation/` | Sequence dataset, alignment, IQ-TREE 3 tree and active-site conservation tables |

## Software

AlphaFold 3 Server · GNINA (version in `03_docking/gnina_version.txt`) · AmberTools 19 (tleap, antechamber, pdb4amber) · AmberTools 26 (MMPBSA.py, ante-MMPBSA.py) · ParmEd · GROMACS 2025.4 (GPU) · MDAnalysis 2.10.0 · ProLIF 2.2.3 · RDKit · fpocket 4.0 (mdpocket) · BLASTp · MAFFT · IQ-TREE 3 · TM-align · Python 3 (NumPy, pandas, matplotlib).

## MD systems

Four systems (wild type, R53A, R272L, R53A/R272L) × three independent 100-ns replicates. Folder names: `go19_v2` (wild type), `go19_v2_R53A`, `go19_v2_R272L`, `go19_v2_R53A_R272L`; replicate folders `prod_rep1–3`. All MD averages are computed over the full 0–100 ns production runs.

The trajectories are not stored in this repository because of their size. Solvent-stripped trajectories (protein + FAD + glyphosate, every 100 ps) are available in the Zenodo trajectory record above.

## Licences

- Scripts: MIT licence (`LICENSE`).
- Data files produced by the authors: CC BY 4.0 (`LICENSE-DATA.md`).
- **AlphaFold 3 outputs and files derived from them** (folder `02_structure_prediction/`, and the receptor and system coordinates built from the GO19 model) are subject to the [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms) (non-commercial use only); see the `terms_of_use.md` files.
- `02_structure_prediction/structural_comparison/1RYI.pdb` is from the Protein Data Bank (Settembre et al., 2003).

## Contact

Javier Carrillo-Campos — Facultad de Zootecnia y Ecología, Universidad Autónoma de Chihuahua, Mexico — jccarrillo@uach.mx
