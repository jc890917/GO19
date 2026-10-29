# 06_md_analysis/scripts

**TO ADD before release** — copy from the Ubuntu analysis PC (`~/GO19_MD/`); see `../../tools/collect_from_analysis_pc.sh`:

- `00_fix_pbc` processing (`-pbc mol -center` → `-pbc nojump` → backbone fit)
- `01_rmsd_backbone.py`, `02_rmsd_fad_lig.py`, `03_rg.py`, `04_contactos_lig_heavy.py`, `05_interacciones_prolif.py`, `06_dssp.py`, `07_pca_v2.py`, `10_hydride_v2.py`, `11_r53_sidechain_rmsd.py`, `run_one_analysis.sh`
- MM-GBSA: `MMPBSA.py` input file(s), `ante-MMPBSA.py` commands, `merge_mmgbsa_0_100.py`
- MDpocket: `quantify_mdpocket.py`, `o2_access.py`

Already here: `nac_sensitivity.py`.
