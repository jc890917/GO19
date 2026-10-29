"""R2 comment 3: TM-score / RMSD of GO19 AF3 model vs 1RYI, 1NG4 and the other GO AF3 models."""
import re, gemmi, json, glob, os, itertools, numpy as np, csv
from tmtools import tm_align
B = os.path.expanduser('~/mnt/GOx')
def ca(path, chain=None):
    st = gemmi.read_structure(path); st.setup_entities(); m = st[0]
    ch = m[chain] if chain else m[0]
    xyz, seq, b = [], [], []
    for r in ch:
        if r.het_flag == 'H' or 'CA' not in [a.name for a in r]: continue
        a = r['CA'][0]; xyz.append(a.pos.tolist()); b.append(a.b_iso)
        seq.append(gemmi.find_tabulated_residue(r.name).one_letter_code.upper() if gemmi.find_tabulated_residue(r.name) else 'X')
    return np.array(xyz), ''.join(seq), np.array(b)
S = {}
for i in range(5): S[f'GO19_m{i}'] = ca(f'{B}/ThiO_XER92051_FAD/fold_thio_xer92051_fad_model_{i}.cif')
for d,n in [('go_bcereus','GO_Bcereus'),('go_blicheniformis','GO_Blicheniformis'),('gow_crhyseo','GOW_Chryseobacterium'),('gox_caballeronia','GOX_Caballeronia'),('gox_ochro','GOX_Ochrobactrum')]:
    S[n] = ca(glob.glob(f'{B}/gos/{d}/*_model_0.cif')[0])
for c in 'ABCD': S[f'1RYI_{c}'] = ca(f'{B}/gos/1RYI.pdb', c)
S['1NG4_A'] = ca(f'{B}/1NG4.pdb', 'A')
def cmp(a, b):
    (x1,s1,_),(x2,s2,_) = S[a], S[b]
    r = tm_align(x1, x2, s1, s2)
    al = [(p,q) for p,q in zip(r.seqxA, r.seqyA) if p!='-' and q!='-']
    ident = sum(p==q for p,q in al)
    return dict(A=a, B=b, lenA=len(s1), lenB=len(s2), TM_normA=round(r.tm_norm_chain1,3), TM_normB=round(r.tm_norm_chain2,3),
                RMSD=round(r.rmsd,2), aligned=len(al), seqid_aligned=round(100*ident/len(al),1))
rows = []
for ref in ['1RYI_A','1RYI_B','1RYI_C','1RYI_D','1NG4_A','GO_Bcereus','GO_Blicheniformis','GOW_Chryseobacterium','GOX_Caballeronia','GOX_Ochrobactrum']:
    rows.append(cmp('GO19_m0', ref))
for i in range(1,5): rows.append(cmp('GO19_m0', f'GO19_m{i}'))           # AF3 sample consistency
for a,b in itertools.combinations('ABCD',2): rows.append(cmp(f'1RYI_{a}', f'1RYI_{b}'))  # crystal baseline
for n in ['GO_Bcereus','GO_Blicheniformis','GOW_Chryseobacterium','GOX_Caballeronia','GOX_Ochrobactrum','1NG4_A']: rows.append(cmp(n,'1RYI_A'))
with open('tm_rmsd_table.csv','w',newline='') as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
# confidence of GO19 models
print('\nGO19 AF3 confidences')
for i in range(5):
    sc = json.load(open(f'{B}/ThiO_XER92051_FAD/fold_thio_xer92051_fad_summary_confidences_{i}.json'))
    b = S[f'GO19_m{i}'][2]
    print(i, 'pTM', sc['ptm'], 'ipTM', sc['iptm'], 'rank', sc['ranking_score'], 'clash', sc['has_clash'], 'disord', sc['fraction_disordered'],
          'pLDDT mean %.1f' % b.mean(), '>90 %.1f%%' % (100*(b>90).mean()), '<70 %.1f%%' % (100*(b<70).mean()),
          'low-res idx', [k+1 for k in np.where(b<70)[0]][:40])

# --- which AF3 sample was used for MD/docking (go19/protein_raw.pdb)?
raw = ca(f'{B}/go19/protein_raw.pdb')
for i in range(5):
    d = np.sqrt(((raw[0]-S[f'GO19_m{i}'][0])**2).sum(1)).max()
    print(f'protein_raw.pdb vs model_{i}: max CA deviation (no fit) {d:.3f} A')
# --- FAD placement vs 1RYI after Calpha superposition (TM-align rotation)
def fad(path, chain=None):
    st = gemmi.read_structure(path); m = st[0]; out = {}
    for ch in m:
        for r in ch:
            if r.name == 'FAD' and (chain is None or ch.name == chain) and not out:
                out = {a.name: np.array(a.pos.tolist()) for a in r if a.element.name != 'H'}
    return out
f19 = fad(f'{B}/ThiO_XER92051_FAD/fold_thio_xer92051_fad_model_0.cif')
r = tm_align(S['GO19_m0'][0], S['1RYI_A'][0], S['GO19_m0'][1], S['1RYI_A'][1])
for c in 'ABCD':
    fr = fad(f'{B}/gos/1RYI.pdb', c)
    rc = tm_align(S['GO19_m0'][0], S[f'1RYI_{c}'][0], S['GO19_m0'][1], S[f'1RYI_{c}'][1])
    common = sorted(set(f19) & set(fr))
    X = np.array([f19[k] for k in common]) @ rc.u.T + rc.t
    Y = np.array([fr[k] for k in common])
    dd = np.sqrt(((X-Y)**2).sum(1))
    iso = [k for k in common if k in {'N1','C2','O2','N3','C4','O4','C4X','N5','C5X','C6','C7','C7M','C8','C8M','C9','C9A','N10','C10'}]
    di = np.sqrt(((np.array([f19[k] for k in iso]) @ rc.u.T + rc.t - np.array([fr[k] for k in iso]))**2).sum(1))
    print(f'FAD vs 1RYI_{c}: {len(common)} matched heavy atoms, RMSD {np.sqrt((dd**2).mean()):.2f} A (no refit); isoalloxazine {len(iso)} atoms RMSD {np.sqrt((di**2).mean()):.2f} A; N5 shift {np.linalg.norm((f19["N5"]@rc.u.T+rc.t)-fr["N5"]):.2f} A')
# --- GO19 residues unaligned to 1RYI
i=j=0; un=[]
for p,q in zip(r.seqxA, r.seqyA):
    if p!='-': i+=1
    if p!='-' and q=='-': un.append(i)
print('GO19 residues with no 1RYI counterpart:', un)
