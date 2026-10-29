import os, sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
RDLogger.DisableLog('rdApp.*')

D = os.path.expanduser("~/mnt/GOx/for frontiers/CP2_redocking")
SMI = "[NH2+](CC(=O)[O-])CP(=O)([O-])[O-]"

# ---- reference: correct protonation, canonical atom order ----------------
ref = Chem.AddHs(Chem.MolFromSmiles(SMI))
Chem.AssignStereochemistry(ref, cleanIt=True, force=True)
order = [(a.GetIdx(), a.GetSymbol(), a.GetFormalCharge()) for a in ref.GetAtoms()]
print("reference formula:", Chem.rdMolDescriptors.CalcMolFormula(ref))
print("net charge:", Chem.GetFormalCharge(ref), " atoms:", ref.GetNumAtoms())

# ---- conformer A : the validated docked pose (heavy atoms fixed) ---------
pose = Chem.MolFromMolFile(os.path.join(D, "seed102", "pose1.sdf"), removeHs=False)
print("docked pose atoms (GNINA strips nonpolar H):", pose.GetNumAtoms())
poseH = Chem.AddHs(pose, addCoords=True)
print("after AddHs:", poseH.GetNumAtoms(), " charge:", Chem.GetFormalCharge(poseH))

# optimise ONLY the added hydrogens, heavy atoms and polar H frozen
mp = AllChem.MMFFGetMoleculeProperties(poseH)
if mp is None:
    print("!! MMFF typing failed on the pose"); sys.exit(1)
ff = AllChem.MMFFGetMoleculeForceField(poseH, mp)
nfix = 0
for a in poseH.GetAtoms():
    if a.GetAtomicNum() > 1:
        ff.AddFixedPoint(a.GetIdx()); nfix += 1
    elif a.GetNeighbors()[0].GetAtomicNum() == 7:      # polar H on N, keep docked
        ff.AddFixedPoint(a.GetIdx()); nfix += 1
ff.Minimize(maxIts=2000)
print(f"conf A: relaxed {poseH.GetNumAtoms()-nfix} added H, {nfix} atoms frozen")

# map the pose onto the reference atom order
match = poseH.GetSubstructMatch(ref, useChirality=False)
if not match:
    match = ref.GetSubstructMatch(poseH)
    print("!! reverse match used:", bool(match))
print("substruct match length:", len(match))

def write_gjf(mol, confId, name, title):
    conf = mol.GetConformer(confId)
    lines = []
    for i, a in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        lines.append(f" {a.GetSymbol():<2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}")
    body = "\n".join(lines)
    txt = f"""%chk={name}.chk
%mem=8GB
%nprocshared=8
#p HF/6-31G* Opt SCF=Tight SCRF=(IEFPCM,Solvent=Water)

{title} - step 1: PCM(water) optimisation (holds the zwitterion)

-2 1
{body}

--Link1--
%chk={name}.chk
%mem=8GB
%nprocshared=8
#p HF/6-31G* SCF=Tight Pop=MK IOp(6/33=2,6/41=10,6/42=17) Geom=AllCheck Guess=Read

"""
    open(f"{name}.gjf", "w").write(txt)
    print("  wrote", f"{name}.gjf")

# reorder conf A into reference order
confA = Chem.RWMol(ref)
confA.RemoveAllConformers()
c = Chem.Conformer(ref.GetNumAtoms())
src = poseH.GetConformer()
for refIdx, poseIdx in enumerate(match):
    c.SetAtomPosition(refIdx, src.GetAtomPosition(poseIdx))
confA.AddConformer(c, assignId=True)
confA = confA.GetMol()
Chem.MolToMolFile(confA, "confA_bound.sdf")
write_gjf(confA, 0, "confA", "glyphosate dianion, BOUND conformation (seed102 pose 1)")

# ---- conformers B, C : diverse free conformers ---------------------------
free = Chem.Mol(ref)
cids = AllChem.EmbedMultipleConfs(free, numConfs=60, randomSeed=20260927,
                                  useSmallRingTorsions=True)
print("embedded", len(cids), "conformers")
# pick the two most different from A and from each other (heavy-atom RMSD)
probe = Chem.Mol(confA)
def rms(m1, c1, m2, c2):
    return rdMolAlign.GetBestRMS(Chem.RemoveHs(Chem.Mol(m1)), Chem.RemoveHs(Chem.Mol(m2)),
                                 prbId=c1, refId=c2)
scored = sorted(((rms(free, cid, probe, 0), cid) for cid in cids), reverse=True)
pickB = scored[0][1]
print(f"conf B: cid {pickB}  RMSD to A = {scored[0][0]:.2f} A")
rest = [(min(rms(free, cid, probe, 0), rms(free, cid, free, pickB)), cid)
        for _, cid in scored[1:]]
pickC = max(rest)[1]
print(f"conf C: cid {pickC}  min RMSD to A/B = {max(rest)[0]:.2f} A")

for nm, cid in (("confB", pickB), ("confC", pickC)):
    m = Chem.Mol(free); m.RemoveAllConformers()
    m.AddConformer(free.GetConformer(cid), assignId=True)
    Chem.MolToMolFile(m, f"{nm}_free.sdf")
    write_gjf(m, 0, nm, f"glyphosate dianion, free conformer {nm[-1]}")

# ---- mol2 of conf A for the AM1-BCC arm (same geometry, fair comparison) --
Chem.MolToMolFile(confA, "LIG_confA.sdf")
print("\nAtom order (shared by all conformers):")
for i, a in enumerate(confA.GetAtoms(), 1):
    print(f"  {i:2d} {a.GetSymbol():<2s} q={a.GetFormalCharge():+d}  "
          f"nbrs={''.join(sorted(n.GetSymbol() for n in a.GetNeighbors()))}")
