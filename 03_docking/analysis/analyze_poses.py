#!/usr/bin/env python3
import sys, math, glob, os
D = os.path.expanduser("~/mnt/GOx/for frontiers/CP2_redocking")
REC = os.path.join(D, "go19_receptor_fixed.pdb")
FLEXRES = set()  # filled per-seed from the flex PDB
def dist(a,b): return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))
def angle(a,b,c):
    v1=[a[i]-b[i] for i in range(3)]; v2=[c[i]-b[i] for i in range(3)]
    n1=math.sqrt(sum(x*x for x in v1)); n2=math.sqrt(sum(x*x for x in v2))
    cs=sum(v1[i]*v2[i] for i in range(3))/(n1*n2)
    return math.degrees(math.acos(max(-1.0,min(1.0,cs))))
rigid=[]; fad={}
for ln in open(REC):
    if not ln.startswith(("ATOM","HETATM")): continue
    name=ln[12:16].strip(); resn=ln[17:20].strip()
    try: resi=int(ln[22:26])
    except ValueError: continue
    xyz=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
    el=ln[76:78].strip() or name[0]
    if resn=="FAD": fad[name]=xyz
    if el=="H": continue
    rigid.append((resn,resi,name,el,xyz))
N5,C4X=fad["N5"],fad["C4X"]
y257_rigid_OH=next(a[4] for a in rigid if a[1]==257 and a[2]=="OH")
def read_sdf(path):
    poses=[]; lines=open(path).read().splitlines(); i=0
    while i < len(lines):
        if i+3 < len(lines) and "V2000" in lines[i+3]:
            natoms=int(lines[i+3][0:3]); atoms=[]
            for j in range(i+4, i+4+natoms):
                L=lines[j]
                atoms.append((L[31:34].strip(),(float(L[0:10]),float(L[10:20]),float(L[20:30]))))
            poses.append(atoms)
            while i < len(lines) and lines[i].strip()!="$$$$": i+=1
        i+=1
    return poses
def read_flex(path):
    models=[]; cur=None
    for ln in open(path):
        if ln.startswith("MODEL"): cur=[]
        elif ln.startswith("ENDMDL"):
            models.append(cur); cur=None
        elif cur is not None and ln.startswith(("ATOM","HETATM")):
            name=ln[12:16].strip(); resn=ln[17:20].strip(); resi=int(ln[22:26])
            xyz=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
            el=ln[76:78].strip() or name[0]
            cur.append((resn,resi,name,el,xyz))
    return models
ref_path=os.path.join(D,"seed102","pose1.sdf")
ref=read_sdf(ref_path)[0] if os.path.exists(ref_path) else None
# GNINA strips non-polar H: 12 atoms
# 0 N+  1 HN  2 HN  3 C_react  4 C_coo  5 O  6 O-  7 C_P  8 P  9 O 10 O- 11 O-
IN,IC,ICOO,ICP,IP=0,3,4,7,8
HEAVY=[0,3,4,5,6,7,8,9,10,11]
CARBOX_O=(5,6); PHOS_O=(9,10,11)
def build_H(c,n1,n2):
    """two tetrahedral H positions on carbon c bonded to n1,n2"""
    import math
    def sub(a,b): return [a[i]-b[i] for i in range(3)]
    def add(a,b): return [a[i]+b[i] for i in range(3)]
    def sc(a,k): return [x*k for x in a]
    def nrm(a):
        m=math.sqrt(sum(x*x for x in a)); return [x/m for x in a]
    def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    u1=nrm(sub(n1,c)); u2=nrm(sub(n2,c))
    bis=nrm(add(u1,u2))          # bisector of the two existing bonds
    perp=nrm(cross(u1,u2))       # normal to the C(n1)(n2) plane
    ang=math.radians(54.75)      # half of 109.5
    h1=[c[i]-bis[i]*1.09*math.cos(ang)+perp[i]*1.09*math.sin(ang) for i in range(3)]
    h2=[c[i]-bis[i]*1.09*math.cos(ang)-perp[i]*1.09*math.sin(ang) for i in range(3)]
    return h1,h2
def contacts(xyz,recatoms,cutoff=3.5,polar_only=True):
    out=[]
    for resn,resi,name,el,c in recatoms:
        if polar_only and el not in ("N","O"): continue
        d=dist(xyz,c)
        if d<=cutoff: out.append((d,f"{resn}{resi}:{name}"))
    return sorted(out)
def rmsd_heavy(p,q): return math.sqrt(sum(dist(p[i][1],q[i][1])**2 for i in HEAVY)/len(HEAVY))
seeds=sorted(int(os.path.basename(f).replace("poses_seed","").replace(".sdf","")) for f in glob.glob(os.path.join(D,"poses_seed*.sdf")))
print("Receptor : go19_receptor_fixed.pdb   FAD N5 = %.3f %.3f %.3f"%N5)
print("Ligand   : glifosato_m2.sdf (net -2)"); print()
allrows={}
for s in seeds:
    poses=read_sdf(os.path.join(D,f"poses_seed{s}.sdf")); flex=read_flex(os.path.join(D,f"flex_seed{s}.pdb"))
    rows=[]
    for k,p in enumerate(poses,1):
        fl=flex[k-1] if k-1<len(flex) else []
        fres={a[1] for a in fl}
        fnames={(a[1],a[2]) for a in fl}
        # keep every rigid atom except those the flex model replaces
        rec=[a for a in rigid if (a[1],a[2]) not in fnames]+[a for a in fl if a[3]!="H"]
        c=p[IC][1]; dC=dist(c,N5); dC4=dist(c,C4X); dP=dist(p[IP][1],N5)
        h1,h2=build_H(c,p[IN][1],p[ICOO][1])
        dH,hbest=min((dist(h1,N5),h1),(dist(h2,N5),h2))
        aCHN=angle(c,hbest,N5)
        minheavy=min(dist(p[i][1],N5) for i in HEAVY)
        clash=min((dist(p[i][1],a[4]),f"{a[0]}{a[1]}:{a[2]}") for i in HEAVY for a in rec)
        oh=next((a[4] for a in fl if a[1]==257 and a[2]=="OH"),None)
        dOH=dist(oh,y257_rigid_OH) if oh else float("nan")
        r=rmsd_heavy(p,ref) if ref else float("nan")
        rows.append(dict(pose=k,dC=dC,dC4=dC4,dP=dP,dH=dH,aCHN=aCHN,minheavy=minheavy,clash=clash,dOH=dOH,rmsd=r,p=p,rec=rec))
    allrows[s]=rows
for s in seeds:
    print("="*100); print(f"SEED {s}"); print("="*100)
    print(" pose  C-N5   C-C4X   H-N5  C-H..N5   P-N5  minHvy  RMSDref  minClash              Y257OH")
    for r in allrows[s]:
        flag="  <<<" if r["dC"]<4.0 else ""
        print(f"  {r['pose']:2d}  {r['dC']:5.2f}  {r['dC4']:5.2f}  {r['dH']:5.2f}  {r['aCHN']:6.1f}  {r['dP']:5.2f}  {r['minheavy']:5.2f}   {r['rmsd']:5.2f}   {r['clash'][0]:4.2f} {r['clash'][1]:<16s} {r['dOH']:5.3f}{flag}")
    b=min(allrows[s],key=lambda x:x["dC"])
    print(f"  -> closest reacting C: pose {b['pose']}  {b['dC']:.2f} A"); print()
print("="*100); print("CONSENSUS - poses with reacting C within 4.0 A of FAD N5"); print("="*100)
print("seed   n<4A   bestpose   bestC-N5   P-N5@best   CNNrank_of_best")
for s in seeds:
    good=[r for r in allrows[s] if r["dC"]<4.0]; b=min(allrows[s],key=lambda x:x["dC"])
    print(f"{s:<6} {len(good):^6} {b['pose']:^10} {b['dC']:^10.2f} {b['dP']:^11.2f} {b['pose']:^14}")
print(); print("="*100); print("CONTACTS (<=3.5 A polar) FOR THE CLOSEST-C POSE OF EACH SEED"); print("="*100)
for s in seeds:
    b=min(allrows[s],key=lambda x:x["dC"]); p,rec=b["p"],b["rec"]
    print(f"--- seed {s}, pose {b['pose']}  (C-N5 {b['dC']:.2f} A) ---")
    cx=[]; px=[]
    for i in CARBOX_O: cx+=contacts(p[i][1],rec)
    for i in PHOS_O: px+=contacts(p[i][1],rec)
    nx=contacts(p[IN][1],rec)
    print("  carboxylate :", ", ".join(f"{n} {d:.2f}" for d,n in sorted(set(cx))) or "none")
    print("  phosphonate :", ", ".join(f"{n} {d:.2f}" for d,n in sorted(set(px))) or "NONE <=3.5 A")
    if not px:
        near=[]
        for i in PHOS_O: near+=contacts(p[i][1],rec,cutoff=5.5)
        print("                nearest >3.5:", ", ".join(f"{n} {d:.2f}" for d,n in sorted(set(near))[:6]))
    print("  amine N     :", ", ".join(f"{n} {d:.2f}" for d,n in sorted(set(nx))) or "none")
    print()
print("="*100); print("PAIRWISE HEAVY-ATOM RMSD (A) OF CLOSEST-C POSES"); print("="*100)
bests={s:min(allrows[s],key=lambda x:x["dC"]) for s in seeds}
print("       "+"".join(f"{s:>8}" for s in seeds))
for s1 in seeds:
    print(f"{s1:>6} "+"".join(f"{rmsd_heavy(bests[s1]['p'],bests[s2]['p']):8.2f}" for s2 in seeds))
if ref:
    print(); print("vs seed102 pose1 (reference):", "  ".join(f"s{s}={rmsd_heavy(bests[s]['p'],ref):.2f}" for s in seeds))
