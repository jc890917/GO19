exec(open("analyze_poses.py").read().split('for s in seeds:\n    print("="*100)')[0])
import statistics as st
cnn={}
for s in seeds:
    v=[]
    for ln in open(os.path.join(D,f"log_seed{s}.txt")):
        t=ln.split()
        if len(t)==5 and t[0].isdigit(): v.append((float(t[1]),float(t[3])))
    cnn[s]=v
print("RANK-1 POSE OF EVERY SEED  (--flexres A:257, corrected ligand)")
print(f"{'seed':<6}{'CNN':>8}{'aff':>7}{'C-N5':>7}{'C-C4X':>7}{'H-N5':>7}{'C-H..N5':>9}{'P-N5':>7}{'RMSDref':>9}{'Y257OH':>8}{'minClash':>10}  verdict")
good=[]
for s in seeds:
    x=allrows[s][0]; c,a=cnn[s][0][1],cnn[s][0][0]
    v="CATALYTIC" if x['dC']<4.0 else "decoy"
    if x['dC']<4.0: good.append(x)
    print(f"{s:<6}{c:>8.4f}{a:>7.2f}{x['dC']:>7.2f}{x['dC4']:>7.2f}{x['dH']:>7.2f}{x['aCHN']:>9.1f}{x['dP']:>7.2f}{x['rmsd']:>9.2f}{x['dOH']:>8.3f}{x['clash'][0]:>10.2f}  {v}")
print()
print(f"CNN>=0.70 : {sum(1 for s in seeds if cnn[s][0][1]>=0.70)}/10 seeds, ALL catalytic")
print(f"CNN<=0.60 : {sum(1 for s in seeds if cnn[s][0][1]<=0.60)}/10 seeds, NONE catalytic")
print(f"gap       : {min(cnn[s][0][1] for s in seeds if cnn[s][0][1]>=0.70):.4f} vs {max(cnn[s][0][1] for s in seeds if cnn[s][0][1]<=0.60):.4f}")
print()
def ms(v): return f"{st.mean(v):.2f} +/- {st.stdev(v):.2f}"
print("CONSENSUS ENSEMBLE (the 6 catalytic rank-1 poses):")
print("  C(react)-N5   ", ms([x['dC'] for x in good]), " A   range %.2f-%.2f"%(min(x['dC'] for x in good),max(x['dC'] for x in good)))
print("  C(react)-C4X  ", ms([x['dC4'] for x in good]), " A")
print("  H...N5        ", ms([x['dH'] for x in good]), " A  (H rebuilt tetrahedrally)")
print("  C-H...N5 angle", ms([x['aCHN'] for x in good]), " deg")
print("  P-N5          ", ms([x['dP'] for x in good]), " A")
print("  Y257 OH shift ", ms([x['dOH'] for x in good]), " A")
r=[rmsd_heavy(good[i]['p'],good[j]['p']) for i in range(len(good)) for j in range(i+1,len(good))]
print("  pairwise RMSD ", ms(r), " A   range %.2f-%.2f"%(min(r),max(r)))
print()
print("PHOSPHONATE PARTNER IN EACH CATALYTIC RANK-1 POSE (<=3.5 A, polar)")
for s in seeds:
    x=allrows[s][0]
    if x['dC']>=4.0: continue
    hits=set()
    for i in PHOS_O:
        for resn,resi,name,el,c in x['rec']:
            if el in "NO" and dist(x['p'][i][1],c)<=3.5: hits.add(f"{resn}{resi}:{name} {dist(x['p'][i][1],c):.2f}")
    print(f"  seed{s:<5}"+(", ".join(sorted(hits)) if hits else "NONE"))
print()
print("R53 NH1/NH2/NE distance to nearest phosphonate O, catalytic rank-1 poses")
for s in seeds:
    x=allrows[s][0]
    if x['dC']>=4.0: continue
    best=min((dist(x['p'][i][1],c),name) for i in PHOS_O for resn,resi,name,el,c in x['rec'] if resi==53 and name in("NE","NH1","NH2"))
    print(f"  seed{s:<5} ARG53:{best[1]} {best[0]:.2f} A")
print()
print("DECOY SEEDS - where the rank-1 pose actually binds")
for s in seeds:
    x=allrows[s][0]
    if x['dC']<4.0: continue
    hits=set()
    for i in [ICOO]+list(CARBOX_O)+list(PHOS_O)+[IP]:
        for resn,resi,name,el,c in x['rec']:
            if el in "NO" and dist(x['p'][i][1],c)<=3.5: hits.add(f"{resn}{resi}")
    print(f"  seed{s:<5} C-N5 {x['dC']:.2f} A  CNN {cnn[s][0][1]:.4f}  contacts: "+",".join(sorted(hits)))
