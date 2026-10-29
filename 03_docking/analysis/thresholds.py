exec(open("analyze_poses.py").read().split('for s in seeds:\n    print("="*100)')[0])
import json
cnn={}
for s in seeds:
    vals=[]
    for ln in open(os.path.join(D,f"log_seed{s}.txt")):
        t=ln.split()
        if len(t)==5 and t[0].isdigit(): vals.append(float(t[3]))
    cnn[s]=vals
print("THRESHOLD SCAN - number of poses (out of 20) with reacting C within X A of FAD N5")
print(f"{'seed':<6}{'<4.0':>6}{'<4.5':>6}{'<5.0':>6}{'<5.5':>6}{'<6.0':>6}{'best':>8}{'rank':>6}{'CNN':>8}{'RMSDref':>9}")
for s in seeds:
    r=allrows[s]; b=min(r,key=lambda x:x['dC'])
    print(f"{s:<6}"+"".join(f"{sum(1 for x in r if x['dC']<t):>6}" for t in (4.0,4.5,5.0,5.5,6.0))
          +f"{b['dC']:>8.2f}{b['pose']:>6}{cnn[s][b['pose']-1]:>8.4f}{b['rmsd']:>9.2f}")
print()
print("CATALYTIC-ORIENTATION POSES (carboxylate O within 3.5 A of R272/R314/Y257 side chain)")
print(f"{'seed':<6}{'pose':>5}{'C-N5':>7}{'CNN':>8}{'RMSDref':>9}  anchors")
ANCH={(272,'NE'),(272,'NH1'),(272,'NH2'),(314,'NE'),(314,'NH1'),(314,'NH2'),(257,'OH')}
for s in seeds:
    for x in allrows[s]:
        hits=set()
        for i in CARBOX_O:
            for resn,resi,name,el,c in x['rec']:
                if (resi,name) in ANCH and dist(x['p'][i][1],c)<=3.5: hits.add(f"{resn}{resi}:{name}")
        if hits:
            print(f"{s:<6}{x['pose']:>5}{x['dC']:>7.2f}{cnn[s][x['pose']-1]:>8.4f}{x['rmsd']:>9.2f}  "+",".join(sorted(hits)))
print()
print("WHERE DOES THE PHOSPHONATE SIT? (residues within 3.5 A of any phosphonate O, rank-1 poses)")
for s in seeds:
    x=allrows[s][0]; hits=set()
    for i in PHOS_O:
        for resn,resi,name,el,c in x['rec']:
            if el in "NO" and dist(x['p'][i][1],c)<=3.5: hits.add(f"{resn}{resi}")
    print(f"  seed{s} pose1 (CNN {cnn[s][0]:.4f}, C-N5 {x['dC']:.2f}): "+(",".join(sorted(hits)) or "none"))
