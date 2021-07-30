#for each ibd segment, output it with probability proportional to the proportion of it that is assigned to a given ancestry
# reads in both gapfilled and non-gapfilled and uses the hap indices from the non-merged part, but the length from the gapfilled
# adds a column with the ancestry assignment, then can filter on this and cut it off to get the ancestry-specific files

# usage: python filter_gapfilled_ibd_ancestry.py nongapfilledibd gapfilledibd viterbifile nancestries
import sys,random,gzip


if sys.argv[1][-3:]==".gz":
    nonmergedfile = gzip.open(sys.argv[1])
else:
    nonmergedfile = open(sys.argv[1])
if sys.argv[2][-3:]==".gz":
    mergedfile = gzip.open(sys.argv[2])
else:
    mergedfile = open(sys.argv[2])
if sys.argv[3][-3:]==".gz":
    viterbifile = gzip.open(sys.argv[3])
else:
    viterbifile = open(sys.argv[3])
nanc = int(sys.argv[4])

pos = []
vt = []
ids = viterbifile.readline().split()[1::2]
for line in viterbifile:
    bits = line.split()
    pos.append(int(bits[0]))
    vt.append([int(x)-1 for x in bits[1:]])

ibd = dict()
# read in the nonmerged data
for line in nonmergedfile:
    bits = line.split()
    id1 = bits[0]; id2 = bits[2]
    if id1 not in ids or id2 not in ids: continue
    hap1 = int(bits[1]); hap2 = int(bits[3])
    pos1 = int(bits[5]); pos2 = int(bits[6])
    if hap1 == 0 or hap2 == 0: continue
    if id1 not in ibd:
        ibd[id1] = dict()
    if id2 not in ibd[id1]: 
        ibd[id1][id2] = []
    ibd[id1][id2].append((pos1,pos2,hap1,hap2))

for line in mergedfile:
    bits = line.split()
    id1 = bits[0]; id2 = bits[2]
    if id1 not in ids or id2 not in ids: continue
    pos1 = int(bits[5]); pos2 = int(bits[6])
    acount1 = [0]*nanc; acount2 = [0]*nanc; n = 0
    for segment in ibd[id1][id2]:
        p1,p2,hap1,hap2 = segment
        if pos1 > p1 or pos2 < p2: continue
        index1 = (ids.index(id1)*2+hap1-1)
        index2 = (ids.index(id2)*2+hap2-1)
        oldpos = p1
        for i,x in enumerate(pos):
            if x>=p1 and x<=p2:
   
                a1 = vt[i][index1]
                a2 = vt[i][index2]
                acount1[a1] += x-oldpos
                acount2[a2] += x-oldpos
                n += x - oldpos
                oldpos = x
    aprop1 = [x/float(n) for x in acount1]
    aprop2 = [x/float(n) for x in acount2]
    aprops = [(x+y)/2.0 for x,y in zip(aprop1,aprop2)]
    if sum(aprops)<0.99 or sum(aprops)>1.01:
        print >> sys.stderr, "code error",aprops
        raise SystemExit
    # choose an ancestry, randomly according to aprop
    myrand = random.random()
    cumsum = 0
    for anc in range(nanc):
        cumsum += aprops[anc]
        if cumsum > myrand:
            thisanc = anc
            break
    for x in bits:
        print x,
    print thisanc+1
