DATASET=$1
ANC=$2
ADMIDS=$3

# combine ancestries - old code from Brownings script compatible with HapIBD
# for chr in `seq 1 22`; do cat results/IBD-segs/${DATASET}.chr${chr}.phased.filled.allanc.ibd | grep -i "[[:space:]]${ANC}"'$' | cut -f1-8 -d' ' > results/IBD-segs/${DATASET}.chr${chr}.phased.filled.anc${ANC}.ibd ; done

# fixed code - compatible with refined-ibd
for chr in `seq 1 22`; do cat results/IBD-segs/${DATASET}.chr${chr}.phased.filled.allanc.ibd | grep -i "[[:space:]]${ANC}"'$' | cut -f1-7,9 -d' ' > results/IBD-segs/${DATASET}.chr${chr}.phased.filled.anc${ANC}.ibd ; done

# combine chroms
cat results/IBD-segs/${DATASET}.chr*.phased.filled.anc${ANC}.ibd > results/IBD-segs/${DATASET}.phased.filled.anc${ANC}.ibd

# calculate adjusted number of pairs of haplotypes
cat results/RFmix/${DATASET}.chr*.vit.tsv | java -jar progs/filtercolumns.jar 1 ${ADMIDS} | python scripts/adjust_npairs.py ${ANC} > results/IBDne/${DATASET}.anc${ANC}_npairs
