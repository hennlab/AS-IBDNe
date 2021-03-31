# usage: bash reformat_vit.sh americans_subset_remove.chr18.msp.tsv results/phasing/americans_subset_remove.chr18.bim results/RFmix/americans_subset_remove.chr18.vit.tsv 18
MSP=$1
BIM=$2
VIT=$3
CHR=$4

grep '#' ${MSP} | tail -n 1 | cut -f7- | tr '\t' '\n' | cut -f1 -d'.' | tr '\n' ' ' | sed 's/^/POS /' > results/RFmix/chr${CHR}.vit_header
cut -f4 ${BIM} > results/RFmix/chr${CHR}.bim_positions ; paste -d' ' results/RFmix/chr${CHR}.bim_positions ${VIT} > chr${CHR}.newfile ; mv chr${CHR}.newfile ${VIT}
cat results/RFmix/chr${CHR}.vit_header ${VIT} > chr${CHR}.newvit ; mv chr${CHR}.newvit ${VIT}


#cut -f4 results/phasing/americans_subset_remove.chr18.bim > results/RFmix/chr18.bim_positions ; paste -d' ' results/RFmix/bim_positions results/RFmix/americans_subset_remove.chr18.vit.tsv > chr18.newfile ; mv chr18.newfile results/RFmix/americans_subset_remove.chr18.vit.tsv


#python scripts/filter_gapfilled_ibd_ancestry.py results/IBD-segs/americans_subset_remove.chr18.phased.ibd.gz results/IBD-segs/americans_subset_remove.chr18.phased.filled.ibd chr18_new_vit_fixed.tsv 3 > {output}
