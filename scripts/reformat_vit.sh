# usage: bash reformat_vit.sh americans_subset_remove.chr18.msp.tsv americans_subset_remove.chr18.bim chr18_new.viterbi.tsv 18
MSP=$1
BIM=$2
VIT=$3
CHR=$4

grep '#' ${MSP} | tail -n 1 | cut -f7- | tr '\t' '\n' | cut -f1 -d'.' | tr '\n' ' ' | sed 's/^/POS /' > chr${CHR}.vit_header
echo >> chr${CHR}.vit_header
cut -f4 ${BIM} > chr${CHR}.bim_positions
paste -d' ' chr${CHR}.bim_positions ${VIT} > chr${CHR}.vit.bimpos
cat chr${CHR}.vit_header chr${CHR}.vit.bimpos > chr${CHR}.vit.fixed
mv chr${CHR}.vit.fixed ${VIT}.fixed
