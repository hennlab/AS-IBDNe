## Author: Shyamalika Gopalan

INPUT=$1
MAP=$2
LOG=$3
OUTPUT=$4
REF_HAP=$5
REF_LEG=$6
REF_SAMP=$7
THREADS=$8

# ShapeIT phasing and conversion to plink files
#----------------------------------------------
# Divide plink file by chromosome, trim SNPs that lie outside the boundaries of the genetic map
# Divide plink file by chromosome, trim SNPs that lie outside the boundaries of the genetic map

Rscript scripts/trim_bimfile.R ${INPUT} ${MAP}
plink --bfile ${INPUT} --exclude ${INPUT}_removeSNPs --make-bed --out ${INPUT}_endsTrimmed
#cat ${FILE}.chr${chr}_endsTrimmed.bim >> dat_QCed_endsTrimmed.bim

# Run SHAPEIT
shapeit -B ${INPUT}_endsTrimmed -M ${MAP} --input-ref ${REF_HAP} ${REF_LEG} ${REF_SAMP} --duohmm -W 5 -O ${OUTPUT} --output-log ${LOG}.log -T ${THREADS}

files=( ${LOG}.snp.strand )
if (( ${#files[@]} )); then
  N=`ls -l ${LOG}.snp.strand | wc -l`
  while [ $N -gt 0 ]
  do
      grep Missing ${LOG}.snp.strand | cut -f4 > ${INPUT}_missing_from_ref.txt
      grep Strand ${LOG}.snp.strand | awk '$11==1 {print $4}' > ${INPUT}_flip_in_dataset.txt
      grep Strand ${LOG}.snp.strand | awk '$11==0 {print $4}' > ${INPUT}_cannot_flip_in_dataset.txt
      cat ${INPUT}_missing_from_ref.txt ${INPUT}_cannot_flip_in_dataset.txt > ${INPUT}_remove_from_dataset.txt
      rm -f ${LOG}.snp.strand*
      rm -f ${INPUT}_missing_from_ref.txt ${INPUT}_cannot_flip_in_dataset.txt
      plink --bfile ${INPUT}_endsTrimmed --flip ${INPUT}_flip_in_dataset.txt --exclude ${INPUT}_remove_from_dataset.txt --make-bed --out ${INPUT}_endsTrimmed_misalignRepaired
      shapeit -B ${INPUT}_endsTrimmed_misalignRepaired -M ${MAP} --input-ref ${REF_HAP} ${REF_LEG} ${REF_SAMP} --duohmm -W 5 -O ${OUTPUT} --output-log ${LOG}.log -T ${THREADS}
      N=`ls -l ${LOG}.snp.strand | wc -l`
  done
fi
