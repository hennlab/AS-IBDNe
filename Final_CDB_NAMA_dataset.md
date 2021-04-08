### 1. Create new plink files with only the reference individuals.
```
plink --bfile nama_tgp_qc_pops --keep reference.keep --make-bed --out ref_samples
```
### 2. Run king on reference plink files to check for relatedness
```
module load king
king -b ref_samples.bed --ibdseg --cpus 15 --prefix ref_samples_king
king -b ref_samples.bed --ibs --cpus 15 --prefix ref_samples_king
king -b ref_samples.bed --related --cpus 15 --prefix ref_samples_king
king -b ref_samples.bed --unrelated --degree 2
king -b ref_samples.bed --unrelated --degree 3
```

### 3. Keep the unrelated individuals in the reference population

```
plink --bfile ref_samples --keep kingunrelated.txt --make-bed --out ref_unrelated
```


### 4. Combine with CDB data

```
plink --bfile ref_unrelated --bmerge /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NCTB_H3Africa/CDB/SNP-QC/cdb_80.bed /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NCTB_H3Africa/CDB/SNP-QC/cdb_80.bim /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NCTB_H3Africa/CDB/SNP-QC/cdb_80.fam --make-bed --out ref_CDB_merge

cut -f2 ref_CDB_merge.bim | wc -l # 3245933 overlapping SNPs

plink --bfile ref_CDB_merge --geno 0.05 --make-bed --out ref_CDB_merge_geno0.05
cut -f2 ref_CDB_merge_geno0.05.bim | wc -l
```
After removing relatives,
85 LWK, 88 GBR, 103 CHB and 26 NAMA

SA2002,SA2076,SA2088 need to be removed from Namas, but otherwise we will use all of the reference NAMA and 50 of each population in the reference set.

Taking first 50 of each population from the `kingunrelated.txt` file, and all the Namas minus SA2002,SA2076,SA2088.

```
grep "CHB" kingunrelated.txt | head -n 50 > balanced_ref_unrelated_smps.txt
grep "GBR" kingunrelated.txt | head -n 50 >> balanced_ref_unrelated_smps.txt
grep "LWK" kingunrelated.txt | head -n 50 >> balanced_ref_unrelated_smps.txt
grep "NAMA" kingunrelated.txt >> balanced_ref_unrelated_smps.txt
grep "NAMA" kingunrelated_toberemoved.txt >> balanced_ref_unrelated_smps.txt
# manually delete SA2002,SA2076,SA2088 with vim
```

List of reference samples to keep
```
plink --bfile ref_samples --keep balanced_ref_unrelated_smps.txt --make-bed --out ref_balanced_unrelated
```

Re-extract 1kg samples (150)

/share/genomes/1000genomes/phase3

```
for chr in `seq 1 22`; do plink --vcf /share/genomes/1000genomes/phase3/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr ${chr} --keep balanced_ref_unrelated_smps.txt --out 1kg_chr${chr} ; done
```

combine chromosomes
```
plink --bfile 1kg_chr1 --merge-list merge_chroms.txt --make-bed --out 1kg_all # tri-allelic errors
for chr in `seq 1 22` ; do plink --bfile 1kg_chr${chr} --exclude 1kg_all-merge.missnp --make-bed --out 1kg_chr${chr}_ex ; done
plink --bfile 1kg_chr1_ex --merge-list merge_chroms_ex.txt --make-bed --out 1kg_all # tri-allelic errors
```

merge with admixed CDB samples
```
plink --bfile 1kg_all --bmerge /share/hennlab/projects/IBDNe/as-ibdne/version_2.3/data/cdb_nama_0.05.admix.bed /share/hennlab/projects/IBDNe/as-ibdne/version_2.3/data/cdb_nama_0.05.admix.bim /share/hennlab/projects/IBDNe/as-ibdne/version_2.3/data/cdb_nama_0.05.admix.fam --make-bed --out 1kg_recombine
```

# remove indels
```
grep ";" 1kg_recombine.bim | cut -f2 > indels
plink --bfile 1kg_recombine --geno 0.05 --snps-only --exclude indels --make-bed --out 1kg_recombine_restart

78793121 variants removed due to missing genotype data (--geno).
535149 variants and 574 people pass filters and QC.
```
# add in nama reference individuals
```
plink --bfile /share/hennlab/projects/IBDNe/as-ibdne/dataset/NAMA_ref_only --bmerge /share/hennlab/projects/IBDNe/as-ibdne/version_2.3/data/1kg_recombine_restart.bed /share/hennlab/projects/IBDNe/as-ibdne/version_2.3/data/1kg_recombine_restart.bim /share/hennlab/projects/IBDNe/as-ibdne/version_2.3/data/1kg_recombine_restart.fam --make-bed --out 1kg_cdb_nama_final
NA19138
plink --bfile 1kg_cdb_nama_final --geno 0.05 --snps-only  --make-bed --out 1kg_cdb_nama_final_qc


414976 variants and 615 people pass filters and QC.
```
Running pipeline
```
nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile 1kg_recombine.yaml -j 30 -n -R msp_to_bed
```

# Verifying rfmix results match admixture results

```
cat data/admixed.txt | while read line ; do echo results/plots/bedfiles/1kg_cdb_nama_final_qc.${line}.0.BED results/plots/bedfiles/1kg_cdb_nama_final_qc.${line}.1.BED >> bedlist.txt ; done
```


```
python scripts/lai_global.py \
--bed_list bedlist.txt \
--ind_list data/admixed.txt \
--pops CHB,GBR,LWK,NAMA \
--out lai_global.txt
```

```
cat data/admixed.txt | while read line ; do grep $line lai_just_values ; grep $line admixture/4Q_with_ids ; done >> admixture_rfmix_compare
```

```
SA3169  0.0929  0.2938  0.0478  0.5655
SA3169  0.046198 0.082615 0.345124 0.526063
SA3202  0.1539  0.312   0.2196  0.3145
SA3202  0.238239 0.149250 0.338758 0.273753
SA3053  0.0523  0.2617  0.1109  0.5751
SA3053  0.104655 0.034265 0.308953 0.552127
SA3066  0.0824  0.3016  0.1362  0.4798
SA3066  0.134890 0.048822 0.354876 0.461411
SA3116  0.045   0.3809  0.0491  0.525
SA3116  0.053477 0.038554 0.413888 0.494081
SA3162  0.1038  0.2346  0.0867  0.575
SA3162  0.098579 0.084648 0.278129 0.538644
NC292   0.0162  0.0423  0.5091  0.4324
NC292   0.503549 0.007219 0.054356 0.434875
NC118   0.0909  0.2376  0.0594  0.6121
NC118   0.050769 0.072869 0.296741 0.579621
NC147   0.0436  0.1054  0.2464  0.6046
NC147   0.283985 0.045340 0.154925 0.515750
NC165   0.0132  0.2561  0.2225  0.5082
NC165   0.230303 0.003442 0.296643 0.469612
NC187   0.0328  0.1799  0.1275  0.6598
NC187   0.117568 0.022331 0.194933 0.665169

```
