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

# add in nama reference individuals

plink --bfile /share/hennlab/projects/IBDNe/as-ibdne/dataset/NAMA_ref_only --bmerge /share/hennlab/projects/IBDNe/as-ibdne/version_2.3/data/1kg_recombine_restart.bed /share/hennlab/projects/IBDNe/as-ibdne/version_2.3/data/1kg_recombine_restart.bim /share/hennlab/projects/IBDNe/as-ibdne/version_2.3/data/1kg_recombine_restart.fam --make-bed --out 1kg_cdb_nama_final
NA19138
plink --bfile 1kg_cdb_nama_final --geno 0.05 --snps-only  --make-bed --out 1kg_cdb_nama_final_qc


414976 variants and 615 people pass filters and QC.

nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile 1kg_recombine.yaml -j 30 -n
