### 1. Running King to check for relatedness between individuals

```bash
module load king
king -b nama_tgp_qc.bed --ibdseg --cpus 15 --prefix nama_tgp_qc

king -b nama_tgp_qc.bed --ibs --cpus 15 --prefix nama_tgp_qc

KING starts at Tue Dec 29 14:39:18 2020
Loading genotype data in PLINK binary format...
Read in PLINK fam file nama_tgp_qc.fam...
  PLINK pedigrees loaded: 395 samples
Read in PLINK bim file nama_tgp_qc.bim...

  Genotype data consist of 1721935 autosome SNPs
  PLINK maps loaded: 1721935 SNPs
Read in PLINK bed file nama_tgp_qc.bed...
  PLINK binary genotypes loaded.
  KING format genotype data successfully converted.
Autosome genotypes stored in 26906 words for each of 395 individuals.

Options in effect:
        --ibs
        --cpus 15
        --prefix nama_tgp_qc

Total length of chromosomal segments usable for IBD segment analysis is 2668.7 MB
  Information of these chromosomal segments can be found in file nama_tgp_qcallsegs.txt

Within-family IBS data saved in file nama_tgp_qc.ibs

Relationship summary (total relatives: 0 by pedigree, 798 by inference)
  Source        MZ      PO      FS      2nd     3rd     OTHER
  ===========================================================
  Pedigree      0       0       0       0       0       5151
  Inference     0       33      3       366     396     4353

IBS and relationship inference across families starts at Tue Dec 29 14:39:28 2020
15 CPU cores are used.
                                         ends at Tue Dec 29 14:39:31 2020
Between-family IBS data saved in file nama_tgp_qc.ibs0
```

To solve the issue with the chromosomes all looking like the same color because of relatedness in the namas:

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

```
Family clustering starts at Tue Jan  5 11:37:51 2021
Autosome genotypes stored in 26906 words for each of 337 individuals.
Sorting autosomes...
Total length of chromosomal segments usable for IBD segment analysis is 2668.7 MB
  Information of these chromosomal segments can be found in file kingallsegs.txt

32 CPU cores are used to compute the pairwise kinship coefficients...
Clustering up to 2nd-degree relatives in families...
Individual IDs are unique across all families.
No families were found to be connected.

A list of 302 unrelated individuals saved in file kingunrelated.txt
An alternative list of 35 to-be-removed individuals saved in file kingunrelated_toberemoved.txt

Extracting a subset of unrelated individuals ends at Tue Jan  5 11:38:02 2021
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

Keep new balanced reference samples
```
plink --bfile ref_samples --keep balanced_ref_unrelated_smps.txt --make-bed --out ref_balanced_unrelated
```

Combine with CDB data
```
# first add population to FID column
awk '{print "CDB"" "$2" "$3" "$4" "$5" "$6}' /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NCTB_H3Africa/CDB/SNP-QC/cdb_80.fam > cdb_80_pops.fam

plink --bfile ref_balanced_unrelated --bmerge cdb_80_pops.bed cdb_80_pops.bim cdb_80_pops.fam --make-bed --out ref_bal_CDB_merge
```

Check overlapping SNPs
```
cut -f2 ref_bal_CDB_merge.bim | wc -l # 3245933 SNPs
plink --bfile ref_bal_CDB_merge --geno 0.05 --make-bed --out ref_bal_CDB_merge_geno0.05 # 435582 SNPs
```

Check PCA
```
plink --bfile ref_bal_CDB_merge_geno0.05 --pca 4 --out ref_bal_CDB_merge_geno0.05
```
PCA looks as expected. proceeding with running as-ibdne.


### 5. Run as-ibdne on the new dataset

```
sed 's/\t/ /g' balanced_ref_unrelated_smps.txt > balanced_ref_unrelated_smps.ref.keep
cut -f2 cdb_80_pops.fam > cdb_admix.keep
```

config file
```yaml
gmap: /share/hennlab/reference/recombination_maps/genetic_map_HapMapII_GRCh37/

dataset: ref_bal_CDB_merge_geno0.05

ref: /share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3

ids_ref: data/balanced_ref_unrelated_smps.ref.keep

ids_admix: data/cdb_admix.keep

pops: data/cdb_merge_pops.txt
```

Run snakemake

```bash
/share/hennlab/progs/miniconda3/bin/snakemake --configfile config/config_ref_CDB_merge.yaml -j 20
```

/share/hennlab/progs/miniconda3/bin/snakemake --configfile config/config_ref_CDB_merge.yaml -j 20 --rulegraph | dot -Tpng > rulegraph.png

Check on karyograms
```
python scripts/collapse_ancestry.py \
--rfmix results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.2.Viterbi.txt \
--snp_locations results/RFMIX/ref_bal_CDB_merge_geno0.05_chr1.snp_locations \
--fbk results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.2.ForwardBackward.txt \
--fbk_threshold 0.9 \
--ind SA3097 \
--ind_info results/RFMIX/ref_bal_CDB_merge_geno0.05.sample \
--pop_labels GBR,CHB,LWK,NAMA \
--out SA3097
```

python scripts/plot_karyogram.py \
--bed_a SA3097_A.bed \
--bed_b SA3097_B.bed \
--ind SA3097 \
--pop_order GBR,CHB,LWK,NAMA \
--out SA3097.png

## Plot 6 karyograms for CDB and send to Austin & Brenna

SA3107, SA3140, SA3154, SA3124, SA3089, SA3097

SA3060 SA3188 SA3164

### Testing step 7 script

## implementing new rfmix version 2.3

./rfmix -f <query VCF/BCF file>
	-r <reference VCF/BCF file>
	-m <sample map file>
	-g <genetic map file>
	-o <output basename>
	--chromosome=<chromosome to analyze>


python scripts/RunRFMix.py -e 2 -w 0.2 --num-threads {threads} --use-reference-panels-in-EM --forward-backward PopPhased {input.all} {input.cl} {input.snp} -o {params}


Trying to figure out why RFmix was rerun and what changed between the good runs and bad runs:

cut -d" " -f2 balanced_ref_unrelated_smps.ref.keep > reference.samples.keep
cut -d" " -f2 cdb_admix.keep > admix.samples.keep
python ../scripts/shapeit2rfmix.py --shapeit_hap_ref ref_bal_CDB_merge_geno0.05.chr10.phased.haps --shapeit_hap_admixed ref_bal_CDB_merge_geno0.05.chr10.phased.haps --shapeit_sample_ref ref_bal_CDB_merge_geno0.05.chr10.phased.sample --shapeit_sample_admixed ref_bal_CDB_merge_geno0.05.chr10.phased.sample --ref_keep reference.samples.keep  --admixed_keep admix.samples.keep --chr 10 --genetic_map /share/hennlab/reference/recombination_maps/genetic_map_HapMapII_GRCh37/chr10.gmap.txt --out script_out


Looks like the list of ref pops included CDB
problem run fix classes rule inputs (less 2021-01-19T160309.758281.snakemake.log) : results/RFMIX/ref_bal_CDB_merge_geno0.05.GBR.keep, results/RFMIX/ref_bal_CDB_merge_geno0.05.CHB.keep, results/RFMIX/ref_bal_CDB_merge_geno0.05.LWK.keep, results/RFMIX/ref_bal_CDB_merge_geno0.05.NAMA.keep, results/RFMIX/ref_bal_CDB_merge_geno0.05.
CDB.keep, results/RFMIX/ref_bal_CDB_merge_geno0.05.sample

fix classes rule inputs for the run that worked: results/RFMIX/nama_tgp_qc_pops.GBR.keep, results/RFMIX/nama_tgp_qc_pops.CHB.keep, results/RFMIX/nama_tgp_qc_pops.LWK.keep, results/RFMIX/nama_tgp_qc_pops.NAMA.keep, results/RFMIX/nama_tgp_qc_pops.sample

python ../scripts/classes.py --ref ref_bal_CDB_merge_geno0.05.GBR.keep,ref_bal_CDB_merge_geno0.05.CHB.keep,ref_bal_CDB_merge_geno0.05.LWK.keep,ref_bal_CDB_merge_geno0.05.NAMA.keep --sample ref_bal_CDB_merge_geno0.05.sample --out new_fix.classes


When the good version of rfmix was run, it used this file results/RFMIX/new.classes , which looks like
```
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 0 0 0 0 0 0 4 4 4 4 4 4 0 0 0 0 0 0 4 4 4 4 0 0 0 0 4 4 0 0 4 4 0 0 4 4 0 0 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 4 4 4 4 4 4 0 0 0 0 0 0 4 4 4 4 0 0 0 0 4 4 4 4 4 4 0 0 4 4 0 0 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 4 4 0 0 0 0 4 4 0 0 4 4 4 4 0 0 0 0 0 0 4 4 0 0 0 0 0 0 4 4 0 0 4 4 0 0 0 0 0 0 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

while the file output by the classes.py script looks like this
```
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```
the old new.classes file was leftover from the previous run where  

python ../scripts/classes.py --ref ref_bal_CDB_merge_geno0.05.GBR.keep,ref_bal_CDB_merge_geno0.05.CHB.keep,ref_bal_CDB_merge_geno0.05.LWK.keep,ref_bal_CDB_merge_geno0.05.NAMA.keep --sample new_tail.sample --out new_tail_fix.classes


Testing running rfmix on dataset with classes file size 2

rule run_rfmix:
    input: results/RFMIX/ref_bal_CDB_merge_geno0.05_chr1.alleles, results/RFMIX/ref_bal_CDB_merge_geno0.05_chr1.snp_locations, results/RFMIX/new.classes
    output: results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.0.ForwardBackward.txt, results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.1.ForwardBackward.txt, results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.2.ForwardBackward.txt, results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.0.Viterbi.txt, results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.1.Viterbi.txt, results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.2.Viterbi.txt, results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.allelesRephased0.txt, results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.allelesRephased1.txt, results/RFMIX/ref_bal_CDB_merge_geno0.05.chr1.rfmix.allelesRephased2.txt
    jobid: 1
    wildcards: dataset=ref_bal_CDB_merge_geno0.05, chrnum=1
    threads: 10

running in snakemake screen

python RunRFMix.py -e 2 -w 0.2 --num-threads 20 --use-reference-panels-in-EM --forward-backward PopPhased results/RFMIX/ref_bal_CDB_merge_geno0.05_chr1.alleles debugging/small.classes results/RFMIX/ref_bal_CDB_merge_geno0.05_chr1.snp_locations -o debugging/small_classes_file_test

Also, run admixture plot for our data with K =4
