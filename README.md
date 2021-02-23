# AS-IBDNe
Snakemake pipeline for running ancestry-specific IBDNe


1. Set Up:

Creating conda environment:
conda install -c bioconda rfmix

```bash
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
export PATH=/share/hennlab/progs/RFMix/RFMix_v1.5.4/:$PATH
ln -s /share/hennlab/progs/RFMix
conda activate IBDne-env
```


2. QC beforehand
```bash
plink --bfile nama_tgp --geno 0.05 --mind 0.1 --make-bed --out nama_tgp_qc
# with populations in family ID column --> nama_tgp_qc_pops
```

3. Set up config file
The config file will contain all the file paths that change with every run of the snakemake pipeline.

Example file:
```bash
gmap: /share/hennlab/reference/recombination_maps/genetic_map_HapMapII_GRCh37/

dataset: nama_tgp_qc_pops

ref: /share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3

ids_ref: data/reference.keep

ids_admix: data/admix.keep

pops: data/population.ids.txt
```

ids_ref and ids_admix must be a file with a single column - the list of sample IDs

4. Set up working directory
data - contains input bed/bim/fam files. must be named according to {dataset}.bed /.bim /.fam
  must include both reference and admixed individuals

5. Run Snakemake

```bash
nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile config/config.yaml -j 20
nice /share/hennlab/progs/miniconda3/bin/snakemake -R rfmix_input --configfile config/config.yaml -j 20
nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile config_ref_CDB_merge.yaml -j 20
```

6. Plot karyograms after running RFMix

```bash
python scripts/collapse_ancestry.py \
--rfmix results/RFMIX/nama_tgp_qc_pops.chr1.rfmix.2.Viterbi.txt \
--snp_locations results/RFMIX/nama_tgp_qc_pops_chr1.snp_locations \
--fbk results/RFMIX/nama_tgp_qc_pops.chr1.rfmix.2.ForwardBackward.txt \
--fbk_threshold 0.9 \
--ind SA2046 \
--ind_info results/RFMIX/nama_tgp_qc_pops.sample \
--pop_labels GBR,CHB,LWK,NAMA \
--out SA2046

python scripts/collapse_ancestry.py \
--rfmix results/RFMIX/nama_tgp_qc_pops.chr1.rfmix.2.Viterbi.txt \
--snp_locations results/RFMIX/nama_tgp_qc_pops_chr1.snp_locations \
--fbk results/RFMIX/nama_tgp_qc_pops.chr1.rfmix.2.ForwardBackward.txt \
--fbk_threshold 0.9 \
--ind SAK050 \
--ind_info results/RFMIX/nama_tgp_qc_pops.sample \
--pop_labels GBR,CHB,LWK,NAMA \
--out SAK050

python scripts/collapse_ancestry.py \
--rfmix results/RFMIX/nama_tgp_qc_pops.chr1.rfmix.2.Viterbi.txt \
--snp_locations results/RFMIX/nama_tgp_qc_pops_chr1.snp_locations \
--fbk results/RFMIX/nama_tgp_qc_pops.chr1.rfmix.2.ForwardBackward.txt \
--fbk_threshold 0.9 \
--ind SA2106 \
--ind_info results/RFMIX/nama_tgp_qc_pops.sample \
--pop_labels GBR,CHB,LWK,NAMA \
--out SA2106

python scripts/collapse_ancestry.py \
--rfmix results/RFMIX/nama_tgp_qc_pops.chr1.rfmix.2.Viterbi.txt \
--snp_locations results/RFMIX/nama_tgp_qc_pops_chr1.snp_locations \
--fbk results/RFMIX/nama_tgp_qc_pops.chr1.rfmix.2.ForwardBackward.txt \
--fbk_threshold 0.9 \
--ind SA2085 \
--ind_info results/RFMIX/nama_tgp_qc_pops.sample \
--pop_labels GBR,CHB,LWK,NAMA \
--out SA2085


python scripts/plot_karyogram.py \
--bed_a SA3107_A.bed \
--bed_b SA3107_B.bed \
--ind SA3107 \
--pop_order GBR,CHB,LWK,NAMA \
--out SA3107.png

python scripts/plot_karyogram.py \
--bed_a SA2046_A.bed \
--bed_b SA2046_B.bed \
--ind SA2106 \
--pop_order GBR,CHB,LWK,NAMA \
--out SA2106.png


```

## Getting versions 1.7 running:

in version 1.5: --use-reference-panels-in-EM flag (not using this flag means reference panels discarded after initial inference step) isn't in version 1.7

Version 1.7 does have --include-reference (tree building includes reference haplotypes in tip node counts, not just admixed. (see manual))

```
nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile config_ref_CDB_merge.yaml -j 20
```


Deciphering the output files of three different versions:

RFMIX 1.5.4:
- ForwardBackward.txt
    - one row per SNP and one column per ancestry.
    - Value in each column is posteriorprobability of that ancestryt that SP n that haplotype.
- Viterbi.txt
    - One rowper SNP and one columner amixed alotype
- allelesRephased0.txt
    - one rowper SNP and one column peraplotype. Values are 0or1, wher each SNP has ad t allesconverte to binry frmat

RFMMIX 1.7:
- fb-probs.tsv
- fb-sis.tsv
- rfmix.Q
- viterbi-msp.tsv

RFMix2.0:
- msp.tsv: the most likely assignment of subpopulations per CRF point (Viterbi)
    - rows corresponding to genomic position and columns corresponding to haplotypes
    - condensed such that CRF windows are combined if all query samples are in the sample subpopulations for successive windows. Thus, each line might represent several CRF points.
- fb.tsv: marginal probabilities of each subpopulation being the ancestral population of the corresponding CRF point (FwBw)
    - rows corresponding to genomic position and columns corresponding to haplotypes
    - haplotypes are tab delimited, but the array of probabilities for each haplotype at each window (row) is a set of space delimited columns within each tab delimited haplotype   column
- rfmix.Q: Global diploid ancestry estimates
