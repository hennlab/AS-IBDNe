# AS-IBDNe
Snakemake pipeline for running ancestry-specific IBDNe

1. Set Up Environment :

Required packages:
 - rfmix version 2.3
 - shapeit

Required scripts in `scripts` folder:
-

```bash
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate IBDne-env
module load bcftools
module load samtools
module load htslib
```

2. QC beforehand
```bash
plink --bfile cdb_nama.ref --geno 0.05 --make-bed --out cdb_nama_0.05.ref
plink --bfile cdb_nama.admix --geno 0.05 --make-bed --out cdb_nama_0.05.admix
```

3. Required files and folder structure
- data
   - {dataset}.bim /bed /fam # plink dataset containing reference and admixed individuals
- config.yaml
- genetic map file : "The genetic map file is tab delimited text containing at least 3 columns. The first 3 columns are intepreted as chromosome, physical position in bp, genetic position in cM. Any number of columns or other information may follow, it is ignored. The chromosome column is a string token (which may be an string of digits) that must match those used in the VCF/BCF inputs. The genetic map file should contain the map for the entire genome (all chromosomes). Blank lines and lines beginning with a '#' are ignored."
- sample map file: tab delimited text with at least two columns. The first column gives the sample name or identifier, which must match the one used in the reference VCF/BCF. The second column is a string naming a subpopulation and may contain spaces (e.g., "European", or "East African")
- Snakefile

config file:
- dataset: {dataset}.bim /bed/fam
- genmap: genetic map file for all chromosomes
- smpmap: sample map file
- gmap: location of genetic map files for individual chromosomes.
- ref: path to reference haps/sample/legend files, including prefix
- admix_samples: text file containing list of individuals in dataset that are admixed, no ref individuals (1 column only)


4. Set up config file

```yaml
dataset: ref_bal_CDB_merge_geno0.05

genmap: /share/hennlab/reference/recombination_maps/genetic_map_AfricanAmerian/AAmap_rfmix_sort

smpmap: data/ref_samplefile.txt

admix_samples: data/admixed.txt

ref: /share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3

gmap: /share/hennlab/reference/recombination_maps/genetic_map_HapMapII_GRCh37/
```

5. Run pipeline

```
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate IBDne-env
nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile ref_balanced.config.yaml -j 20 -n
nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile cdb_nama_sepphase.yaml -j 20 -n
nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile cdb_nama_0.05.config.yaml -j 20 -R run_rfmix -n
```

# need to go back and sort chr 13 in the genetic map file

6. Plotting rfmix results - scripts from Gerald Van der Eeden

1. Convert msp to bed format

```
Order of populations in this script: whatever is defined by rfmix in the top of the msp file eg CHB=0       GBR=1   LWK=2   NAMA=3
Rscript msp_to_bed.R results/RFmix/msp_files results/plots/bedfiles/ 20 CHB GBR LWK NAMA
```

2. Use alicia martins karyogram script



python plot_karyogram.py \
--bed_a results/plots/bedfiles/SA3042.0.BED \
--bed_b results/plots/bedfiles/SA3042.1.BED \
--ind SA3042 \
--pop_order CHB,GBR,LWK,NAMA \
--out results/plots/simple_SA3042.png


python plot_karyogram.py \
--bed_a results/plots/bedfiles/SA3103.0.BED \
--bed_b results/plots/bedfiles/SA3103.1.BED \
--ind SA3103 \
--pop_order CHB,GBR,LWK,NAMA \
--out results/plots/SA3103.png

python plot_karyogram.py \
--bed_a results/plots/bedfiles/SA3180.0.BED \
--bed_b results/plots/bedfiles/SA3180.1.BED \
--ind SA3180 \
--pop_order CHB,GBR,LWK,NAMA \
--out results/plots/SA3180.png

python plot_karyogram.py \
--bed_a results/plots/bedfiles/SA3071.0.BED \
--bed_b results/plots/bedfiles/SA3071.1.BED \
--ind SA3071 \
--pop_order CHB,GBR,LWK,NAMA \
--out results/plots/SA3071.png

python scripts/plot_karyogram.py \
--bed_a results/plots/bedfiles_newmap/SA3126.0.BED \
--bed_b results/plots/bedfiles_newmap/SA3126.1.BED \
--ind SA3126 \
--pop_order CHB,GBR,LWK,NAMA \
--out results/plots/newmap_SA3126.png

SA3071

python scripts/plot_karyogram.py \
--bed_a results/plots/bedfiles_test3/SA3159.0.BED \
--bed_b results/plots/bedfiles_test3/SA3159.1.BED \
--ind SA3159 \
--pop_order CHB,GBR,LWK,NAMA \
--out results/plots/test3_SA3159.png



python scripts/plot_karyogram.py \
--bed_a results/plots/bedfiles/cdb_nama_0.05SA3146.0.BED \
--bed_b results/plots/bedfiles/cdb_nama_0.05SA3146.1.BED \
--ind SA3146 \
--pop_order GBR,CHB,LWK,NAMA \
--out results/plots/SA3146.png


- put in one plink file
- put in one text file for each set
- use vcf tools to separate vcf files
- bcf is binary version of vcf
- use bcftools or vcftools to pull out individuals not plink because plink messes up phasing






python scripts/plot_karyogram.py \
--bed_a results/plots/bedfiles_phasedtogetherSA3042.0.BED \
--bed_b results/plots/bedfiles_phasedtogetherSA3042.1.BED \
--ind SA3042 \
--pop_order GBR,CHB,LWK,NAMA \
--out results/plots/ref_bal_SA3042.png


python scripts/plot_karyogram.py \
--bed_a results/plots/bedfiles_phasedtogetherSA3146.0.BED \
--bed_b results/plots/bedfiles_phasedtogetherSA3146.1.BED \
--ind SA3146 \
--pop_order GBR,CHB,LWK,NAMA \
--out results/plots/ref_bal_SA3146.png

python scripts/plot_karyogram.py \
--bed_a results/plots/bedfiles_phasedtogetherSA3103.0.BED \
--bed_b results/plots/bedfiles_phasedtogetherSA3103.1.BED \
--ind SA3103 \
--pop_order GBR,CHB,LWK,NAMA \
--out results/plots/ref_bal_SA3103.png



for file in *.cram ; do NAME=`echo "$file" | cut -d'.' -f1` ; samtools fastq $file | gzip -c > ../fastq/${NAME}.fastq.gz


Meeting with Gerald:


rfmix does internal phasing, messes with order of haplotypes from the vcf. need to get the phase / haplotype order back. rfmix 1.7 has allelesrephased, but 2.0 doesnt
no genotype information in output files
change code to output genotype information

one issue in the pipeline:

merge ibd segments: for some reason the last two columns in output are supposed to have full stops but they have commas


for i in $(seq 1 22); do singularity exec /home/gerald/Documents/PhD/papers/singularity_containers/rfmix_old.sif RFMix_PopPhased -f /home/gerald/Documents/PhD/papers/paper2/results/rfmix/rfmix_input/rfmix_chr${i}_comb_ref_nama_40 -o

/home/gerald/Documents/PhD/papers/paper2/results/rfmix/rfmix_old_chr${i}_comb_ref_nama_40 --n-threads=4 --minimum-snps=15 --maximum-snps=15 -w 0.015; done

singularity shell


## Simple test: just shapeit phasing --> rfmix




shapeit -B results/phasing/cdb_nama_0.05.ref.chr1 -M /share/hennlab/reference/recombination_maps/genetic_map_HapMapII_GRCh37/chr1.gmap.txt --input-ref /share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz /share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz /share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3.sample --duohmm -W 5 -O results/phasing/cdb_nama_0.05.ref.chr1.phased --output-log shapeit_log/cdb_nama_0.05.ref_chr1.phased -T 10
