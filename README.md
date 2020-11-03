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
nice /share/hennlab/progs/miniconda3/bin/snakemake -R parse_pops --configfile config/config.yaml -j 10
```
