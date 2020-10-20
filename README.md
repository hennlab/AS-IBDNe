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

3. Run Snakemake

```bash
nice /share/hennlab/progs/miniconda3/bin/snakemake --config gmap=/share/hennlab/reference/recombination_maps/genetic_map_HapMapII_GRCh37/ dataset=nama_tgp_qc ref=/share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3 ids_ref=data/ref_nama.inds -j 10
```
