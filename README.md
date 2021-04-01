# AS-IBDNe
Snakemake pipeline for running ancestry-specific IBDNe


## 1. Set Up:

Creating conda environment (first time only):
```
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda create -n IBDne-env
conda activate IBDne-env
conda install -c bioconda rfmix
conda install pandas
```
Activating conda environment (before running pipeline)
```bash
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate IBDne-env
```

## 2. Required input files and structure of pipeline directory

Please create a working directory to run the Snakefile in. The structure of the directory should be as follows. Items marked by asterisk are required by the snakefile to be located where they are specified. The location of non-asterisked items are passed to the Snakefile by the user, so their location in the working directory is flexible.
```bash
|__ Snakefile*
|__ config.yaml
|__ data*
    |__ {dataset}.bim*
    |__ {dataset}.bed*
    |__ {dataset}.fam*
    |__ ref_smpmap.txt
    |__ admixed_ids.txt
|__ scripts*
    |__ shapeit_to_germline.py*
    |__ msp_to_vit.py*
    |__ filter_gapfilled_ibd_ancestry.py*
    |__ reformat_vit.sh*
    |__ sep_anc_comb_chr.sh*
    |__ adjust_npairs.py*
    |__ plot_karyogram.R*
    |__ plot_ibdne.R*
|__ progs*
    |__ ibdne.05May18.1c3.jar*
    |__ filtercolumns.jar*
    |__ refined-ibd.17Jan20.102.jar*
```

The input dataset to this pipeline is a single bim/bed/fam fileset containing both the admixed individuals and the reference individuals.

## 3. QC beforehand
```bash
plink --bfile nama_tgp --geno 0.05 --mind 0.1 --make-bed --out nama_tgp_qc
# add population label to bim file too
```

## 4. Config file
The config file will contain all the file paths that change with every run of the snakemake pipeline.

Example file:
```yaml
dataset: americans_subset_remove

rfmix_genmap: /share/hennlab/reference/recombination_maps/rfmix_combined_b37.map

smpmap: data/ref_smpmap.txt

admix_samples: data/admixed_ids.txt

ref: /share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3

chr_gmap: /share/hennlab/reference/recombination_maps/genetic_map_HapMapII_GRCh37/

nanc: 3
```

Explanation of config input parameters:
- **dataset**: Name of the dataset. needs to match prefix of input bim/bed/fam files : {dataset}.bed /.bim /.fam
- **rfmix_genmap**: /share/hennlab/reference/recombination_maps/rfmix_combined_b37.map
- **smpmap**: contains all reference individuals. space delimited file with IDs in first column and ancestry in second.
- **admix_samples**: list of admixed sample IDs. 1 column
- **ref**: prefix not including chromosome number and path to 1000 genome reference .hap.gz, .legend.gz and .sample files. for example, provide this:
```
ref: /share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3
```
if the files are named as:
```
/share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz
```
- **chr_gmap**: path to genetic map files for individual chromosomes.
- **nanc**: number of ancestries in reference panel



## 5. Run Snakemake

```bash
# Dry run: always run first with -n flag to make sure the workflow will execute properly
nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile config.yaml -j 20 -n
# Generate DAG
/share/hennlab/progs/miniconda3/bin/snakemake --configfile config.yaml -j 20 -n --rulegraph | dot -Tpng > rulegraph.png
# Run pipeline
nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile config.yaml -j 20

```

## 6. Plot karyograms after running RFMix

The snakefile will automatically generate karograms for the first three individuals in your admixed data text file. If you want to generate these plots for other individuals manually, re-run the script in rule "plot_rfmix".

## 7. Acknowledgements and sources:

- Rfmix version 2.3: https://github.com/slowkoni/rfmix/blob/master/MANUAL.md
- Scripts for processing rfmix 1.5 output to IBDne input: https://faculty.washington.edu/sguy/asibdne/
- IBDne: http://faculty.washington.edu/browning/ibdne.html
- filtercolumns.jar https://faculty.washington.edu/browning/beagle_utilities/utilities.html
