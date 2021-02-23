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
```

2. QC beforehand
```bash
plink --bfile cdb_nama.ref --geno 0.05 --make-bed --out cdb_nama_0.05.ref
plink --bfile cdb_nama.admix --geno 0.05 --make-bed --out cdb_nama_0.05.admix
```

3. Required files and folder structure
- data
   - {dataset}.ref.bim /bed /fam
   - {dataset}.admix.bim /bed /fam
- config.yaml
- genetic map file : "The genetic map file is tab delimited text containing at least 3 columns. The first 3 columns are intepreted as chromosome, physical position in bp, genetic position in cM. Any number of columns or other information may follow, it is ignored. The chromosome column is a string token (which may be an string of digits) that must match those used in the VCF/BCF inputs. The genetic map file should contain the map for the entire genome (all chromosomes). Blank lines and lines beginning with a '#' are ignored."
- sample map file: tab delimited text with at least two columns. The first column gives the sample name or identifier, which must match the one used in the reference VCF/BCF. The second column is a string naming a subpopulation and may contain spaces (e.g., "European", or "East African")

config file:
- dataset: {dataset}.ref.admix
- genmap: genetic map file for all chromosomes
- smpmap: sample map file
- gmap: location of genetic map files for individual chromosomes.
- ref: path to reference haps/sample/legend files, including prefix


4. Set up config file

```yaml
dataset: cdb_nama_0.05

genmap: /share/hennlab/reference/recombination_maps/genetic_map_AfricanAmerian/AAmap.superChr.map.rfmix

smpmap: cdb_nama_0.05.samplefile.txt

ref: /share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3/1000GP_Phase3

gmap: /share/hennlab/reference/recombination_maps/genetic_map_HapMapII_GRCh37/
```

5. Run pipeline

```
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate IBDne-env
nice /share/hennlab/progs/miniconda3/bin/snakemake --configfile cdb_nama_0.05.config.yaml -j 20 --keep-going
```

# need to go back and sort chr 13 in the genetic map file

6. Plotting rfmix results - scripts from Gerald Van der Eeden

1. Convert msp to bed format

```
Rscript scripts/msp_to_bed.R results/RFmix results/plots/bedfiles 10 GBR CHB LWK NAMA CDB
```

2. Use alicia martins karyogram script

python scripts/plot_karyogram.py \
--bed_a results/plots/bedfiles/cdb_nama_0.05SA3042.1.BED \
--bed_b results/plots/bedfiles/cdb_nama_0.05SA3042.0.BED \
--ind SA3042 \
--pop_order GBR,CHB,LWK,NAMA \
--out results/plots/SA3042.png

cdb_nama_0.05SA3103.0.BED
python scripts/plot_karyogram.py \
--bed_a results/plots/bedfiles/cdb_nama_0.05SA3103.0.BED \
--bed_b results/plots/bedfiles/cdb_nama_0.05SA3103.1.BED \
--ind SA3103 \
--pop_order GBR,CHB,LWK,NAMA \
--out results/plots/SA3103.png


python scripts/plot_karyogram.py \
--bed_a results/plots/bedfiles/cdb_nama_0.05SA3159.0.BED \
--bed_b results/plots/bedfiles/cdb_nama_0.05SA3159.1.BED \
--ind SA3159 \
--pop_order GBR,CHB,LWK,NAMA \
--out results/plots/SA3159.png



python scripts/plot_karyogram.py \
--bed_a results/plots/bedfiles/cdb_nama_0.05SA3146.0.BED \
--bed_b results/plots/bedfiles/cdb_nama_0.05SA3146.1.BED \
--ind SA3146 \
--pop_order GBR,CHB,LWK,NAMA \
--out results/plots/SA3146.png
