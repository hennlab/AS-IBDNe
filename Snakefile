"""
Title: AS-IBDne
Description: Custom pipeline to infer ancestry-specific historical effective population sizes using IBDne
Authors: Mira Mastoras, Austin Reynolds
"""

import os
import pandas as pd
import numpy as np

# -------
# SET-UP
# -------

DATASET = config['dataset']
GENMAP = config['rfmix_genmap'] # The genetic map file should contain the map for the entire genome (all chromosomes).
GMAP = config['chr_gmap'] # path to genetic map files for each chromosome
SMPMAP = config['smpmap']
ADMIX = config['admix_samples']
REF=config['ref'] # path to reference haps/sample/legend files, including prefix.
MINCM=config['mincM']
COLORS= config['colors']
CHR = [i for i in range(1,23)]

# Get string containing list of reference pops for scripts msp_to_bed.R and plot_karyogram.R
samples = pd.read_csv(SMPMAP, sep='\t', header = None)
ref_pops = np.unique(np.array(samples.iloc[:,1]))
ref_pops = sorted(ref_pops.tolist())
NANC=len(ref_pops) # get number of ref ancestries
ref_pops = ",".join(ref_pops)



# get sequence of 1 : num ancestries to define ibdne input
anc_list = [i for i in range(1,NANC+1)]

# convert admixed samples to list for plotting karyograms
samples = pd.read_csv(ADMIX, header = None)
admix_toplot = np.array(samples)
admix_toplot = admix_toplot.tolist()

# Input functions
def chrom_combine_inputs(wildcards):
  files = expand("results/phasing/{dataset}.chr{chrnum}.phased.vcf", chrnum=CHR, dataset=DATASET)
  return files

def prep_ibdne_input(wildcards):
  files = expand("results/IBD-segs/{dataset}.chr{chrnum}.phased.filled.allanc.ibd", dataset=DATASET, chrnum=CHR)
  return files

# ---------
#  Targets
# ---------

rule all:
  input:
    expand("results/IBDne/{dataset}.anc{anc}_{cm}cM.ibdne.ne", dataset = DATASET, anc=anc_list, cm=MINCM),
    expand("results/plots/{dataset}.{adm}.rfmix.karogram.png",  adm = admix_toplot, dataset = DATASET),
    expand("results/plots/{dataset}.ibdne.pdf",  dataset = DATASET)

# ---------
#  Phasing
# ---------

rule break_chrom:
  input:
    multiext("data/{dataset}", ".bed", ".bim", ".fam")
  output:
    multiext("results/phasing/{dataset}.chr{chrnum}", ".bed", ".bim", ".fam")
  params:
    "results/phasing/{dataset}.chr{chrnum}"
  shell:
    """
    plink --bfile data/{wildcards.dataset} --chr {wildcards.chrnum} --make-bed --out {params}
    """

rule trim_phase:
  input:
    multiext("results/phasing/{dataset}.chr{chrnum}", ".bed", ".bim", ".fam")
  output:
    multiext("results/phasing/{dataset}.chr{chrnum}.phased", ".haps", ".sample"),
    "shapeit_log/{dataset}_chr{chrnum}.phased.log"
  params:
    inputmap = GMAP+"chr{chrnum}.gmap.txt",
    in_pre="results/phasing/{dataset}.chr{chrnum}",
    out_pre="results/phasing/{dataset}.chr{chrnum}.phased",
    log_pre="shapeit_log/{dataset}_chr{chrnum}.phased"
  threads: 10
  shell:
    """
    bash scripts/shapeit_iterate.sh {params.in_pre} {params.inputmap} {params.log_pre} {params.out_pre} {REF}_chr{wildcards.chrnum}.hap.gz {REF}_chr{wildcards.chrnum}.legend.gz {REF}.sample {threads}
    """


# ----------
# Run RFmix
# ----------

rule make_vcf:
  input:
    multiext("results/phasing/{dataset}.chr{chrnum}.phased", ".haps", ".sample")
  output:
    vcf = "results/phasing/{dataset}.chr{chrnum}.phased.vcf"
  params:
    in_pre="results/phasing/{dataset}.chr{chrnum}.phased",
  shell:
    """
    shapeit -convert --input-haps {params.in_pre} --output-vcf {output.vcf}
    """

rule split_vcf:
  input:
    "results/phasing/{dataset}.chr{chrnum}.phased.vcf"
  output:
    ref = "results/phasing/{dataset}.chr{chrnum}.phased.ref.vcf",
    admix = "results/phasing/{dataset}.chr{chrnum}.phased.admix.vcf"
  shell:
    """
    bgzip -c {input} > {input}.gz
    tabix -p vcf {input}.gz
    bcftools view -S {ADMIX} {input}.gz > {output.admix}
    bcftools view -S ^{ADMIX} {input}.gz > {output.ref}
    """

rule run_rfmix:
  input:
    ref = "results/phasing/{dataset}.chr{chrnum}.phased.ref.vcf",
    admix = "results/phasing/{dataset}.chr{chrnum}.phased.admix.vcf"
  output:
    "results/RFmix/{dataset}.chr{chrnum}.msp.tsv",
    "results/RFmix/{dataset}.chr{chrnum}.fb.tsv"
  params:
    "results/RFmix/{dataset}.chr{chrnum}"
  shell:
    """
    rfmix -f {input.admix} -r {input.ref} -m {SMPMAP} -g {GENMAP} -o {params} --chromosome={wildcards.chrnum}
    """

# ----------------
# Get IBD segments
# ----------------

rule convert_germline:
  input:
    multiext("results/phasing/{dataset}.chr{chrnum}.phased", ".haps", ".sample")
  output:
    multiext("results/phasing/{dataset}.chr{chrnum}.phased", ".ped", ".map")
  params:
    prefix = "results/phasing/{dataset}.chr{chrnum}.phased",
    map = GMAP+"chr{chrnum}.gmap.txt"
  shell:
    """
    python scripts/shapeit_to_germline.py {params.prefix} {params.map}
    """

rule refined_ibd:
  input:
    vcf = "results/phasing/{dataset}.chr{chrnum}.phased.vcf",
    map = "results/phasing/{dataset}.chr{chrnum}.phased.map"
  output:
    multiext("results/IBD-segs/{dataset}.chr{chrnum}.phased", ".log", ".hbd.gz", ".ibd.gz")
  params:
    "results/IBD-segs/{dataset}.chr{chrnum}.phased"
  threads: 10
  shell:
    """
    java -jar progs/refined-ibd.17Jan20.102.jar gt={input.vcf} map={input.map} out={params} nthreads={threads}
    """

rule fill_gaps:
  input:
    ibd = "results/IBD-segs/{dataset}.chr{chrnum}.phased.ibd.gz",
    vcf = "results/phasing/{dataset}.chr{chrnum}.phased.vcf",
    map = "results/phasing/{dataset}.chr{chrnum}.phased.map"
  output:
    "results/IBD-segs/{dataset}.chr{chrnum}.phased.filled.ibd"
  shell:
    """
    gunzip -c {input.ibd} | java -jar progs/merge-ibd-segments.16May19.ad5.jar {input.vcf} {input.map} 0.6 1 > {output}
    """

# ------------------------
# Add Ancestry from RFMIX
# ------------------------

rule convert_msp:
  input:
    msp = "results/RFmix/{dataset}.chr{chrnum}.msp.tsv",
    bim = "results/phasing/{dataset}.chr{chrnum}.bim"
  output:
    vit = "results/RFmix/{dataset}.chr{chrnum}.vit.tsv",
    fixed = "results/RFmix/{dataset}.chr{chrnum}.vit.tsv.fixed"
  shell:
    """
    python scripts/msp_to_vit.py {input.msp} {input.bim} {output.vit}
    bash scripts/reformat_vit.sh {input.msp} {input.bim} {output.vit} {wildcards.chrnum}
    """

rule add_ancestry:
  input:
    ibd = "results/IBD-segs/{dataset}.chr{chrnum}.phased.ibd.gz",
    filled_ibd = "results/IBD-segs/{dataset}.chr{chrnum}.phased.filled.ibd",
    vit = "results/RFmix/{dataset}.chr{chrnum}.vit.tsv.fixed"
  output:
    "results/IBD-segs/{dataset}.chr{chrnum}.phased.filled.allanc.ibd"
  shell:
    """
    python scripts/filter_gapfilled_ibd_ancestry.py {input.ibd} {input.filled_ibd} {input.vit} {NANC} > {output}
    """

rule prep_ibdne:
  input:
    prep_ibdne_input
  output:
    "results/IBD-segs/{dataset}.phased.filled.anc{anc}.ibd",
    "results/IBDne/{dataset}.anc{anc}_npairs"
  shell:
    """
    bash scripts/sep_anc_comb_chr.sh {DATASET} {wildcards.anc} {ADMIX}
    """

# ---------------------------
#  Make map files for IBDne
# ---------------------------

rule combine_chrom:
  input:
    chrom_combine_inputs
  output:
    "results/phasing/{dataset}.allchr.phased.vcf"
  shell:
    """
    sed -n -e '/^#/ p' results/phasing/{wildcards.dataset}.chr1.phased.vcf > {output} #initialize allchr.vcf with headers
    for i in {input}; do sed '/^#/d' $i >> {output} ; done
    """

rule make_allchr_map:
  input:
    "results/phasing/{dataset}.allchr.phased.vcf"
  output:
    "results/phasing/{dataset}.allchr.phased.map"
  params:
    "results/phasing/{dataset}.allchr.phased"
  shell:
    """
    plink --vcf {input} --recode --out {params}
    """

rule add_cm:
  input:
    "results/phasing/{dataset}.allchr.phased.map"
  output:
    "results/phasing/{dataset}.allchr.phased.cm.map"
  params:
    gmap = GMAP+"chr@.gmap.txt",
    in_pre = "results/phasing/{dataset}.allchr.phased",
    out_pre = "results/phasing/{dataset}.allchr.phased.cm"
  shell:
    """
    plink --file {params.in_pre} --cm-map {params.gmap} --recode --out {params.out_pre}
    """

# --------
#  IBDne
# --------

rule ibdne:
  input:
    ibd = "results/IBD-segs/{dataset}.phased.filled.anc{anc}.ibd",
    npairs = "results/IBDne/{dataset}.anc{anc}_npairs",
    map = "results/phasing/{dataset}.allchr.phased.cm.map"
  output:
    "results/IBDne/{dataset}.anc{anc}_{cm}cM.ibdne.ne"
  params:
    "results/IBDne/{dataset}.anc{anc}_{cm}cM.ibdne"
  shell:
    """
    cat {input.ibd} | java -jar progs/ibdne.07May18.6a4.jar map={input.map} nthreads=12 mincm={MINCM} npairs=`cat {input.npairs}` filtersamples=false out={params}
    """
# ----------
#  Plotting
# ----------

rule plot_ibdne:
  input:
    expand("results/IBDne/{dataset}.anc{anc}_{cm}cM.ibdne.ne", dataset=DATASET, anc=anc_list, cm=MINCM)
  output:
    "results/plots/{dataset}.ibdne.pdf"
  params:
    prefix = "results/IBDne/{dataset}.anc",
    suffix = expand("_{cm}cM.ibdne.ne", cm=MINCM)
  shell:''
    """
    Rscript scripts/plot_ibdne.R {params.prefix} {params.suffix} {output} {ref_pops} {NANC} {COLORS}
    """

rule msp_to_bed:
  input:
    expand("results/RFmix/{dataset}.chr{chrnum}.msp.tsv", dataset=DATASET, chrnum=CHR)
  output:
    "results/plots/bedfiles/{dataset}.{adm}.0.BED",
    "results/plots/bedfiles/{dataset}.{adm}.1.BED"
  threads: 10
  params:
    input = "results/RFmix/",
    out = "results/plots/bedfiles/{dataset}."
  shell:
    """
    Rscript scripts/msp_to_bed.R {params.input} {params.out} {threads} {ref_pops}
    """

rule plot_rfmix:
  input:
    a = "results/plots/bedfiles/{dataset}.{adm}.0.BED",
    b = "results/plots/bedfiles/{dataset}.{adm}.1.BED"
  output:
    "results/plots/{dataset}.{adm}.rfmix.karogram.png"
  shell:
    """
    python scripts/plot_karyogram.py --bed_a {input.a} --bed_b {input.b} --ind {wildcards.adm} --pop_order {ref_pops} --colors {COLORS} --out {output}
    """
