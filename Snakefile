"""
Title: AS-IBDne
Description: Custom pipeline to infer ancestry-specific historical effective population sizes using IBDne
Authors: Austin Reynolds, Mira Mastoras
"""

import os

# -------
# SET-UP
# -------

DATASET = config['dataset']
GMAP = config['gmap']
REF=config['ref']
REF_IDS=config['ids_ref']


CHR = [i for i in range(1,23)]
IBDne_THREADs = 10
RFMIX_THREADs = 10

rule all:
  input:
    expand("results/RFMIX/{dataset}.chr{chrnum}.rfmix", dataset = DATASET, chrnum = CHR)

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
    log = "shapeit_log/{dataset}_chr{chrnum}.phased.log"
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

rule make_vcf:
  input:
    multiext("results/phasing/{dataset}.chr{chrnum}.phased", ".haps", ".sample")
  output:
    vcf = "results/phasing/{dataset}.chr{chrnum}.phased.vcf",
    log = "shapeit_log/{dataset}_{chrnum}.vcf.log"
  params:
    in_pre="results/phasing/{dataset}.chr{chrnum}.phased",
  shell:
    """
    shapeit -convert --input-haps {params.in_pre} --output-vcf {output.vcf} --output-log {output.log}
    """

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
    multiext("results/IBD-segs/{dataset}.chr{chrnum}", ".log", ".hbd.gz", ".ibd.gz")
  params:
    "results/IBD-segs/{dataset}.chr{chrnum}"
  threads: 10
  shell:
    """
    java -jar progs/refined-ibd.17Jan20.102.jar gt={input.vcf} map={input.map} out={params} nthreads={threads}
    """

rule fill_gaps:
  input:
    ibd = "results/IBD-segs/{dataset}.chr{chrnum}.ibd.gz",
    vcf = "results/phasing/{dataset}.chr{chrnum}.phased.vcf",
    map = "results/phasing/{dataset}.chr{chrnum}.phased.map"
  output:
    "results/IBD-segs/{dataset}.chr{chrnum}.filled.ibd"
  shell:
    """
    gunzip -c {input.ibd} | java -jar progs/merge-ibd-segments.16May19.ad5.jar {input.vcf} {input.map} 0.6 1 > {output}
    """

rule rfmix_input:
  input:
    haps = "results/phasing/{dataset}.chr{chrnum}.phased.haps",
    sample = "results/phasing/{dataset}.chr{chrnum}.phased.sample",
    ref = REF_IDS
  output:
    multiext("results/RFMIX/{dataset}_chr{chrnum}", ".alleles", ".snp_locations", ".map")
  params:
    map = GMAP+"chr{chrnum}.gmap.txt",
    out = "results/RFMIX/{dataset}"
  shell:
    """
    cut -d" " -f2 results/phasing/{wildcards.dataset}.chr{wildcards.chrnum}.phased.sample | sed '1d' > results/RFMIX/ref.keep # all ids
    cut -f2 {input.ref} > results/RFMIX/{wildcards.dataset}.ref # reference ids
    grep -v -f results/RFMIX/{wildcards.dataset}.ref results/RFMIX/ref.keep > results/RFMIX/admix.keep # find lines unique to ref.keep - ids not in the ref set
    python scripts/shapeit2rfmix.py --shapeit_hap_ref {input.haps} --shapeit_hap_admixed {input.haps} --shapeit_sample_ref {input.sample} --shapeit_sample_admixed {input.sample} --ref_keep results/RFMIX/ref.keep --admixed_keep results/RFMIX/admix.keep --chr {wildcards.chrnum} --genetic_map {params.map} --out {params.out}
    """

rule run_rfmix:
  input:
    all = "results/RFMIX/{dataset}_chr{chrnum}.alleles",
    snp = "results/RFMIX/{dataset}_chr{chrnum}.snp_locations"
  output:
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix"
  params:
    "results/RFMIX/{dataset}.classes"
  threads: 4
  shell:
    """
    python scripts/RunRFMix.py -e 2 -w 0.2 --num-threads {threads} --use-reference-panels-in-EM --forward-backward PopPhased {input.all} {params} {input.snp} -o {output}
    """
