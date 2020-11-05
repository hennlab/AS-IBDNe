"""
Title: AS-IBDne
Description: Custom pipeline to infer ancestry-specific historical effective population sizes using IBDne
Authors: Austin Reynolds, Mira Mastoras
"""

import os

# -------
# SET-UP
# -------

DATASET = config['dataset'] # Name of the dataset. Must match input files data/{dataset}.bed .bim .fam
GMAP = config['gmap'] # path to genetic map files for each chromosome
REF=config['ref'] # path to reference haps/sample/legend files, including prefix.
REF_IDS=config['ids_ref'] # file containing reference sample IDS (col2) and population (col1)
ADMIX_IDS=config['ids_admix'] # file containing admixed sample IDs (col2) and population (col1)
POP_FILE= config['pops'] # file listing all reference populations (exclude the admixed population)

POPS= [line.strip() for line in open(POP_FILE, 'r')]
CHR = [i for i in range(1,23)]
IBDne_THREADs = 10
RFMIX_THREADs = 10

def generate_pop_input(wildcards):
    files = expand("results/RFMIX/{dataset}.{pops}.keep", pops=POPS, dataset=DATASET)
    return files

LIST = ["results/RFMIX/" + DATASET + "."+ p + ".keep" for p in POPS]
s = ","
POP_INPUT_LIST = s.join(LIST)


rule all:
  input:
    expand("results/RFMIX/{dataset}.chr{chrnum}.rfmix.0.ForwardBackward.txt", dataset = DATASET, chrnum = CHR),
    expand("results/RFMIX/{dataset}.chr{chrnum}.rfmix.1.ForwardBackward.txt", dataset = DATASET, chrnum = CHR),
    expand("results/RFMIX/{dataset}.chr{chrnum}.rfmix.2.ForwardBackward.txt",dataset = DATASET, chrnum = CHR),
    expand("results/RFMIX/{dataset}.chr{chrnum}.rfmix.0.Viterbi.txt",dataset = DATASET, chrnum = CHR),
    expand("results/RFMIX/{dataset}.chr{chrnum}.rfmix.1.Viterbi.txt", dataset = DATASET, chrnum = CHR),
    expand("results/RFMIX/{dataset}.chr{chrnum}.rfmix.2.Viterbi.txt", dataset = DATASET, chrnum = CHR),
    expand("results/RFMIX/{dataset}.chr{chrnum}.rfmix.allelesRephased0.txt",dataset = DATASET, chrnum = CHR),
    expand("results/RFMIX/{dataset}.chr{chrnum}.rfmix.allelesRephased1.txt",dataset = DATASET, chrnum = CHR),
    expand("results/RFMIX/{dataset}.chr{chrnum}.rfmix.allelesRephased2.txt",dataset = DATASET, chrnum = CHR),
    expand("results/RFMIX/{dataset}.{pops}.keep", pops=POPS, dataset=DATASET)


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
    sample = "results/phasing/{dataset}.chr{chrnum}.phased.sample"
  output:
    multiext("results/RFMIX/{dataset}_chr{chrnum}", ".alleles", ".snp_locations", ".map")
  params:
    map = GMAP+"chr{chrnum}.gmap.txt",
    out = "results/RFMIX/{dataset}"
  shell:
    """
    cut -d" " -f2 {REF_IDS} > data/reference.samples.keep
    cut -d" " -f2 {ADMIX_IDS} > data/admix.samples.keep
    python scripts/shapeit2rfmix.py --shapeit_hap_ref {input.haps} --shapeit_hap_admixed {input.haps} --shapeit_sample_ref {input.sample} --shapeit_sample_admixed {input.sample} --ref_keep data/reference.samples.keep  --admixed_keep data/admix.samples.keep --chr {wildcards.chrnum} --genetic_map {params.map} --out {params.out}
    """

rule parse_pops:
  input:
    REF_IDS
  output:
    "results/RFMIX/{dataset}.{pops}.keep"
  shell:
    """
    grep "{wildcards.pops}" {input} > {output}
    """

rule fix_classes:
  input:
    pop = generate_pop_input,
    samp = "results/RFMIX/{dataset}.sample"
  output:
    "results/RFMIX/{dataset}.fix.classes"
  shell:
    """
    python scripts/classes.py --ref {POP_INPUT_LIST} --sample {input.samp} --out {output}
    """

rule run_rfmix:
  input:
    all = "results/RFMIX/{dataset}_chr{chrnum}.alleles",
    snp = "results/RFMIX/{dataset}_chr{chrnum}.snp_locations",
    cl = "results/RFMIX/{dataset}.fix.classes"
  output:
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix.0.ForwardBackward.txt",
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix.1.ForwardBackward.txt",
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix.2.ForwardBackward.txt",
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix.0.Viterbi.txt",
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix.1.Viterbi.txt",
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix.2.Viterbi.txt",
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix.allelesRephased0.txt",
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix.allelesRephased1.txt",
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix.allelesRephased2.txt",
  params:
    "results/RFMIX/{dataset}.chr{chrnum}.rfmix"
  threads: 4
  shell:
    """
    python scripts/RunRFMix.py -e 2 -w 0.2 --num-threads {threads} --use-reference-panels-in-EM --forward-backward PopPhased {input.all} {input.cl} {input.snp} -o {params}
    """
