import os
import pandas as pd

## Load config values
#configfile: "config/defaults.yaml"
configfile: "config/events.yaml"
#configfile: "config/parameters.yaml"
#configfile: "config/annotationPaths.yaml"

# path to the reference genome fasta
#genome_path = config["genome"][config["genome_name"]]
#genome_name = {"hg19": "GRCh37.75", "hg38": "GRCh38.100", "hg38_no_alt_TCGA_HTMCP_HPVs": "GRCh38.100"}[config["genome_name"]]

# Ensembl 100 genes
#gene_anno = config["annotation"][config["genome_name"]]

# samples
events_dict = config["events"]
event_ids = events_dict.keys()

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
		expand("output/scratch/{event}/region_cna.txt", event=event_ids)

### -------------------------------------------------------------------
### get depth of regions before and after integration sites
### -------------------------------------------------------------------

rule subset_vcf:
    input:
        vcf = lambda w: config["events"][w.event]["vcf"],
        reads = lambda w: config["events"][w.event]["reads"],
    output:
        vcf = "output/scratch/{event}/subset.vcf",
        bed = "output/scratch/{event}/regions.bed"
    conda: "config/conda.yaml"
    shell:
        "python scripts/lookForSVs.py {input.vcf} {input.reads} {output.vcf} {output.bed}"

rule filtRegions:
    input:
        bed = "output/scratch/{event}/regions.bed",
    output:
        "output/scratch/{event}/regionsFilt.bed"
    conda: "config/conda.yaml"
    shell:
        "awk '{{ if ($2 != $3) print $0 }}' {input.bed} > {output}"

rule positionDepth:
    input:
        bed = "output/scratch/{event}/regionsFilt.bed",
        bam = lambda w: config["events"][w.event]["bam"]
    output:
        "output/scratch/{event}/position_depth.bed"
    conda: "config/conda.yaml"
    shell:
        "samtools depth {input.bam} -b {input.bed} | awk '{{print $1\"\\t\"$2\"\\t\"$2+1\"\\t\"$3}}' > {output}"

rule select_chr:
    input:
        bed = "output/scratch/{event}/regionsFilt.bed",
        depth = "output/scratch/{event}/position_depth.bed"
    output:
        bed = "output/scratch/{event}/regions_chr.bed",
        depth = "output/scratch/{event}/position_depth_chr.bed",
    conda: "config/conda.yaml"
    shell:
        """
        cat {input.bed} | grep 'chr' > {output.bed}
        cat {input.depth} | grep 'chr' > {output.depth}
        """

rule regionMeanDepth:
    input:
        bed = "output/scratch/{event}/regionsFilt.bed",
        depth = "output/scratch/{event}/position_depth.bed",
    output:
        "output/scratch/{event}/regionsMeanDepth.bed"
    conda: "config/conda.yaml"
    shell:
        "bedtools map -a {input.bed}  -b {input.depth} -c 4 -o mean > {output}"

rule bedpe:
    input:
        vcf = "output/scratch/{event}/subset.vcf",
    output: 
        "output/scratch/{event}/subset.bedpe"
    conda: "config/conda.yaml"
    shell:
        "/gsc/software/linux-x86_64-centos7/survivor-1.0.3/bin/SURVIVOR vcftobed {input.vcf} 5 -1 {output}"

rule regionCN:
    input:
        ploidy = lambda w: config["events"][w.event]["ploidy"],
        cna = lambda w: config["events"][w.event]["cna"],
        depth = "output/scratch/{event}/regionsMeanDepth.bed",
        vcf = "output/scratch/{event}/subset.bedpe"
    output: 
        "output/scratch/{event}/region_cna.txt"
    conda: "config/conda.yaml"
    shell:
        "scripts/calculateCN.R -p {input.ploidy} -c {input.cna} -d {input.depth} -o {output}"


