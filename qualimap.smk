__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "11/11/2020"
__version__ = "1.0"

import time, os

configfile: "/work/sduvarcall/G73-2017-Hoftedysplasi/all_samples.yaml"

SAMPLES = [sample for sample in config["all"]]
# SAMPLES = SAMPLES[0]
print(SAMPLES)

# Resources
ref = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"
bed = "/work/sduvarcall/Resources/GRCh37/MedExome_target_regions/MedExome_GRCh37_capture_targets.bed"


###############################################################################
### Rule all																###
###############################################################################
rule all:
	input:
		expand("qualimap/{sample}/qualimapReport.html", sample=SAMPLES)
		# [expand("qualimap/{sample}/qualimapReport.html", sample=sample) for sample in SAMPLES]
		



###############################################################################
### Qualimap																###
###############################################################################
'''
Map reads to reference genome with bwa
'''
rule Qualimap:
	input:
		bam="bam/{sample}_recal.bam"
	output:
		report="qualimap/{sample}/qualimapReport.html"
	params:
		outdir="qualimap/{sample}"
	threads: 24
	shell:
		"""
		qualimap bamqc \
		-bam {input.bam} \
		-gff {bed} \
		-nt {threads} \
		-c \
		-sd \
		-outdir {params.outdir} \
		--java-mem-size=12G
		"""