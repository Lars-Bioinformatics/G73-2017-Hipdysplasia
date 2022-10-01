__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "03/07/2022"
__version__ = "1.0"

import time, os

configfile: "/work/sduvarcall/G73-2017-Hoftedysplasi/fam.yaml"

# Resources
ref = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"
bed = "/work/sduvarcall/Resources/target_regions/MedExome_target_regions/MedExome_GRCh37_capture_targets.bed"
interval = "/work/sduvarcall/Resources/target_regions/MedExome_target_regions/MedExome_GRCh37_capture_targets.interval_list"


###############################################################################
### Rule all                                                                ###
###############################################################################
rule all:
	input:
		# [expand("vcf_gnomAD/{fam_name}_variants_gnomAD.vcf.gz", fam_name=fam) for fam in config]
		[expand("vcf_gnomAD/{fam_name}_variants_gnomAD_filtered.vcf.gz", fam_name=fam) for fam in config]



###############################################################################
### Bcftools annotate and filter gnomAD                                     ###
###############################################################################
'''
Annotate VCF with gnomAD AF and AC
'''
rule annotate:
	input:
		gnomad="/work/sduvarcall/resources/GRCh37/gnomad/af-only-gnomad.raw.sites.b37.vcf.gz",
		vcf="vcf_files/{sample}_variants.vcf.gz"
	output:
		vcf="vcf_gnomAD/{sample}_variants_gnomAD.vcf.gz"
	threads: 4
	shell:
		"""
		bcftools annotate -a {input.gnomad} -c gnomAD_AF:=AF,gnomAD_AC:=AC --threads {threads} -Oz -o {output.vcf} {input.vcf}
		"""
		# conda run -n bcftools \

'''
Filter VCF with gnomAD AF > 1%
'''
rule filter:
	input:
		gnomad="/work/sduvarcall/resources/GRCh37/gnomad/af-only-gnomad.raw.sites.b37.vcf.gz",
		vcf="{sample}_variants_gnomAD.vcf.gz"
	output:
		vcf="{sample}_variants_gnomAD_filtered.vcf.gz"
	threads: 4
	shell:
		"""
		bcftools filter -i "INFO/gnomAD_AF[0]>0.01" --threads {threads} -Oz -o {output.vcf} {input.vcf}
		"""
		# conda run -n bcftools \