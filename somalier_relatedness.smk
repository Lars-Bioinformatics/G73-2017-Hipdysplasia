import glob

# SAMPLES, = glob_wildcards("bam/{sample}.bam")
project="G73-2017-Hoftedysplasi"
configfile: "/work/sduvarcall/"+project+"/all_samples.yaml"
# flat_list = [item for sublist in t for item in sublist]
# SAMPLES = [sample for s in config for sample in config[s]["all"]]
SAMPLES = [s.split("_",1)[0] for s in config["all"]]
print(SAMPLES)

ref = "/work/sduvarcall/resources/hg38/Homo_sapiens_assembly38.fasta"

# NGSCheckMate install dir
somalier_dir = "/work/sduvarcall/G73-2017-Hoftedysplasi/somalier"

onstart:
    shell("mkdir -p somalier_output")

rule all:
    input:
        expand("somalier_output/{project}_somalier.html", project=project)


rule extract_sites:
    input:
        bam=lambda wildcards: glob.glob("bam/{sample}*.bam".format(sample=wildcards.sample)),
        sites=somalier_dir+"/sites.hg38.vcf.gz"
    output:
        vcf="somalier_output/extracted/{sample}.somalier"
    shell:
        """
        {somalier_dir}/somalier extract -d somalier_output/extracted --sites {input.sites} -f {ref} {input.bam}
        """

rule relate:
    input:
        vcf=expand("somalier_output/extracted/{sample}.somalier", sample=SAMPLES)
    output:
        "somalier_output/{project}_somalier.html"
    params:
        prefix="somalier_output/{project}_somalier"
    shell:
        """
        {somalier_dir}/somalier relate -o {params.prefix} {input.vcf}
        """