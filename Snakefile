import pandas as pd
import os

#############################################################################
# J. Baran-Gale, M. Morgan
# WTSA - SmartSeq2 pipeline
# This module take FASTQ files & trims, aligns and quant them using the 
# STAR aligner. Note: there are several choices for quant (1) STAR quant, 
# 	(2) featureCounts, (3) salmon
#  snakemake --cluster "other"
#  snakemake --cluster "qsub"
#  snakemake -c 'qsub -V  -q igmm_long -pe sharedmem 8 -l h_vmem=8G -j y -cwd' --jobs=100
#############################################################################
configfile: "siteProfiles/configEBI.yaml"
configfile: "config.yaml"


starIndexPrefix = os.path.join(config['ref_mm10']['referenceFolder'],config['ref_mm10']['starIndex'])
annotation = os.path.join(config['ref_mm10']['referenceFolder'],config['ref_mm10']['annotation'])
transcriptome = os.path.join(config['ref_mm10']['referenceFolder'],config['ref_mm10']['transcriptome'])


BASE_DIR = "/nfs/leia/research/marioni/mikemorgan/"
WKDIR = BASE_DIR + "Thymus/"

DIRS = ['Trimmed/','Aligned/','fastqs/QC','FC_QUANT','salmon_QUANT', 'temp.dir', 'Dedup.dir']

workdir: WKDIR


SAMPLES,dummy = glob_wildcards("fastqs/rawdata/{sample}/{sample1}_R1_001.fastq.gz")

#############################################################################
# Define a set of tasks to run locally
############################################################################
localrules: all, dirs

#############################################################################
# Input Rule
#############################################################################
rule all:
	input: 
		DIRS,
		expand("fastqs/rawdata/{sample}/{sample}_R1_001.fastq.gz",sample=SAMPLES),
		expand("fastqs/rawdata/{sample}/{sample}_R2_001.fastq.gz",sample=SAMPLES),
		expand("fastqs/QC/{sample}_R1_001_fastqc.html",sample=SAMPLES),
		expand("fastqs/QC/{sample}_R2_001_fastqc.html",sample=SAMPLES)
		#expand("Aligned/{sample}.Aligned.sortedByCoord.out.bam",sample=SAMPLES)
		#expand("salmon_QUANT/{sample}/quant.sf",sample=SAMPLES)

#############################################################################
# DIR Rule
#############################################################################
rule dirs:
	output: DIRS
	shell: "mkdir -p "+' '.join(DIRS)

#############################################################################
# fastqc
#############################################################################
rule run_fastqc:
	"""
	Run FastQC to generate a html report with a selection of QC modules
	"""
	input:
		R1="fastqs/rawdata/{sample}/{sample}_R1_001.fastq.gz",
		R2="fastqs/rawdata/{sample}/{sample}_R2_001.fastq.gz"
	output:
		"fastqs/QC/{sample}_R1_001_fastqc.html",
		"fastqs/QC/{sample}_R2_001_fastqc.html"
	params:
		qcEX=config['FQC']

	threads: 12
	shell:
		"{params.qcEX}  -o fastqs/QC/  --noextract  --threads {threads}  --dir temp.dir/  {input}"
#############################################################################
# Trim adaptor
#############################################################################
rule trim_adaptor:
	"""
	Use trimmomatic to remove adaptor contamination and quality trim reads to a minimum length
	NB: this should be user-configured to handle different read lengths and trimming parameter values
	"""
	input:
		R1="fastqs/rawdata/{sample}/{sample}_R1_001.fastq.gz",
		R2="fastqs/rawdata/{sample}/{sample}_R2_001.fastq.gz"
	output:
		R1 = "Trimmed/{sample}_R1_001.fastq.gz",
		R2 = "Trimmed/{sample}_R2_001.fastq.gz"
	threads: 12
	params: 
		trEX=config['TRIM'],
		trPA=config['Tparam']
	log:
		"logs/trim_adaptor/{sample}.log"
	shell:
	# trimmomatic needs to be in the PATH variable ideally, and the adaptor sequences in a common location
	# I think to make this as universal as possible I'll need to create a script that executes trimmomatic on
	# the command line with the relevant options - still needs to know where the trimmomatic binary is though
	# is this a job for conda?  Could we deploy this pipeline as a conda environment, that way we can
	# control the paths to software binaries...
		"""java -jar {params.trEX} PE -threads {threads} -phred33 {input.R1} {input.R2} {output.R1} {output.R1}.unpaired {output.R2} {output.R2}.unpaired {params.trPA}"""
#############################################################################
# Align Rule
#############################################################################
rule star_align:
	input:
		R1="Trimmed/{sample}_R1_001.fastq.gz",
		R2="Trimmed/{sample}_R2_001.fastq.gz"
	output:
		"Aligned/{sample}.Aligned.sortedByCoord.out.bam"
	params:
		starEX=config['STAR'],
		prefix = "Aligned/{sample}.",
		readFilesCommand = config['params']['star']['readFilesCommand'], 
		outSAMtype = config['params']['star']['outSAMtype'],
		outSAMattributes = config['params']['star']['outSAMattributes'],
		outSAMunmapped = config['params']['star']['outSAMunmapped'],
		quantMode = config['params']['star']['quantMode']
	threads: 12
	shell: 
		"{params.starEX} --runThreadN {threads}  --genomeDir {starIndexPrefix} --readFilesIn {input.R1} {input.R2} --readFilesCommand {params.readFilesCommand} --outFileNamePrefix {params.prefix} --outSAMtype {params.outSAMtype} --outSAMattributes {params.outSAMattributes} --outSAMunmapped {params.outSAMunmapped} --quantMode {params.quantMode} "
		#	os.system(command)

#############################################################################
# Deduplicate positional duplicates with PicardTools
#############################################################################
rule dedup_bams:
    input: bam="Aligned{sample}.Aligned.sortedByCoord.out.bam"
    output: bam="Dedup.dir/{sample}.dedup.bam",
            metrics="Dedup.dir/{sample}.metrics.txt"
    threads: 12
    params: PicardEX=config['PICARD']
    shell : """
            java -jar {params.PicardEX} MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH
            """

#############################################################################
# salmon
#############################################################################
rule salmon_counts:
	input: bam="Dedup.dir/{sample}.dedup.bam" 
	output: "salmon_QUANT/{sample}/quant.sf"
	threads: 12
	params: 
		salEX=config['SALMON'],
		salStrand="IU",
		anno=config['ref_mm10']['transcriptome']
	shell: """
		{params.salEX} quant -t {anno} -l {params.salStrand} -p {threads} -a {input.bam} -o salmon_QUANT/{sample}
	"""
#############################################################################
# featureCounts
#############################################################################
rule feature_counts:
	input: anno="{annotation}", bam="Dedup.dir/{sample}.dedup.bam" 
	output: "FC_QUANT/gene_counts.txt"
	threads: 12
	params: 
		fcEX=config['FC']
	shell: """
		{params.fcEX} -p -s 0 -T {threads} -t exon -g gene_id -a {input.anno} -o {output[0]} {input.bam} &> {output[0]}.log
	"""

