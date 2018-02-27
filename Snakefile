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
#  snakemake -c 'qsub -V  -q igmm_long -pe sharedmem 8 -l h_vmem=8G -j y -cwd' --jobs=100 >runLog_DATE 2>&1 &
#############################################################################
configfile: "siteProfiles/configIGMM.yaml"
configfile: "config.yaml"


starIndexPrefix = os.path.join(config['ref_mm10']['referenceFolder'],config['ref_mm10']['starIndex'])
annotation = os.path.join(config['ref_mm10']['referenceFolder'],config['ref_mm10']['annotation'])
transcriptome = os.path.join(config['ref_mm10']['referenceFolder'],config['ref_mm10']['transcriptome'])


BASE_DIR = "/exports/eddie/scratch/jbarang/"
WKDIR = BASE_DIR + "TEST/" 

DIRS = ['Trimmed/','Aligned/','fastqs/QC','FC_QUANT','salmon_QUANT']

workdir: WKDIR


SAMPLES,dummy = glob_wildcards("fastqs/rawdata/{sample}/{sample1}_R1_001.fastq.gz")

#############################################################################
# Input Rule
#############################################################################
rule all:
	input: 
		DIRS,
		expand("fastqs/rawdata/{sample}/{sample}_R1_001.fastq.gz",sample=SAMPLES),
		expand("fastqs/rawdata/{sample}/{sample}_R2_001.fastq.gz",sample=SAMPLES),
		expand("fastqs/QC/{sample}_R1_001_fastqc.html",sample=SAMPLES),
		expand("fastqs/QC/{sample}_R2_001_fastqc.html",sample=SAMPLES),
		expand("Aligned/{sample}.Aligned.toTranscriptome.out.bam",sample=SAMPLES),
		expand("salmon_QUANT/{sample}_sf",sample=SAMPLES)

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
		"fastqs/rawdata/{sample}/{sample}_R1_001.fastq.gz",
		"fastqs/rawdata/{sample}/{sample}_R2_001.fastq.gz"
	output:
		"fastqs/QC/{sample}_R1_001_fastqc.html",
		"fastqs/QC/{sample}_R2_001_fastqc.html"
	params:
		qcEX=config['FQC'],
		dir='/tmp'

	threads: 8
	shell:
		"{params.qcEX}  -o fastqs/QC/  --noextract  --threads {threads}  --dir {params.dir}  {input}"
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
	threads: 8
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
		"Aligned/{sample}.Aligned.toTranscriptome.out.bam"
	params:
		starEX=config['STAR'],
		prefix = "Aligned/{sample}.",
		readFilesCommand = config['params']['star']['readFilesCommand'], 
		outSAMtype = config['params']['star']['outSAMtype'],
		outSAMattributes = config['params']['star']['outSAMattributes'],
		outSAMunmapped = config['params']['star']['outSAMunmapped'],
		quantMode = config['params']['star']['quantMode']
	threads: 8
	shell: 
		"{params.starEX} --runThreadN {threads}  --genomeDir {starIndexPrefix} --readFilesIn {input.R1} {input.R2} --readFilesCommand {params.readFilesCommand} --outFileNamePrefix {params.prefix} --outSAMtype {params.outSAMtype} --outSAMattributes {params.outSAMattributes} --outSAMunmapped {params.outSAMunmapped} --quantMode {params.quantMode} "

#############################################################################
# salmon
#############################################################################
rule salmon_counts:
	input: 
		bam="Aligned/{sample}.Aligned.toTranscriptome.out.bam" 
	output: "salmon_QUANT/{sample}_sf"
	threads: 8
	params: 
		salEX=config['SALMON'],
		salStrand="IU"
	shell: """
		{params.salEX} quant -t {transcriptome} -l {params.salStrand} -p {threads} -a {input.bam} -o {output}
	"""
#############################################################################
# featureCounts
#############################################################################
rule feature_counts:
	input: anno="{annotation}", bam="Aligned/{sample}.Aligned.toTranscriptome.out.bam" 
	output: "FC_QUANT/gene_counts.txt"
	threads: 8
	params: 
		fcEX=config['FC']
	shell: """
		{params.fcEX} -p -s 0 -T {threads} -t exon -g gene_id -a {input.anno} -o {output[0]} {input.bam} &> {output[0]}.log
	"""

