#!/bin/bash

# path to directories
WORK_DIR="/home/dcf-04/NGS"
GENOME_DIR="$WORK_DIR/genome"
SOFTWARE_DIR="$WORK_DIR/software"
FASTQC="$SOFTWARE_DIR/FastQC/fastqc"
ANALYSIS_DIR="$WORK_DIR/Analysis"
SAMTOOLS="$SOFTWARE_DIR/samtools-1.21/samtools"
BOWTIE2="$SOFTWARE_DIR/bowtie2-2.5.3-linux-x86_64/bowtie2"
GATK="$SOFTWARE_DIR/gatk-4.6.1.0/gatk"
FASTQ_DIR="$WORK_DIR/fastqc"
IGV_DIR="$SOFTWARE_DIR/IGV_Linux_2.16.0"
FASTQ_RESULTS="$WORK_DIR/fastqc/fastqc_results"


#Input files
GENOME_REF="$GENOME_DIR/hg38.fa"
GENOME_PREFIX="$GENOME_DIR/hg38"
DISEASE_1="$FASTQ_DIR/Disease_1.fastq.gz"
DISEASE_2="$FASTQ_DIR/Disease_2.fastq.gz"
NORMAL_1="$FASTQ_DIR/Normal_1.fastq.gz"
NORMAL_2="$FASTQ_DIR/Normal_2.fastq.gz"

#Fastqc
"$FASTQC" "$FASTQ_DIR"/*.fastq.gz -o "$FASTQ_RESULTS"
echo "Fastqc complete"

#bowtie2
DISEASE_SAM="$ANALYSIS_DIR/DisR1R2.sam"
NORMAL_SAM="$ANALYSIS_DIR/R1R2.sam"

"$BOWTIE2" -x "$GENOME_PREFIX" -1 "$NORMAL_1" -2 "$NORMAL_2" -S "$NORMAL_SAM"
echo "bowtie2 alignment complete"
"$BOWTIE2" -x "$GENOME_PREFIX" -1 "$DISEASE_1" -2 "$DISEASE_2" -S "$DISEASE_SAM"
echo "bowtie2 disease alignment complete"

# sam to bam
DISEASE_BAM="$ANALYSIS_DIR/DisR1R2.bam"
NORMAL_BAM="$ANALYSIS_DIR/R1R2.bam"

"$SAMTOOLS" view -S -b "$DISEASE_SAM" -o "$DISEASE_BAM"
echo "sam to bam for disease completed"
"$SAMTOOLS" view -S -b "$NORMAL_SAM" -o "$NORMAL_BAM"
echo "sam to bam for NORMAL completed"

# read groups
DISEASE_RG_BAM="$ANALYSIS_DIR/DisR1R2_rg.bam"
NORMAL_RG_BAM="$ANALYSIS_DIR/R1R2_rg.bam"

"$SAMTOOLS" addreplacerg -r "@RG\tID:dis\tSM:disease" "$DISEASE_BAM" -o "$DISEASE_RG_BAM"
"$SAMTOOLS" addreplacerg -r "@RG\tID:ctrl1\tSM:normal" "$NORMAL_BAM" -o "$NORMAL_RG_BAM"

#sorting
DIS_SORTED_BAM="$ANALYSIS_DIR/DisR1R2_rg_sorted.bam"
NORMAL_SORTED_BAM="$ANALYSIS_DIR/R1R2_rg_sorted.bam"

"$SAMTOOLS" sort "$DISEASE_RG_BAM" -o "$DIS_SORTED_BAM"
"$SAMTOOLS" sort "$NORMAL_RG_BAM" -o "$NORMAL_SORTED_BAM"
echo "sorting complete"

#indexing
"$SAMTOOLS" index "$DIS_SORTED_BAM"
"$SAMTOOLS" index "$NORMAL_SORTED_BAM"
echo "indexing complete"

#mutattion callinG
DIS_VCF="$ANALYSIS_DIR/DisR1R2_rg_sorted.vcf.gz"
NORMAL_VCF="$ANALYSIS_DIR/R1R2_rg_sorted.vcf.gz"

"$GATK" HaplotypeCaller -R "$GENOME_REF" -I "$DIS_SORTED_BAM" -O "$DIS_VCF"
"$GATK" HaplotypeCaller -R "$GENOME_REF" -I "$NORMAL_SORTED_BAM" -O "$NORMAL_VCF"
echo "mutation calling complete"

#somatic callinG

SOMATIC_VCF="$ANALYSIS_DIR/Somatic_calling.vcf.gz"

"$GATK" Mutect2 -R "$GENOME_REF" -I "$DIS_SORTED_BAM" -I "$NORMAL_SORTED_BAM" -O "$SOMATIC_VCF"
echo "somatic calling complete"
