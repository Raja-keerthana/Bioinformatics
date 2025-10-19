# NGS Variant Calling Pipeline (Course Project)

This repository contains a Bash script (`variant_calling.sh`) developed as part of my college NGS practical coursework.

## Description
The script automates the following steps:
1. Quality control of sequencing data using FASTQC
2. Alignment to reference genome using Bowtie2
3. Conversion, sorting, and indexing of alignment files using Samtools
4. Variant calling using GATK HaplotypeCaller
5. Somatic variant calling using GATK Mutect2

## Tools Required
- Bash shell
- FastQC
- Bowtie2
- Samtools
- GATK (requires Java 11+)

## Notes
- Input FASTQ and genome files are not included due to size.
- Outputs were generated and verified during college NGS lab sessions.
