#!/bin/bash -le


#SBATCH --job-name=mapping
#SBATCH --partition=shared
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --output=slurm_output/slurm_%A_%a.out
#SBATCH --error=slurm_output/slurm_%A_%a.out

## input parameters
sample=$1
fastq=$2
reference=$3
output_dir=$4

## parameters
mapping_quality=1
bowtie2_threads=32
samtools_threads=16


############################################################################################################
### Part 1: Mapping
############################################################################################################

output_dir_raw_bams=${output_dir}raw_bams/
mkdir -p ${output_dir_raw_bams}

## loading modules
ml -q bowtie2/2.5.4 samtools/1.20

## checking if the file already exists
if [ -f ${output_dir_raw_bams}${sample}_sorted.bam ]; then
    printf "\t\t Skipping mapping, file already exists: ${output_dir_raw_bams}${sample}_sorted.bam\n"
    exit 0
fi

## Script to map reads to the reference genome using Bowtie2 and sort the output using samtools
printf "STARTING mapping.................................................................................. ${sample}\n"
bowtie2 --sensitive -p ${bowtie2_threads} -x ${reference} -U ${fastq} | \
samtools view -F 4 -q 1 -@ ${bowtie2_threads} -bh - | \
samtools sort -@ ${bowtie2_threads} -o ${output_dir_raw_bams}${sample}_sorted.bam -
samtools index ${output_dir_raw_bams}${sample}_sorted.bam
printf "${sample}......................................................................................... completed\n\n"


############################################################################################################
### Part 2: Filter for mapping quality
############################################################################################################

output_dir_filtered_bams=${output_dir}filtered_bams/
mkdir -p ${output_dir_filtered_bams}

ml -q samtools/1.20

printf "Filtering for mapping quality............................................................ ${sample}\n"
samtools view -F 4 -q ${mapping_quality} -@ ${samtools_threads} -bh ${output_dir_raw_bams}${sample}_sorted.bam -o ${output_dir_filtered_bams}${sample}_mq${mapping_quality}_filtered.bam
samtools index ${output_dir_filtered_bams}${sample}_mq${mapping_quality}_filtered.bam
printf "${sample}...........................................................................................  completed\n\n"

