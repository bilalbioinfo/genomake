#!/bin/bash -le


#SBATCH --job-name=bcftools_call
#SBATCH --partition=shared
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --output=slurm_output/slurm_%A_%a.out
#SBATCH --error=slurm_output/slurm_%A_%a.out


ml bcftools/1.20

### input parameters
sample=$1
bamfile=$2
reference=$3
output_dir=$4

## parameters
mapQ=1
baseQ=30

### parameters
min_depth=2
max_depth=2
SNP_gap_indels=3
bcftools_threads=32



date
printf "Calling Variants using bcftools ............................\n"
bcftools mpileup -q ${mapQ} -Q ${baseQ} -B -f ${reference} ${bamfile} --ignore-RG --threads ${bcftools_threads} -Ou | \
bcftools call -mv --threads ${bcftools_threads} -Ob -o ${output_dir}${sample}.bcf
bcftools sort -O b -o ${output_dir}${sample}_sorted.bcf ${output_dir}${sample}.bcf
printf "......................................................completed\n\n"

printf "Filtering variants using bcftools............................\n"
## filter variants based on depth
bcftools filter -i "DP>=${min_depth} & DP<=${max_depth}" -O b --threads ${bcftools_threads} \
    -o ${output_dir}${sample}_DP${min_depth}-${max_depth}.bcf ${output_dir}${sample}_sorted.bcf

## filter variants based on SNP gap from indels
bcftools filter -g ${SNP_gap_indels} -O b --threads ${bcftools_threads} -o ${output_dir}${sample}_DP${min_depth}-${max_depth}_g${SNP_gap_indels}.bcf \
    ${output_dir}${sample}_DP${min_depth}-${max_depth}.bcf

## filter out indels
bcftools filter -i 'INDEL=0' -O b --threads ${bcftools_threads} -o ${output_dir}${sample}_DP${min_depth}-${max_depth}_g${SNP_gap_indels}_no_indels.bcf \
    ${output_dir}${sample}_DP${min_depth}-${max_depth}_g${SNP_gap_indels}.bcf

## Only keep homozygous alternate variants
bcftools filter -i 'GT="1/1"' -O b --threads ${bcftools_threads} -o ${output_dir}${sample}_DP${min_depth}-${max_depth}_g${SNP_gap_indels}_no_indels_homalt.bcf \
    ${output_dir}${sample}_DP${min_depth}-${max_depth}_g${SNP_gap_indels}_no_indels.bcf

## Index the final filtered BCF file
bcftools index --threads ${bcftools_threads} ${output_dir}${sample}_DP${min_depth}-${max_depth}_g${SNP_gap_indels}_no_indels_homalt.bcf

## Get stats for the filtered variants
bcftools stats --threads ${bcftools_threads} ${output_dir}${sample}_DP${min_depth}-${max_depth}_g${SNP_gap_indels}_no_indels_homalt.bcf \
    > ${output_dir}${sample}_DP${min_depth}-${max_depth}_g${SNP_gap_indels}_no_indels_homalt.bcf.stats

printf "......................................................DONE\n\n"
date
