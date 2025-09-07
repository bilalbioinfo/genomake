#!/bin/bash -le


#SBATCH --job-name=bcftools_consensus
#SBATCH --partition=shared
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=slurm_output/slurm_%A_%a.out
#SBATCH --error=slurm_output/slurm_%A_%a.out


ml bcftools/1.20

### input parameters
output_filename=$1
bcf=$2
reference=$3
output_dir=$4


date
printf "Calling consensus using bcftools ............................\n"
bcftools consensus -f ${reference} -o ${output_dir}${output_filename} ${bcf}
printf "......................................................DONE\n\n"
date
