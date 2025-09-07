#!/bin/bash -le


#SBATCH --job-name=chopping_fasta
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_output/slurm_%A_%a.out
#SBATCH --error=slurm_output/slurm_%A_%a.out


## input parameters
script=$1
fasta=$2
readlen=$3
step_size=$4
output_filename=$5
output_dir=$6


ml -q bioinfo-tools biopython/1.80-py3.10.8

printf "\nSTARTING chopping fasta...............................${output_filename}\n"
python3 ${script} ${fasta} ${readlen} ${step_size} ${output_dir}${output_filename}_${readlen}_${step_size}.fq
printf "DONE chopping fasta...............................${output_filename}\n\n"

