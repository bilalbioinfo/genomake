#!/bin/bash -le

### Author: Bilal Sharif <bilal.bioinfo@gmail.com>
### Description: Genomake - A pipeline to chop a reference genome into reads, map them to another reference genome, call variants, and generate a consensus sequence.
### Date: 2025-7-15

## project account for SLURM jobs - replace with your own account ID and update the line with #SLURM -A <ACCOUNT_ID>
ACCOUNT_ID="naiss2025-22-1155"
#SBATCH -A naiss2025-22-1155
#SBATCH --job-name=5_genomake
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm_output/slurm_%A_%a.out
#SBATCH --error=slurm_output/slurm_%A_%a.out


######################################################################################
### Set variables - these directories and files should be updated according to your project structure
######################################################################################
#region

## Project directories
DIR_WORKING=/cfs/klemming/projects/snic/sllstore2017093/homotherium/bilal/5_genomake/
DIR_SCRIPTS=/cfs/klemming/projects/snic/sllstore2017093/homotherium/bilal/0_scripts/5_genomake/
DIR_RESULTS=${DIR_WORKING}results/
DIR_DATA=${DIR_WORKING}data/

## Variables
## provide path to the fasta file that will be chopped to create fastq reads
FASTA=/cfs/klemming/projects/snic/sllstore2017093/homotherium/reference/GCA_046254995.1_Pleo_Luke/GCA_046254995.1_Pleo_Luke_hap1Y_genomic+X+MT.fasta
## provide path to the reference file that will be used for mapping and variant calling, and consensus calling
REF_FILE=/cfs/klemming/projects/snic/sllstore2017093/homotherium/reference/GCF_018350175.1_F.catus/GCF_018350175.1_F.catus.fasta


# # ## for testing purposes, - ignore this
# FASTA=/cfs/klemming/projects/snic/sllstore2017093/homotherium/reference/OK513017.1/OK513017.1.fasta
# REF_FILE=/cfs/klemming/projects/snic/sllstore2017093/homotherium/reference/NC_001700.1/NC_001700.1.fasta

cd ${DIR_WORKING}
source ${DIR_SCRIPTS}0_tools/slurm_util.sh
printf "\n\n\033[1;36m========================= Starting the Genomake pipeline =========================\033[0m\n\n\n"

#endregion
######################################################################################
# Part 1: chop the reference sequence into reads
######################################################################################
#region

## set variables
readlen=150
step_size=75
output_dir_chopped_fasta=${DIR_RESULTS}1_chopped_fasta/

## output file name
output_filename_fastq=$(basename "${FASTA}" | sed -E 's/\.(fasta|fa|fna)$//I')
mkdir -p ${output_dir_chopped_fasta}
cd ${output_dir_chopped_fasta}

printf "1. Chopping fasta file into reads of length ${readlen} with step size ${step_size}............!\n"
## Submit the job to chop the fasta file into reads
job_id_chop_fasta=$(sbatch --parsable -A "$ACCOUNT_ID" \
    ${DIR_SCRIPTS}1_chopping_fasta/0_chopping_fasta.sh \
    ${DIR_SCRIPTS}1_chopping_fasta/chop_fasta.py \
    ${FASTA} \
    ${readlen} \
    ${step_size} \
    ${output_filename_fastq} \
    ${output_dir_chopped_fasta})

wait_for_slurm_job "$job_id_chop_fasta"

printf "Chopping fasta completed successfully.\n\n\n"

#endregion
######################################################################################
### Part 2: map the chopped reads to the reference genome and filter the BAM file
######################################################################################
#region

## set variables
fastq=${output_dir_chopped_fasta}${output_filename_fastq}_${readlen}_${step_size}.fq
output_dir_mapping=${DIR_RESULTS}2_mapping/

## output file name
output_filename_bam=$(basename ${fastq} .fq)
mkdir -p ${output_dir_mapping}
cd ${output_dir_mapping}

printf "2. Mapping reads with Bowtie2 and filtering..................................!\n"
## Submit the job to map the reads to the reference genome
job_id_mapping=$(sbatch --parsable -A "$ACCOUNT_ID" \
    ${DIR_SCRIPTS}2_mapping/0_mapping_bowtie2.sh \
    ${output_filename_bam} \
    ${fastq} \
    ${REF_FILE} \
    ${output_dir_mapping})

wait_for_slurm_job "$job_id_mapping"
printf "Mapping completed successfully.\n\n\n"

#endregion
######################################################################################
### Part 3: call variants using bcftools and filter the BCF file for only alt alleles
######################################################################################
#region

## set variables
filtered_bam=${output_dir_mapping}filtered_bams/${output_filename_bam}_mq*_filtered.bam
output_dir_variants=${DIR_RESULTS}3_variants/

## output file name
output_filename_bcf=$(basename ${filtered_bam} .bam)
mkdir -p ${output_dir_variants}
cd ${output_dir_variants}

printf "3. Calling variants with bcftools call and filtering.........................!\n"
## Submit the job to call variants using bcftools
job_id_variant_calling=$(sbatch --parsable -A "$ACCOUNT_ID" \
    ${DIR_SCRIPTS}3_variant_calling/0_variant_calling_bcftools.sh \
    ${output_filename_bcf} \
    ${filtered_bam} \
    ${REF_FILE} \
    ${output_dir_variants})

wait_for_slurm_job "$job_id_variant_calling"
printf "Variant calling completed successfully.\n\n\n"

#endregion
######################################################################################
### Part 4: call consensus sequence using bcftools
######################################################################################
#region

## set variables
filtered_bcf=${output_dir_variants}${output_filename_bcf}_DP2-2_g3_no_indels_homalt.bcf
output_dir_consensus=${DIR_RESULTS}4_consensus_calling/

## output file name
ref_filename=$(basename "${REF_FILE}" | sed -E 's/\.(fasta|fa|fna)$//I')
fasta_filename=$(basename "${FASTA}" | sed -E 's/\.(fasta|fa|fna)$//I')
output_filename_fasta=${ref_filename}_updated_with_${fasta_filename}.fasta
mkdir -p ${output_dir_consensus}
cd ${output_dir_consensus}

printf "4. Calling consensus sequence with bcftools..................................!\n"
## Submit the job to call consensus sequence using bcftools
job_id_consensus_calling=$(sbatch --parsable -A "$ACCOUNT_ID" \
    ${DIR_SCRIPTS}4_consensus_calling/0_consensus_calling_bcftools.sh \
    ${output_filename_fasta} \
    ${filtered_bcf} \
    ${REF_FILE} \
    ${output_dir_consensus})

wait_for_slurm_job "$job_id_consensus_calling"

## get number of variants applied to the reference genome
alt_variants=$(cat ${output_dir_consensus}/slurm_output/slurm_${job_id_consensus_calling}*.out | grep "Applied" | awk '{print $2}')

printf "Consensus calling completed successfully! ${alt_variants} variants applied to the reference genome.\n\n\n"

#endregion
######################################################################################

printf "\033[1;36m==================== Genomake pipeline completed successfully ====================\033[0m\n\n\n"


