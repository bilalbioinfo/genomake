# Genomake

**Genomake** is a pipeline to **chop a reference genome into reads, map them to another reference genome, call variants, and generate a consensus sequence**.

It is designed to update a reference genome with variants from another species. This makes the reference closer to the target species, allowing reads from a more distant species (with no available reference genome) to be mapped successfully.

---

## Workflow

1. **Chop reference genome** → generates synthetic reads.
2. **Map reads** → align chopped reads to another reference genome.
3. **Variant calling** → identify differences between the genomes.
4. **Consensus calling** → update the reference with variants.

Result: an updated reference genome incorporating the detected variants.

---

## Requirements

- SLURM workload manager
- Installed tools (called from scripts):
  - `bowtie2`
  - `samtools`
  - `bcftools`
- Access to an HPC project account

---

## How to Run

1. Clone the repository and adjust directory paths in `0_genomake_main.sh` to match your project structure.

   Key variables to update inside the script:
   - `ACCOUNT_ID` → your SLURM project account
   - `DIR_WORKING` → working directory
   - `DIR_SCRIPTS` → path to pipeline scripts
   - `FASTA` → genome to chop into reads
   - `REF_FILE` → reference genome to update

2. Launch the pipeline:

   ```bash
   bash 0_genomake_main.sh
   ```
   The pipeline will:

   - Submit jobs to SLURM
   - Manage dependencies between steps
   - Write logs to slurm_output/

---

## Output

   - 1_chopped_fasta/ → synthetic reads (FASTQ)
   - 2_mapping/ → mapped & filtered BAM files
   - 3_variants/ → called & filtered BCF files
   - 4_consensus_calling/ → updated reference genome (FASTA)

---
## Contact

For **bug reports**, **comments**, or **suggestions**, please open an issue here on GitHub.

If you need to reach me directly:
**Bilal Sharif**
📧 [bilal.bioinfo@gmail.com](mailto:bilal.bioinfo@gmail.com)
