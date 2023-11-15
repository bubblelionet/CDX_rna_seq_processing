#!/bin/sh

#SBATCH -t 5-00:00:00
#SBATCH --mem=60G
#SBATCH -J bwa_alignment
#SBATCH -p himem
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o bwa_alignment_output_%j.txt
#SBATCH -e bwa_alignment_error_%j.txt


# Load modules or set paths (if required)
module load trimgalore


# Directory containing FASTQ files
FASTQ_DIR="/cluster/projects/lokgroup/rotations_students/cdx_rna_seq"

# Directory for output
OUTPUT_DIR="/cluster/projects/lokgroup/rotations_students/victoria_gao/CDX_rna_fastq_trimmed"

# Loop through each FASTQ file in the directory
for FASTQ_FILE in "$FASTQ_DIR"/*.fastq.gz; do
    # Extract base name of file for naming output
    BASE_NAME=$(basename "$FASTQ_FILE" .fastq.gz)

    # Define output file names
    OUTPUT_SAM="$OUTPUT_DIR/cdx_rna_seq_results/${BASE_NAME}.sam"
    OUTPUT_BAM="$OUTPUT_DIR/cdx_rna_seq_results/${BASE_NAME}.bam"
    SORTED_BAM="$OUTPUT_DIR/cdx_rna_seq_results/${BASE_NAME}_sorted.bam"

    # Remove adapters from fastq
    trim_galore --quality 20 --gzip --fastqc $FASTQ_FILE -o $OUTPUT_DIR

    echo "Trimmed $BASE_NAME"
done

echo "All FASTQ adapters trimmed!"
