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
module load bwa
module load samtools
module load igenome-human/hg19


# Define path to the reference genome
# No need if loading the module igenome-human/hg19
# REF_GENOME="/path/to/hg19.fasta"

# Directory containing FASTQ files
FASTQ_DIR="/cluster/projects/lokgroup/rotations_students/victoria_gao/CDX_rna_fastq_trimmed"

# Directory for output
OUTPUT_DIR="/cluster/projects/lokgroup/rotations_students/victoria_gao/cdx_rna_seq_results"

# Index the reference genome (if not already indexed)
# No need if loading the module igenome-human/hg19
# bwa index $REF_GENOME

# Loop through each FASTQ file in the directory
for FASTQ_FILE in "$FASTQ_DIR"/*.fastq.gz; do
    # Extract base name of file for naming output
    BASE_NAME=$(basename "$FASTQ_FILE" .fastq.gz)

    # Define output file names
    OUTPUT_SAM="$OUTPUT_DIR/${BASE_NAME}.sam"


    # Alignment with BWA
    bwa mem $BWAINDEX $FASTQ_FILE > $OUTPUT_SAM

    echo "Created sam file: $BASE_NAME"
done


# Loop through each FASTQ file in the directory
for SAM_FILE in "$OUTPUT_SAM"/*.sam; do
    # Extract base name of file for naming output
    BASE_NAME=$(basename "$OUTPUT_SAM" .sam)

    # Define output file names
    OUTPUT_BAM="$OUTPUT_DIR/${BASE_NAME}.bam"

    # Convert SAM to BAM
    samtools view -Sb $OUTPUT_SAM > $OUTPUT_BAM

    echo "Created bam files: $BASE_NAME"
done

echo "All BWA alignments completed!"

# # Loop through each FASTQ file in the directory
# for BAM_FILE in "$FASTQ_DIR"/*.fastq.gz; do
#     # Extract base name of file for naming output
#     BASE_NAME=$(basename "$FASTQ_FILE" .fastq.gz)

#     # Define output file names
#     OUTPUT_SAM="$OUTPUT_DIR/${BASE_NAME}.sam"
#     OUTPUT_BAM="$OUTPUT_DIR/${BASE_NAME}.bam"
#     SORTED_BAM="$OUTPUT_DIR/${BASE_NAME}_sorted.bam"

#     # Alignment with BWA
#     bwa mem $BWAINDEX $FASTQ_FILE > $OUTPUT_SAM

#     # Convert SAM to BAM
#     samtools view -Sb $OUTPUT_SAM > $OUTPUT_BAM

#     # Sort the BAM file
#     samtools sort $OUTPUT_BAM -o $SORTED_BAM

#     # Index the sorted BAM file
#     samtools index $SORTED_BAM

#     echo "Processed $BASE_NAME"
# done
# echo "All BWA alignments completed!"
