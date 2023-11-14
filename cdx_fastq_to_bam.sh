#!/bin/sh

#SBATCH -t 5-00:00:00
#SBATCH --mem=16G
#SBATCH -J bwa_alignment
#SBATCH -p all
#SBATCH -c = 1
#SBATCH -N 1
#SBATCH -o bwa_alignment_output_%j.txt
#SBATCH -e bwa_alignment_error_%j.txt


# Load modules or set paths (if required)
# module load bwa
# module load samtools

# Define path to the reference genome
# REF_GENOME="/path/to/hg19.fasta"

# Directory containing FASTQ files
FASTQ_DIR="/cluster/projects/lokgroup/rotations_students/cdx_rna_seq"

# Directory for output
OUTPUT_DIR="/cluster/projects/lokgroup/rotations_students/victoria_gao/cdx_rna_seq"

# Index the reference genome (if not already indexed)
# bwa index $REF_GENOME

# Loop through each FASTQ file in the directory
for FASTQ_FILE in "$FASTQ_DIR"/*.fastq.gz; do
    # Extract base name of file for naming output
    BASE_NAME=$(basename "$FASTQ_FILE" .fastq.gz)

    # Define output file names
    OUTPUT_SAM="$OUTPUT_DIR/${BASE_NAME}.sam"
    OUTPUT_BAM="$OUTPUT_DIR/${BASE_NAME}.bam"
    SORTED_BAM="$OUTPUT_DIR/${BASE_NAME}_sorted.bam"

    # Alignment with BWA
    bwa mem $BWAINDEX $FASTQ_FILE > $OUTPUT_SAM

    # Convert SAM to BAM
    samtools view -Sb $OUTPUT_SAM > $OUTPUT_BAM

    # Sort the BAM file
    samtools sort $OUTPUT_BAM -o $SORTED_BAM

    # Index the sorted BAM file
    samtools index $SORTED_BAM

    echo "Processed $BASE_NAME"
done

echo "All BWA alignments completed."
