#!/bin/sh

#SBATCH -t 5-00:00:00
#SBATCH --mem=60G
#SBATCH -J bwa_alignment
#SBATCH -p himem
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
FASTQ_FILE="/cluster/projects/lokgroup/rotations_students/victoria_gao/CDX10.rep1.total.RNA_S3_L002_R1_001_trimmed.fq"

# Directory for output
OUTPUT_DIR="/cluster/projects/lokgroup/rotations_students/victoria_gao/cdx_rna_seq_results"
OUTPUT_SAM="/cluster/projects/lokgroup/rotations_students/victoria_gao/cdx_rna_seq_results/test_output.sam"
OUTPUT_BAM="/cluster/projects/lokgroup/rotations_students/victoria_gao/cdx_rna_seq_results/test_output.bam"
SORTED_BAM="/cluster/projects/lokgroup/rotations_students/victoria_gao/cdx_rna_seq_results/test_sorted_output.bam"


# Alignment with BWA
bwa mem $BWAINDEX $FASTQ_FILE > $OUTPUT_SAM

# Convert SAM to BAM
samtools view -Sb $OUTPUT_SAM > $OUTPUT_BAM

# Sort the BAM file
samtools sort $OUTPUT_BAM -o $SORTED_BAM

# Index the sorted BAM file
samtools index $SORTED_BAM

echo "BWA alignment process completed for CDX10 test."




