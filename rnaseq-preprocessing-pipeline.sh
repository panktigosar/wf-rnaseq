#!/bin/bash

# For requexting the node on the server: Asking for 8 threadsand 1 CPU
# #bsub -W 8:00 -R "rusage[mem=5G] span[hosts=1]" -n 8 -q short "./rnaseq-preprocessing-pipeline.sh"

# the command run while running this file
# #nohup ./rnaseq-preprocessing-pipeline.sh > output.log 2>&1 &
# the nohup and disown commands allow the script to keep running in background even after you exited the shell or the terminal. 
# The output log is to record the log of the events

email="pankti.gosar@gmail.com"

# Function to send an email for job termination
send_email_on_failure() {
    echo "The pipeline was terminated unexpectedly." | mail -s "Job Failed" "$email"
}
# Set up trap †o catch any termination signals
trap 'send_email_on_failure' ERR EXIT

# Directories
basepath="path to experiment directory"
INPUT_DIR="$basepath/rawData"
OUTPUT_DIR="$basepath/output_files"
GENOME_DIR="path to genome_dir/"
ADAPTOR="$basepath/data_pre-processing/adaptor.fa"
ANNOTATION="path to annotation files"
out_log="pipeline$(date +'%d-%B-%Y').log"

# Create output directories
mkdir -p $OUTPUT_DIR/fastqc
mkdir -p $OUTPUT_DIR/fastqc-trimmed
mkdir -p $OUTPUT_DIR/trimmed
mkdir -p $OUTPUT_DIR/alignments
mkdir -p $OUTPUT_DIR/reports

module load star/2.7.10a
module load fastqc/0.11.9
module load trimmomatic/0.32
module load samtools/1.16.1
module load subread/1.6.2


# Tools
FASTQC="fastqc"
TRIMMOMATIC="/share/pkg/trimmomatic/0.32/trimmomatic-0.32.jar"
STAR="/share/pkg/star/2.7.10a/bin/STAR"

# Loop through all the files in the input directory 
for file in $INPUT_DIR/*_R1_001.fastq.gz; do 
    base=$(basename $file "_R1_001.fastq.gz")
    R1=${INPUT_DIR}/${base}_R1_001.fastq.gz
    R2=${INPUT_DIR}/${base}_R2_001.fastq.gz
    echo "Processing sample $base..."

    # Step 1: Quality check
    echo "Step 1: Running FastQC on raw data..."
    $FASTQC $R1 $R2 -o $OUTPUT_DIR/fastqc

    # Step 2: Trimming
    echo "Step 2: Trimming the adaptors and low quality bases..."
    TRIMMED_R1=$OUTPUT_DIR/trimmed/${base}_R1_trimmed_paired.fastq.gz
    TRIMMED_R2=$OUTPUT_DIR/trimmed/${base}_R2_trimmed_paired.fastq.gz
    UNPAIRED_R1=$OUTPUT_DIR/trimmed/${base}_R1_trimmed_unpaired.fastq.gz
    UNPAIRED_R2=$OUTPUT_DIR/trimmed/${base}_R2_trimmed_unpaired.fastq.gz

    java -jar $TRIMMOMATIC PE -phred33 -threads 4 \
    $R1 $R2 \
    $TRIMMED_R1 $UNPAIRED_R1 $TRIMMED_R2 $UNPAIRED_R2 \
    ILLUMINACLIP:${ADAPTOR}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Step 3: Post-trimming Quality Check
    echo "Step 3: Running FastQC on trimmed data..."
    $FASTQC $TRIMMED_R1 $TRIMMED_R2 -o $OUTPUT_DIR/fastqc-trimmed

    # Step 4: Alignment
    echo "Generating the genome files" 

    # STAR --runThreadN 8 \
    #  --runMode genomeGenerate \
    #  --genomeDir $GENOME_DIR \
    #  --genomeFastaFiles genome.fa \
    #  --sjdbGTFfile $ANNOTATION \
    #  --sjdbOverhang 100

    echo "Step 4: Aligning reads to the reference genome..."

    $STAR --runThreadN 4 --genomeDir $GENOME_DIR \
    --readFilesIn $TRIMMED_R1 $TRIMMED_R2 --readFilesCommand zcat \
    --outFileNamePrefix $OUTPUT_DIR/alignments/${base}_ \
    --outSAMtype BAM SortedByCoordinate

    ALIGNMENT=$OUTPUT_DIR/alignments/${base}_Aligned.sortedByCoord.out.bam

    # Step 5: Post-alignment QC
    echo "Step 5: Generating alignment statistics..."
    samtools flagstat $ALIGNMENT > $OUTPUT_DIR/reports/${base}_alignment_report.txt

    # Step 6: Consolidated Report
    echo "Generating consolidated report..."
    cat $OUTPUT_DIR/reports/*_alignment_report.txt > $OUTPUT_DIR/reports/consolidated_report.txt

    # Step 7: creating count matrix:
    featureCounts -T 4 -p -s 2 -t exon -g gene_name -a $ANNOTATION -o $OUTPUT_DIR/endo_grch38v107_counts.txt $OUTPUT_DIR/alignments/*out.bam

    echo "All samples processed. Consolidated report is available at $OUTPUT_DIR/reports/consolidated_report.txt."

done

# If successful, clear the failure trap and send success email
trap - ERR EXIT
echo "The pipeline ran successfully."

echo "The logs have been written to $out_log | mail -s "Job complete" "$email"
