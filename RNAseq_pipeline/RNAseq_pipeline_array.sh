#!/bin/bash
# This is an array job for each pairing of sample. MAKE SURE you are
# in dir that you want the outlogs to output to.
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name="trimmhisatsamtools_ht1080_\${SLURM_ARRAY_TASK_ID}"
#SBATCH --output=array_job_%A_%a.out  # %A is the master job ID, %a is the array task ID
#SBATCH --error=array_job_%A_%a.err   # Error log for each task
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=10G  # Memory per task/CPU
#SBATCH --nodes=4               ## how many computers/nodes do you need (no default)
#SBATCH --ntasks-per-node=4     ## how many cpus or processors do you need on per computer/node 1
#SBATCH --array=0-3     ## there are 30 samples (2 fastq files per sample because paired-end) only run 4 tasks at a time
#SBATCH --mail-type=ALL ## you can receive e-mail alerts from SLURM when your job begins and when your job finishes (completed, failed, etc)
#SBATCH --mail-user=natalie.stegman@northwestern.edu ## your email

#loading the tools I need to submit my job
module purge #gets rid of previously loaded modules 
module load java/jdk11.0.10 #loads java in because trimmomatic uses is
module load trimmomatic/0.39 #loads in trimmmomatic
module load python #loads in python
module load hisat2/2.1.0 #loads in hisat2
module load samtools/1.10.1 #loads in samtools

#trimmomatic required me to call its full path so I defined that here
TRIMMOMATIC_JAR="/software/trimmomatic/0.39/trimmomatic-0.39.jar" # path to trimmomatic 

#this is the path to the adapter files for the trimmomatic command
ADAPTER_PATH="/software/trimmomatic/0.39/adapters" # path to the adapter files for trimmomatic command

# These are the paths to where I want my data to be outputted to. Make sure you 
# have already created these folders so this pipeline does not produce an error
# saying the directory you called isn't found
PAIRED_PATH="/projects/p32549/RNA_seq_analysis_ht1080/ht1080_trimmed_fastq_paired" # outpath for trimmed paired fastqs
UNPAIRED_PATH="/projects/p32549/RNA_seq_analysis_ht1080/ht1080_trimmed_fastq_unpaired" # outpath for trimmed unpaired fastqs
PATH_RAW_READS="/projects/p32549/RNA_seq_analysis_ht1080/ht1080_raw_fastqs" # base path to where unzipped raw reads where
SAM_PATH="/projects/p32549/RNA_seq_analysis_ht1080/ht1080_sam_files" # outpath to put sam files
BAM_PATH="/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files" # outpath to put bam files

# Before I even TOUCHED this script, I spearated the filenames for the R1 and R2 fastqs, sorted them,
# and wrote them out to either a txt file of R1 files or a txt file of R2 files. This made it easier
# for me to iterate through each paring of fastq files per job.
R1_FILES='/projects/p32549/RNA_seq_analysis_ht1080/R1_ht1080_lil.txt'
R2_FILES='/projects/p32549/RNA_seq_analysis_ht1080/R2_ht1080_lil.txt'

# This just isolates one file from the R1_FILE or R2_FILE path based on what job # computer is on.
R1_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $R1_FILES)
R2_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $R2_FILES)
SAM_BAM=$(echo "$R1_FILE" | cut -d'_' -f1-2)

# the echo commands you see are sanity checks so that I know where the computer is at 
# in the outlog
echo 'I am now on '$R1_FILE 'and '$R2_FILE'! and my base name is '$SAM_BAM

# Run Trimmomatic on the zipped fastq files. No need to unzip.
java -jar $TRIMMOMATIC_JAR PE -phred33 -threads 32 $PATH_RAW_READS/$R1_FILE $PATH_RAW_READS/$R2_FILE \
  $PAIRED_PATH/$R1_FILE.fq.gz $UNPAIRED_PATH/$R1_FILE.fq.gz \
  $PAIRED_PATH/$R2_FILE.fq.gz $UNPAIRED_PATH/$R2_FILE.fq.gz \
  ILLUMINACLIP:$ADAPTER_PATH/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
echo 'I am finished runnning trimmmomatic!'

# Run hisat2
# HISAT2: again, only care about trimmed paired reads
hisat2 --dta -q -p 32 -x /projects/p32549/RNA_seq_analysis/human_index/grch38_tran/genome_tran -1 $PAIRED_PATH/$R1_FILE.fq.gz -2 $PAIRED_PATH/$R2_FILE.fq.gz -S $SAM_PATH/$SAM_BAM.sam
echo 'I am finished running hisat2!'

# Run samtools
#samtools view -Su $SAM_PATH/$SAM_BAM.sam | samtools sort -o $BAM_PATH/$SAM_BAM.bam
samtools view -Su $SAM_PATH/$SAM_BAM.sam | samtools sort -@ 2 -m 10G -o $BAM_PATH/$SAM_BAM.bam
echo 'I am finished running samtools!'
