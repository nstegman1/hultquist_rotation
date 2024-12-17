#!/bin/bash
#SBATCH --account=b1042  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=genomics  ### PARTITION (buyin, short, normal, etc)
#SBATCH --mem=10G
#SBATCH --time=10:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --job-name=practice_APA  ## When you run squeue -u NETID this is how you can identify the job
#SBATCH --output=outlog_final_ht1080_repac ## standard out and standard error goes to this file
#SBATCH --mail-type=ALL ## you can receive e-mail alerts from SLURM when your job begins and when your job finishes (completed, failed, etc)
#SBATCH --mail-user=natalie.stegman@northwestern.edu ## your email

module load R/4.4.0
Rscript /projects/p32549/APA_ht1080/APA_script_ht1080.R #Replace to your path to APA R script
