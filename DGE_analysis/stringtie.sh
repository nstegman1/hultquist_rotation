#!/bin/bash
# This is an array job for each pairing of sample
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name="stringtie_\${SLURM_ARRAY_TASK_ID}"
#SBATCH --output=array_job_%A_%a.out  # %A is the master job ID, %a is the array task ID
#SBATCH --error=array_job_%A_%a.err   # Error log for each task
#SBATCH --time=05:00:00
#SBATCH --mem=15G       ### edit this line to reserve more/less memory
#SBATCH --nodes=3               ## how many computers/nodes do you need (no default)
#SBATCH --ntasks-per-node=3     ## how many cpus or processors do you need on per computer/node 1
#SBATCH --array=0-2
#SBATCH --mail-type=ALL ## you can receive e-mail alerts from SLURM when your job begins and when your job finishes (completed, failed, etc)
#SBATCH --mail-user=natalie.stegman@northwestern.edu ## your email

module purge
# module load stringtie I give up I am installing latest version of stringtie bc prepDE.py too old for py3
module load python/3.10.1 # I have to use an older version of python to run prepDE.py

#calls all the input files and folders I want to output to
BAM_PATH="/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files"
IN_FILES="/projects/p32549/diff_gene_exp/input.txt"
RES_PATH="/projects/p32549/diff_gene_exp/result_lists"
ANNO_PATH="/projects/p32549/diff_gene_exp/gencode.v38.annotation_renamed.gtf"
TIE_OUT_PATH="/projects/p32549/diff_gene_exp/stringtie_ht1080_out"
MERGED_GTF="/projects/p32549/diff_gene_exp/merged_gtf"
MATRIX_PATH="/projects/p32549/diff_gene_exp/count_matricies" 
STR_PATH="/projects/p32549/diff_gene_exp/stringtie/stringtie"

#these are the labels of the DGE samples 
labels=("WT_10_DMSO" "N74D_10_DMSO" "WT_10_NVP")

#these are their control file
control_files=("EOK_28.bam" "EOK_29.bam" "EOK_30.bam")
tracker=(0 3 6)

#list of the metadata files to run stringtie
txt_files=("result_list_19.txt" "result_list_22.txt" "result_list_25.txt")

#gets the 3 file names of the results.txt
F1=$(sed -n "$((tracker[SLURM_ARRAY_TASK_ID] + 1))p" $IN_FILES)
F2=$(sed -n "$((tracker[SLURM_ARRAY_TASK_ID] + 2))p" $IN_FILES)
F3=$(sed -n "$((tracker[SLURM_ARRAY_TASK_ID] + 3))p" $IN_FILES)

#this is the name of the outfiles I want them to be called
OUT_NAME1=$(echo "$F1" | cut -d '.' -f 1)
OUT_NAME2=$(echo "$F2" | cut -d '.' -f 1)
OUT_NAME3=$(echo "$F3" | cut -d '.' -f 1)

#My sanity checks so that I know I am using the correct inputs
echo "$((SLURM_ARRAY_TASK_ID))p"
echo 'I am now on '$F1 ', '$F2', and' $F3'!'
echo "I am on this comparison: ${labels[SLURM_ARRAY_TASK_ID]}"
echo "controls: ${control_files[0]}, ${control_files[1]}, ${control_files[2]}"
echo "My results_list.txt is: " $RES_PATH/${txt_files[SLURM_ARRAY_TASK_ID]}
echo "These are the base names for the stringtie output: " $OUT_NAME1' ,'$OUT_NAME2', '$OUT_NAME3
echo -e "\n"

#will be added to stringtie output
ST_OUT="_stringtie_output"

#path to where you want the .gtf file to go
GTF1=$TIE_OUT_PATH/${OUT_NAME1}${ST_OUT}/${OUT_NAME1}.gtf
GTF2=$TIE_OUT_PATH/${OUT_NAME2}${ST_OUT}/${OUT_NAME2}.gtf
GTF3=$TIE_OUT_PATH/${OUT_NAME3}${ST_OUT}/${OUT_NAME3}.gtf

#these are the mock controls that all the conditions are being compared to.
# if you want, you can run these in interactive srun and delete these lines so that it does not 
# get rerun for each job submitted.
stringtie -e -B -G $ANNO_PATH \
-o $TIE_OUT_PATH/EOK_28${ST_OUT}/EOK_28.gtf -l $TIE_OUT_PATH/EOK_28${ST_OUT}/EOK_28\
    $BAM_PATH/EOK_28.bam

stringtie -e -B -G $ANNO_PATH \
-o $TIE_OUT_PATH/EOK_29${ST_OUT}/EOK_29.gtf -l $TIE_OUT_PATH/EOK_29${ST_OUT}/EOK_29\
    $BAM_PATH/EOK_29.bam

stringtie -e -B -G $ANNO_PATH \
-o $TIE_OUT_PATH/EOK_30${ST_OUT}/EOK_30.gtf -l $TIE_OUT_PATH/EOK_30${ST_OUT}/EOK_30\
    $BAM_PATH/EOK_30.bam

# Now time to call stringtie on experimental group (infected samples) files!
$STR_PATH -e -B -G $ANNO_PATH -o $GTF1 -l $TIE_OUT_PATH/${OUT_NAME1}${ST_OUT}/${OUT_NAME1} $BAM_PATH/$F1
$STR_PATH -e -B -G $ANNO_PATH -o $GTF1 -l $TIE_OUT_PATH/${OUT_NAME1}${ST_OUT}/${OUT_NAME1} $BAM_PATH/$F1

$STR_PATH -e -B -G $ANNO_PATH -o $GTF2 -l $TIE_OUT_PATH/${OUT_NAME2}${ST_OUT}/${OUT_NAME2} $BAM_PATH/$F2
$STR_PATH -e -B -G $ANNO_PATH -o $GTF2 -l $TIE_OUT_PATH/${OUT_NAME2}${ST_OUT}/${OUT_NAME2} $BAM_PATH/$F2

$STR_PATH -e -B -G $ANNO_PATH -o $GTF3 -l $TIE_OUT_PATH/${OUT_NAME3}${ST_OUT}/${OUT_NAME3} $BAM_PATH/$F3
$STR_PATH -e -B -G $ANNO_PATH -o $GTF3 -l $TIE_OUT_PATH/${OUT_NAME3}${ST_OUT}/${OUT_NAME3} $BAM_PATH/$F3

echo -e "\n"
echo "running prepDE.py3"
#A note here I had to manually install stringtie and did not just do "module load stringtie"
echo python3 /projects/p32549/diff_gene_exp/stringtie/prepDE.py3 -i $RES_PATH/${txt_files[SLURM_ARRAY_TASK_ID]} -g \
$MATRIX_PATH/${ugh}${labels[SLURM_ARRAY_TASK_ID]}.csv

python3 /projects/p32549/diff_gene_exp/stringtie/prepDE.py3 -i $RES_PATH/${txt_files[SLURM_ARRAY_TASK_ID]} -g \
$MATRIX_PATH/${ugh}${labels[SLURM_ARRAY_TASK_ID]}.csv
