# Hultquist Rotation Scripts

This is the compiled code I made for my Hultquist lab rotation project.

## Table of Contents

- [Description](#Description)
- [RNAseq Pipeline](#RNAseq-Pipeline)
  * [Dependencies](#Dependencies)
  * [Before Running](#Before-Running)
  * [Executing Program](#Executing-Program)
- [APA Analysis](#APA-Analysis)
  * [Dependencies](#Dependencies)
  * [Before Runninng](#Before-Running)
  * [Executing Program](#Executing-Program)
- [Differential Gene Expression Analysis](#Differential-Gene-Expression-Analysis)
  * [Dependencies](#Dependencies)
  * [Before Runninng](#Before-Running)
  * [Executing Program](#Executing-Program)
- [Shared Gene Lists and Correlation Plots](#Shared-Gene-Lists-and-Correlation-Plots)
  * [Dependencies](#Dependencies)
  * [Before Runninng](#Before-Running)
  * [Executing Program](#Executing-Program)
- [Other Data Visualizations](#Other-Data-Visualizations)
  * [Dependencies](#Dependencies)
  * [Before Runninng](#Before-Running)
  * [Executing Program](#Executing-Program)
- [Parse TargetScan](#Parse-TargetScan)
  * [Dependencies](#Dependencies)
  * [Before Runninng](#Before-Running)
  * [Executing Program](#Executing-Program)
  

<!-- toc -->
## Description

This project consists of:
* RNAseq Pipeline
* APA Analysis
* Shared gene lists and correlation plots
* Other data visualizations

The RNAseq, APA, and DGE analyses were performed on Quest; however, the data visualization and creation of shared gene lists were done locally on my Mac machine.

## RNAseq Pipeline

Preprocesses raw bulk RNAseq data

### Dependencies

* Trimmomatic
* Hisat2
* Human Genome transcript: https://daehwankimlab.github.io/hisat2/download/#h-sapiens
* Samtools
* Python
* Bash

### Before Running

* Make the folders you want each tool to output to (a folder to store trimmed fastqs, a folder for .sam files, and a folder for .bam files).
* Rename the paths in the script to yours
* Put all raw paired-end fastq files into one folder 
* Create a text file with all the R1 files and one with all the corresponding R2 files
  + Use .py script provided
* If you are having issues with Trimmomatic, call the full path to the tool (shown in the bash file)
* This is a job array, so adjust memory and number of nodes based on your needs

### Executing program

* ssh into Quest
 
```
sbatch /path/to/bash/RNAseq_array.sh
```
## APA Analysis

Runs REPAC on the preprocessed RNAseq data

### Dependencies

* REPAC (https://github.com/eddieimada/REPAC/tree/main)
* R
  
### Before Running

* Rename paths to the bam files that you want to compare in the R script
* Make folders where you want the REPAC output to be saved to
* There is the R file that runs REPAC and a bash file that submits the job

### Executing program
```
sbatch /path/to/bash/run_APA.sh
```
## Differential Gene Expression Analysis

Runs DGE on the preprocessed RNAseq data

### Dependencies

* StringTie (https://github.com/gpertea/stringtie)
* DeSeq2
* R
* human .gtf file (https://www.gencodegenes.org/human/release_38.html)

### Before Running

* Make a text file containing the metadata of the samples you are comparing (3 infected replicates vs 3 uninfected replicates) (result_list_ex.txt)
* Make a text file for each comparison being tested containing the name of each replicate and the path to that replicate's gtf file (coldata_ex.txt)
* Change paths in scripts to yours and create folders where you want the output to be stored
* You may have to manually install StringTie from its GitHub repo
* I did not submit this as a job, but I still ran it through Quest

### Executing program

#### Run StringTie
```
sbatch /path/to/bash/string_tie.sh
```
#### Run DeSeq2
```
module load R/4.4.0
Rscript /path/to/R/Script/DGE.R
```

## Shared Gene Lists and Correlation Plots

This finds shared genes between:
 * REPAC output from one comparison vs REPAC output from another comparison
 * REPAC output from one comparison vs DGE output from the same comparison
   
This also makes correlation plots between the two variables being compared using their fold changes

### Dependencies

* Python
* Pandas
* REPAC and DGE outputs you want to compare

### Before Running

*  Change the paths in shared_genes.py to where your files are located
*  Change the paths to where you want the plots and lists to be saved to

### Executing program

* I ran this script in the Python IDE Spyder

## Other Data Visualizations

correlation_w_labels.py: Gets the correlation plots but with gene counts of each quadrant

volcano_plot.py: Gets volcano plots from either the DGE output or REPAC output

### Dependencies

* Python
* Pandas
* Matplotlib
* REPAC and/or DGE output

### Before Running

*  Change the paths to the input data and the paths to where you want the plots to be exported to

### Executing program

* I ran this script in the Python IDE Spyder
  
## Parse TargetScan

This script summarizes the TargetScan outputs of a gene list

### Dependencies

* Python
* Pandas

### Before Running

* Run TargetScan (https://www.targetscan.org/vert_80/) through each gene in the gene list of interest
* Save the TargetScan output as a .csv file, and rename file as the gene name.
* put all these csv files into a folder
* Change the input and output path in the python script

* How/where to download your program
* Any modifications needed to be made to files/folders

### Executing program

* I ran this script in the Python IDE Spyder
