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
- [Shared Gene Lists](#Shared-Gene-Lists)
  * [Dependencies](#Dependencies)
  * [Executing Program](#Executing-Program)
- [Data Visualizations](#Data-Visualizations)
  * [Dependencies](#Dependencies)
  * [Executing Program](#Executing-Program)
  

<!-- toc -->
## Description

This project consists of:
* RNAseq Pipeline
* APA Analysis
* Creation of shared gene lists
* Data Visualization

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
* Change paths and create folders where you want the output to be stored
* You may have to manually install StringTie from its GitHub repo
* I did not submit this as a job, but still ran it through Quest

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

## Shared Gene Lists

### Dependencies

### Installing

* How/where to download your program
* Any modifications needed to be made to files/folders

### Executing program
```
print('hello world')
```
## Data Visualizations

### Dependencies

### Installing

* How/where to download your program
* Any modifications needed to be made to files/folders

### Executing program
```
print('hello world')
```


