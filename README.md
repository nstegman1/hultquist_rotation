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
  * [Installing](#Installing)
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

* rename paths to the files that you want to compare in the R script
* Any modifications needed to be made to files/folders

### Executing program
```
print('hello world')
```
## Differential Gene Expression Analysis

### Dependencies

### Installing

* How/where to download your program
* Any modifications needed to be made to files/folders

### Executing program
```
print('hello world')
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


