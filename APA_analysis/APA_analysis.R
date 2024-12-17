# Load necessary packages:
library('REPAC')
library(GenomicRanges)
library(tidyverse)
library(limma)
library(edgeR)
data('hg38_pa')

run_repac <- function (in_files, out_path, separate) {
    # in_files is a vector of the filenames that I want to compare
    # out_path tells computer where to output the final csv, and 
    #separate just tells computer where to split the file path to shorten the colname. 

    print(in_files) # sanity check
    print(out_path) # sanity check

    # in_files is a vector of the filenames that I want to compare. This list of filenames coorindates
    # with the colnames and replicate #'s as shown below. 
    pheno_WT_S <- data.frame(groups=rep(c('Uninfected','Uninfected','Uninfected','Infected','Infected','Infected')),exp=rep(c('R1','R2','R3','R1','R2','R3')))
    
    pheno_WT_S$Sample<-in_files
    pheno_WT_S<-pheno_WT_S%>%separate(Sample,c(NA,"Sample"),sep = separate)
    print(pheno_WT_S)
  
    # This is the modified command that I changed to PE (isPairedEnd = T) and unstranded (strandSpecific = 0).
    # everything else is the same:
    se_WT_S_R2 <- rse_from_bam(in_files, hg38_pa, pheno_WT_S, strandSpecific = 0,  allowMultiOverlap = F,
                            fraction = F, countMultiMappingReads=F, isPairedEnd = T, nthreads = 8,
                            countReadPairs = T,  requireBothEndsMapped = F)

    # The rest of this function are the same as original script.
    keep_WT_S_R2 <- edgeR::filterByExpr(se_WT_S_R2, min.count = 20)
    se_WT_S_R2 <- se_WT_S_R2[keep_WT_S_R2,]
    se_WT_S_3UTR_R2 <- se_WT_S_R2[se_WT_S_R2@rowRanges$annotation == '3UTR',] 

    ## fit model and test for DPU
    groups_WT_S<-c(pheno_WT_S$groups)

    exp_WT_S<- c(pheno_WT_S$exp)
    dMat_WT_S <- model.matrix(~ 0+groups_WT_S+exp_WT_S) # this is missing R1

    colnames(dMat_WT_S) <-  gsub('groups_WT_S', '',colnames(dMat_WT_S))
 

    cMat_WT_S <- limma::makeContrasts(
    levels=colnames(dMat_WT_S),
    InfectedvsUninfected= Infected - Uninfected
    )

    # #Test for DPU
    n <- table(SummarizedExperiment::rowData(se_WT_S_3UTR_R2)$gene_name) > 1
    keep <- names(n[n == TRUE])
    se <- se_WT_S_3UTR_R2[SummarizedExperiment::rowData(se_WT_S_3UTR_R2)$gene_name %in% keep,]
    genes <- unique(SummarizedExperiment::rowData(se_WT_S_3UTR_R2)$gene_name)

    results_WT_S_3UTR_R2 <- fit_repac(se_WT_S_3UTR_R2, dMat_WT_S, cMat_WT_S) 

    # # Resolve how PAS are compared
    results_WT_S_3UTR_R2 <- resolve_comparisons(results_WT_S_3UTR_R2)
    write.csv(results_WT_S_3UTR_R2, out_path)

}

# 7 comparisons for ht1080:

# 4 from MOI 01 bams:
uninfected_01 <- c('/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_13.bam',
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_14.bam',
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_15.bam') #EOK_13,14,15
                    
# Infected WT DMSO vs Uninfected EOK_1,2,3
WT_DMSO_01_in <- c(uninfected_01,
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_1.bam',
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_2.bam',
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_3.bam')

# Infected N74D DMSO vs Uninfected​ EOK_4,5,6
N74D_DMSO_01_in <- c(uninfected_01,
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_4.bam',
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_5.bam',
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_6.bam')

# Infected WT NVP vs Uninfected EOK_7,8,9
WT_NVP_01_in <- c(uninfected_01,
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_7.bam',
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_8.bam',
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_9.bam')

# Infected N74D NVP vs Uninfected EOK_10,11,12
N74D_NVP_01_in <- c(uninfected_01,
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_10.bam',
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_11.bam',
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_12.bam')


# 3 from MOI 10 bams:
uninfected_10 <- c('/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_28.bam',
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_29.bam',
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_30.bam') #EOK_28,29,30

# Infected WT DMSO vs Uninfected EOK_19,20,21
WT_DMSO_10_in <- c(uninfected_10,
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_19.bam',
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_20.bam',
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_21.bam')

# Infected N74D DMSO vs Uninfected​ EOK_22,23,24
N74D_DMSO_10_in <- c(uninfected_10,
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_22.bam',
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_23.bam',
                    '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_24.bam') # /projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/

# Infected WT NVP vs Uninfected EOK_25,26,27
WT_NVP_10_in <- c(uninfected_10,
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_25.bam',
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_26.bam',
                '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_bam_files/EOK_27.bam')

# this is a list of lists of each comparison I wanna make
#this runs all 7 comparisons:
#lol <- list(WT_DMSO_01_in, N74D_DMSO_01_in, WT_NVP_01_in, N74D_NVP_01_in, WT_DMSO_10_in, N74D_DMSO_10_in, WT_NVP_10_in)

#this runs just the 10 MOI comparisons
lol <- list(WT_DMSO_10_in, N74D_DMSO_10_in, WT_NVP_10_in)
print(lol)

# this is a list of out paths where the repac csv output will go for the 7 comparisons
#out_files <- c('WT_DMSO_01','N74D_DMSO_01','WT_NVP_01','N74D_NVP_01','WT_DMSO_10','N74D_DMSO_10','WT_NVP_10')

# this is a list of out paths where the repac csv output will go for the 10 MOI comparisons
out_files <- c('WT_DMSO_10','N74D_DMSO_10','WT_NVP_10')

#Where I want the repac csv output to be exported
out_path <- '/projects/p32549/APA_ht1080/repac_out_ht1080/'

#this runs through all the comparisons I set up in the lol list
for (i in 1:length(lol)) {
    print('run repac')
    print(lol[i])
    #outtie is the full path to where the outfile of repac should be saved and what 
    #the file name should be
    outtie <- paste0(out_path,out_files[i],'.csv')
    print(c('run_repac',lol[[i]],outtie))
    run_repac(lol[[i]],outtie,'ht1080_bam_files/')
    print('DONE WITH THIS ONE\n')
}


