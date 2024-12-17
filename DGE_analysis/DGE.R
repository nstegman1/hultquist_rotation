# RUN ONCE
# Install Bioconductor packages
# BiocManager::install(c(
#   "DESeq2", 
#   "apeglm", 
#   "EnhancedVolcano", 
#   "clusterProfiler", 
#   "org.Hs.eg.db", 
#   "ReactomePA", 
#   "AnnotationDbi"
# ))

# # Install CRAN packages
# install.packages(c("dplyr", "pheatmap"))

# if running in quest, run in terminalk before
# module load R/4.4.0
# and 
# Rscript /path/to/script.R

library(DESeq2)
library(apeglm)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(ReactomePA)
library(pheatmap)
library(AnnotationDbi)

run_DGE <-  function(colData_path, count_matrix_path,title) {
  print(colData_path)
  print(count_matrix_path)
  countData <- read.csv(count_matrix_path,row.names = 1, check.names = FALSE)
  colData <- read.csv(colData_path,row.names = 1)

  new_rn <- c()
  for (rn in row.names(countData)){
    tmp <- strsplit(rn, "|", fixed = TRUE)[[1]][2]
    new_rn <- c(new_rn, tmp)
  }
  countData<-countData[!duplicated(new_rn),]

  #gets rid of duplicates
  new_rn <- new_rn[!duplicated(new_rn)]
  todrop <- which(is.na(new_rn))
  print(todrop)
  new_rn <- new_rn[!is.na(new_rn)]
  countData <- countData[-todrop,]

  #resets the row names 
  row.names(countData) <- new_rn

  #calls DeSeq2 to get gene counts matrix
  dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design= ~ condition)

  dds <- DESeq(dds)
  res <- results(dds, contrast = c('condition','Treatment','Control'))
  summary(res)
  pval_threshold <- 0.05
  fc_threshold <- 1

  total_genes <- nrow(res)
  significant_genes <- sum(res$padj < pval_threshold, na.rm = TRUE)
  upregulated_genes <- sum(res$log2FoldChange > 1 & res$padj < pval_threshold, na.rm = TRUE)
  downregulated_genes <- sum(res$log2FoldChange < -1 & res$padj < pval_threshold, na.rm = TRUE)
  sig_genes <- upregulated_genes + downregulated_genes

  #time to plot volcano plot
  plot <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    FCcutoff = 0.00,
    pCutoff = 0.05,
    xlim = c(-4,4),
    legendLabels = c("Not Significant", "Log2 FC", "p-value", 
    paste0('Upregulated gene count: ',upregulated_genes,'\nDownregulated gene count: ',downregulated_genes, '\ntotal: ',sig_genes)),
    legendPosition = 'right',
    legendLabSize = 10,          # Font size for legend labels
    legendIconSize = 8,         # Size of the legend icons
    title = paste0(title,' no log2FC Threshold'))

    #exports the plot to where you want it to be
    ggsave(paste0("/projects/p32549/diff_gene_exp/no_log2_FC_line_plots/",title,'.pdf'), plot=plot,width = 15, height = 10, dpi = 300)
}

#paths to the metadata
col_data_paths <- c("/projects/p32549/diff_gene_exp/col_data/coldata_19.csv",
                  "/projects/p32549/diff_gene_exp/col_data/coldata_22.csv",
                  "/projects/p32549/diff_gene_exp/col_data/coldata_25.csv")

#paths to the count matricies from stringtie
matrix_path <- c("/projects/p32549/diff_gene_exp/count_matricies/count_matrix_WT_10_DMSO.csv",
                "/projects/p32549/diff_gene_exp/count_matricies/count_matrix_N74D_10_DMSO.csv",
                "/projects/p32549/diff_gene_exp/count_matricies/count_matrix_WT_10_NVP.csv")

for (i in 1:3) {
  base_name <- basename(matrix_path[i])            # Get the file name
  clean_name <- sub("\\.csv$", "", base_name) # Remove the .csv extension

  # Extract the desired substring
  result <- sub("count_matrix_", "", clean_name) # Remove the prefix
  result
  run_DGE(col_data_paths[i], matrix_path[i],result) #run the function
}

