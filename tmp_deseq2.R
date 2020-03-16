options(stringsAsFactors = F)
library(readxl)
library(tidyverse)
library(DESeq2)

## Sample matrix
s <- data.frame(condition = c(0,0,1,1), batch = c(1,2,1,2))
s$condition <- factor(s$condition)
s$batch <- factor(s$batch)

## Load data
fs = list.files('Inputs', pattern='raw', full.names = T)

for (f in fs){
  print (f)
  sid = strsplit(f, '.', fixed = T)[[1]][3]
  
  ## Load data
  d = read.delim(f)
  rownames(d) <- d$gene_id
  d$gene_id <- NULL
  
  ## Create DESeq2
  dds <- DESeqDataSetFromMatrix(countData = d,
                                colData = s,
                                design= ~ batch + condition)
  
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  
  dds <- DESeq(dds)
  resultsNames(dds) 
  res <- lfcShrink(dds, coef="condition_1_vs_0", type='apeglm')
  
  data = data.frame(gene_id = rownames(res), 
                    baseMean = res$baseMean,
                    lfc = res$log2FoldChange,
                    pvalue = res$pvalue,
                    padj = res$padj)
  
  ## Split multiple gene names into one row
  data <- data %>% 
    mutate(gene_id = strsplit(as.character(gene_id), ",")) %>% 
    unnest(gene_id)
  
  data$protein <- gsub('_HUMAN', '', 
                       do.call(rbind.data.frame, strsplit(data$gene_id, '|', fixed = T))[[3]])
  
  data$uniprot <- gsub('_HUMAN', '', 
                       do.call(rbind.data.frame, strsplit(data$gene_id, '|', fixed = T))[[2]])
  
  ## Add Volcano (+ save volcano to pdf)
  
  ## Add GO analysis (+ save GO result to table)
  
  ## Save DESeq2 result to table
}






