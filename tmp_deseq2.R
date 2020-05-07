options(stringsAsFactors = F)
library(readxl)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(gProfileR)
library(dplyr)
date = '20200421'

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
  data$change <- ifelse(data$padj <= 0.05, ifelse(data$lfc > 0, 'Up', 'Down'), 'None')
  data$padj <- ifelse( data$padj < 1e-15, 2e-15, data$padj)
  x_range = max(abs( data$lfc )) * 1.05
  p <- ggplot(data, aes(x = lfc, y = -log10(padj) )) + geom_point(aes(fill = change), 
                                                                  shape=21, alpha = 0.8, na.rm = T, stroke=0.25, size=2) + 
    labs(title = 'DEG analysis: KI') + 
    xlab(expression(log[2]("Celltype Tumor" / "Normal"))) + 
    ylab(expression(-log[10]("adjusted p-value"))) + 
    geom_hline(yintercept = 1.3, colour = "darkgrey") + 
    geom_vline(xintercept = 0, colour = "black") + 
    ggplot2::xlim(-x_range, x_range) + 
    scale_fill_manual(values = c("Up" = "#E64B35", 
                                 "Down" = "#3182bd", 
                                 "None" = "grey")) # change colors
  p
  
  p_out <- paste('~/GoogleDrive/gastric_cancer_samples/Figures/DEG.200421/plot_DEG', sid, date,'pdf', sep = ".")
  ggsave(p_out, p, width = 6, height = 4)
  
  ## Add GO analysis (+ save GO result to table)
  # Add gene information
  data <- data %>% rename('uniprot' = 'uniprot_ids')
  hgnc = read.delim('~/GoogleDrive/resources/hgnc/hgnc_complete_set.txt')
  data = merge(data, hgnc[,c('uniprot_ids', 'ensembl_gene_id', 'symbol', 'hgnc_id')], by='uniprot_ids')
  DEG_total_protein <- data %>% filter(padj < 0.05) %>% pull(symbol) %>% unique()
  DEG_up_protein <- data %>% filter(padj < 0.05 & change == 'Up') %>% pull(symbol) %>% unique()
  DEG_down_protein <- data %>% filter(padj < 0.05 & change == 'Down') %>% pull(symbol) %>% unique()
  DEG_None_protein <- data %>% filter(change == 'None') %>% pull(symbol) %>% unique()
  bg_proteins = data %>% pull(symbol) %>% unique()
  up_protein <- gprofiler(DEG_up_protein, 
                          organism = "hsapiens", 
                          src_filter = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"), 
                          custom_bg = bg_proteins)
  down_protein <- gprofiler(DEG_down_protein, 
                            organism = "hsapiens", 
                            src_filter = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"), 
                            custom_bg = bg_proteins)
  
  write.table(up_protein, paste('Tables/DESeq2_GO.200421/up/table.GO_up', sid, 'txt', sep='.'), sep='\t', quote = F, row.names = F, col.names = T)
  write.table(down_protein, paste('Tables/DESeq2_GO.200421/down/table.GO_down', sid, 'txt', sep='.'), sep='\t', quote = F, row.names = F, col.names = T)
  
  ## Write GO size table
  GO_size <- data.frame('Sample_ID' = sid, 
                        'GO_up_size' = nrow(up_protein),
                        'GO_down_size' = nrow(down_protein))
  
  write.table(GO_size, paste('Tables/DESeq2_GO.200421/size/table.GO_size', sid, 'txt', sep='.'), sep='\t', quote = F, row.names = F, col.names = F)
  
  ## Save DESeq2 result to table
  
  write.table(data, paste('Tables/DESeq2_200421/table.DESeq2', sid, 'txt', sep='.'), sep='\t', quote = F, row.names = F, col.names = T)
  ## Extract DEG and size in order to write size table
  DEG_size <- data.frame("Sample_ID" = sid,
                         "Total_Protein" = length(bg_proteins),
                         "DEG_total_Protein" = length(DEG_total_protein),
                         "DEG_Up_Protein" = length(DEG_up_protein),
                         "DEG_Down_Protein" = length(DEG_down_protein),
                         "DEG_None_Protein" = length(DEG_None_protein))
  
  write.table(DEG_size, paste('Tables/DEG_size.200421/table.DEG_size', sid, 'txt', sep='.'), sep='\t', quote = F, row.names = F, col.names = F)
}

## Write size table

dir = "~/GoogleDrive/gastric_cancer_samples/Tables/DEG_size.200421"
sample_raw <- list.files(dir)
sampleFiles <- file.path(dir, sample_raw)
tmp <- lapply(seq(1,length(sampleFiles)), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

dir_g = "~/GoogleDrive/gastric_cancer_samples/Tables/DESeq2_GO.200421/size"
sample_raw_g <- list.files(dir_g)
sampleFiles_g <- file.path(dir_g, sample_raw_g)
tmp_g <- lapply(seq(1,length(sampleFiles_g)), function(i){
  X <- read.delim(sampleFiles_g[i], header = F, sep = '\t')
  X})

a <- rbind_list(tmp)
colnames(a) <- c('Sample_ID', 'Total_Protein','DEG_total_Protein', 'DEG_Up_Protein', 'DEG_Down_Protein', 'DEG_None_Protein')

b <- rbind_list(tmp_g)
colnames(b) <- c('Sample_ID', 'GO_Up', 'GO_Down')

c <- merge(a, b, by = 'Sample_ID', all.x = T)

write.table(a, paste('Tables/table.DEG_and_GO_size_by_sample', date, 'txt', sep='.'), sep='\t', quote = F, row.names = F, col.names = T)





