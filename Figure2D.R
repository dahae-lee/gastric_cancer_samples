options(stringsAsFactors = F)
library(readxl)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(gProfileR)
library(dplyr)
date = '20200423'

### Save Nitrated Peptide information

fs = list.files('~/GoogleDrive/EOGC_Nitration/', recursive = T, pattern = 'xlsx', full.names = T)

for (f in fs){
  sample = strsplit(strsplit(f, '/', fixed = T)[[1]][8], '_', fixed = T)[[1]][1]
  
  ## Load peptide data
  e = as.data.frame(read_excel(f, sheet='Peptide View'))
  samples = c('norm1', 'tumor1', 'norm2', 'tumor2')
  colnames(e)[5:8] <- samples
  
  ## Parse nitrated peptides
  e$isY45 <- ifelse(grepl('Y[45]', e$Sequence, fixed = T), 1, 0)
  e$isC29 <- ifelse(grepl('C[29]', e$Sequence, fixed = T), 1, 0)
  e$isNitrated <- ifelse(e$isC29==1 | e$isY45 == 1, 1, 0)
  
  ## Remove peptides with missing expression
  e = e[!(e$norm1==0 | e$norm2==0 | e$tumor1==0 | e$tumor2==0),]
  
  #e$Sequence에서 알파벳만 빼서 e$Sequence2에 넣어 준다.
  e$Sequence2 <- gsub("[^a-zA-Z]", "", e$Sequence)
  e$norm1 <- log2(e$norm1 + 1)
  e$norm2 <- log2(e$norm2 + 1)
  e$tumor1 <- log2(e$tumor1 + 1)
  e$tumor2 <- log2(e$tumor2 + 1)
  e$pid <- paste('p', seq(1,nrow(e)), sep='')
  
  ## Make table with peptide, protein and nitrated information
  f_blastp = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/blastp.20200420/table.blastp', sample, '20200420', 'txt.gz', sep='.')
  res = read.table(gzfile(f_blastp), header = F, col.names = c('pid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'))
  res = res[res$pident >= 100,]
  res = na.omit(res)
  m = merge(e[,c('pid', 'isNitrated')], res[,c('pid', 'sseqid')], by='pid')
  m$sseqid <- as.character(m$sseqid)
  m <- m %>% rename('sseqid' = 'gene_id')
  
  ## Save to file
  m_out = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/Nitrated_info/nitration_info', sample, date, 'txt', sep='.')
  write.table(m, m_out, sep='\t', quote = F, row.names = F, col.names = F)
}
  
  

### DEP
fs = list.files('~/GoogleDrive/gastric_cancer_samples/Inputs', pattern='raw', full.names = T)

for (f in fs){
  sid = strsplit(f, '.', fixed = T)[[1]][3]
  
  ## Sample matrix
  s <- data.frame(condition = c(0,0,1,1), batch = c(1,2,1,2))
  s$condition <- factor(s$condition)
  s$batch <- factor(s$batch)
  
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
  
  ## Add DEG change info
  data$change <- ifelse(data$padj <= 0.05, ifelse(data$lfc > 0, 'Up', 'Down'), 'None')
  data$padj <- ifelse( data$padj < 1e-15, 2e-15, data$padj)
  x_range = max(abs( data$lfc )) * 1.05
  
  ### Merge two tables
  m_in <- paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/Nitrated_info/nitration_info', sid, date, 'txt', sep='.')
  m <- read.delim(m_in, col.names = c('pid','isNitrated','gene_id'))
  a <- merge(data[,c('gene_id','protein','uniprot','change')], m, by = 'gene_id')
  a$sample_id <- sid 
  
  ## Save to file
  a_out = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/DEP_and_Nirated_info/DEP_Nitrated_info', sid, date, 'txt', sep='.')
  write.table(a, a_out, sep='\t', quote = F, row.names = F, col.names = F)
}

### Figure 2D
date = '20200423'

##Load data
dir = "~/GoogleDrive/gastric_cancer_samples/Tables/DEP_and_Nirated_info"
sample_raw <- list.files(dir)
sampleFiles <- file.path(dir, sample_raw)
tmp <- lapply(seq(1,length(sampleFiles)), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

a <- rbind_list(tmp)
colnames(a) <- c('gene_id', 'protein','uniprot', 'DEP_change', 'pid', 'isNitrated','sample_id')
a$Nitration <- ifelse(a$isNitrated==1, 'Nitrated', 'NonNitrated')

library(ggmosaic)
p <- a %>% ggplot()+
  geom_mosaic(aes(x=product(DEP_change,Nitration), fill=DEP_change))+
  xlab("Nitrated Peptide")+
  ylab("DEP Change")
p

p_out <- paste('~/GoogleDrive/gastric_cancer_samples/Figures/Fig.2D', date,'pdf', sep = ".")
ggsave(p_out, p, width = 6, height = 6)



  