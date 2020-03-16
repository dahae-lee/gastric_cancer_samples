options(stringsAsFactors = F)
library(readxl)
library(tidyverse)
library(DESeq2)

## list files of summary data
fs = list.files('~/GoogleDrive/EOGC_Nitration/', recursive = T, pattern = 'xlsx', full.names = T)

s <- data.frame(condition = c(0,0,1,1), batch = c(1,2,1,2))
s$condition <- factor(s$condition)
s$batch <- factor(s$batch)

res = data.frame()
for (f in fs){
  if ('Protein View' %in% readxl::excel_sheets(f)){
    d = read_excel(f, sheet = 'Protein View')
    sid = strsplit(strsplit(f, '/', fixed = T)[[1]][8], '_', fixed = T)[[1]][1]
    res1 = data.frame(sid = sid, 
                      file = f,
                      n_protein = nrow(d),
                      total_norm1 = d %>% select(`iTRAQ-114`) %>% sum/100000,
                      total_tumor1 = d %>% select(`iTRAQ-115`) %>% sum/100000,
                      total_norm2 = d %>% select(`iTRAQ-116`) %>% sum/100000,
                      total_tumor2 = d %>% select(`iTRAQ-117`) %>% sum/100000)
    res1$total_size = (res1$total_norm1 + res1$total_tumor1 + 
                         res1$total_norm2 + res1$total_tumor2)
    res = rbind.data.frame(res, res1)
    
    ## Create a normalized count matrix
    d = d %>% mutate(norm1 = as.integer(`iTRAQ-114`/100000),
                     tumor1 = as.integer(`iTRAQ-115`/100000),
                     norm2 = as.integer(`iTRAQ-116`/100000),
                     tumor2 = as.integer(`iTRAQ-117`/100000))
    
    d1 <- d %>% select(norm1, norm2, tumor1, tumor2)
    rownames(d1) <- d$Accession
    
    dds <- DESeqDataSetFromMatrix(countData = d1,
                                  colData = s,
                                  design= ~ batch + condition)
    
    dds <- estimateSizeFactors(dds)
    m = as.data.frame(counts(dds, normalize=T)) %>% 
      mutate(gene_id = rownames(dds)) %>% select(gene_id, norm1, norm2, tumor1, tumor2)
    
    d2 = d1 %>% mutate(gene_id = rownames(dds)) %>% select(gene_id, norm1, norm2, tumor1, tumor2)
    
    write.table(d2, paste('Inputs/table', 'raw', sid, 'txt', sep='.'), sep='\t', quote = F, row.names = F, col.names = T)
    write.table(m, paste('Inputs/table', 'norm', sid, 'txt', sep='.'), sep='\t', quote = F, row.names = F, col.names = T)
  } else {
    readxl::excel_sheets(f)
  }
}

