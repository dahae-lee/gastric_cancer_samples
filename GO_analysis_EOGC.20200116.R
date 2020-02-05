date = '20200116'
options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(tidyr)
library(dplyr)

# f = '~/GoogleDrive/gastric_cancer_samples/Tables/DEP_blastp.20200116./table.DEP_blastp.N13T236.20200116.txt'
f = args[1]
sample = strsplit(f, '.', fixed = T)[[1]][5]

## Load data
m = as.data.frame(read.delim(f))
bg_proteins = m %>% pull(symbol) %>% unique()

## filter the ones which have overexpressed or underexpressed
DEP_up_protein_Nitrated <- m %>% filter(adj.P.Val < 0.05 & t > 0 & isNitrated == '1') %>% pull(symbol) %>% unique()
DEP_down_protein_Nitrated <- m %>% filter(adj.P.Val < 0.05 & t < 0 & isNitrated == '1') %>% pull(symbol) %>% unique()
DEP_up_protein_NonNitrated <- m %>% filter(adj.P.Val < 0.05 & t > 0 & isNitrated == '0') %>% pull(symbol) %>% unique()
DEP_down_protein_NonNitrated <- m %>% filter(adj.P.Val < 0.05 & t < 0 & isNitrated == '0') %>% pull(symbol) %>% unique()

## run gprofiler
library(gProfileR)
up_nitrated <- gprofiler(DEP_up_protein_Nitrated, organism = "hsapiens", src_filter = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"), custom_bg = bg_proteins)
down_nitrated <- gprofiler(DEP_down_protein_Nitrated, organism = "hsapiens", src_filter = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"), custom_bg = bg_proteins)
up_nonnitrated <- gprofiler(DEP_up_protein_NonNitrated, organism = "hsapiens", src_filter = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"), custom_bg = bg_proteins)
down_nonnitrated <- gprofiler(DEP_down_protein_NonNitrated, organism = "hsapiens", src_filter = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"), custom_bg = bg_proteins)

## Save in file
up_nit_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO.20200116./up_nitrated/up_nitrated", sample, 'txt', sep = ".")
down_nit_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO.20200116./down_nitrated/down_nitrated", sample, 'txt', sep = ".")
up_non_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO.20200116./up_nonnitrated/up_nonnitrated", sample, 'txt', sep = ".")
down_non_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO.20200116./down_nonnitrated/down_nonnitrated", sample, 'txt', sep = ".")

write.table(up_nitrated, up_nit_out, sep='\t', quote = F, row.names = F, col.names = T)
write.table(down_nitrated, down_nit_out, sep='\t', quote = F, row.names = F, col.names = T)
write.table(up_nonnitrated, up_non_out, sep='\t', quote = F, row.names = F, col.names = T)
write.table(down_nonnitrated, down_non_out, sep='\t', quote = F, row.names = F, col.names = T)

## Write table about GO size
up_nit_size <- DEP_up_protein_Nitrated %>% length()
down_nit_size <- DEP_down_protein_Nitrated %>% length()
up_non_size <- DEP_up_protein_NonNitrated %>% length()
down_non_size <- DEP_down_protein_NonNitrated %>% length()

u_nit <- data.frame('Sample_ID' = sample, 'GO_Nitrated_Overexpressed' = up_nit_size)
d_nit <- data.frame('Sample_ID' = sample, 'GO_Nitrated_Underexpressed' = down_nit_size)
u_non <- data.frame('Sample_ID' = sample, 'GO_NonNitrated_Overexpressed' = up_non_size)
d_non <- data.frame('Sample_ID' = sample, 'GO_NonNitrated_Underexpressed' = down_non_size)
  
## Save to file
u_nit_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO_size.20200116./up_nitrated/GO_size_nitrated_overexpressed", sample, 'txt', sep = ".")
d_nit_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO_size.20200116./down_nitrated/GO_size_nitrated_overexpressed", sample, 'txt', sep = ".")
u_non_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO_size.20200116./up_nonnitrated/GO_size_nitrated_overexpressed", sample, 'txt', sep = ".")
d_non_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO_size.20200116./down_nonnitrated/GO_size_nitrated_overexpressed", sample, 'txt', sep = ".")

write.table(u_nit, u_nit_out, sep='\t', quote = F, row.names = F, col.names = T)
write.table(d_nit, d_nit_out, sep='\t', quote = F, row.names = F, col.names = T)
write.table(u_non, u_non_out, sep='\t', quote = F, row.names = F, col.names = T)
write.table(d_non, d_non_out, sep='\t', quote = F, row.names = F, col.names = T)




