options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)

f = '~/GoogleDrive/gastric_cancer_samples/Tables/DEP_blastp/table.DEP_blastp.N13T236.txt'
# f = args[1]
sample = strsplit(f, '.', fixed = T)[[1]][3]

## Load data
m = as.data.frame(read.delim(f))

## filter the ones which have overexpressed or underexpressed
library(tidyr)
library(dplyr)
g_name_ov <- m %>% filter(adj.P.Val < 0.05 & t > 0) %>% pull(symbol)
g_name_un <- m %>% filter(adj.P.Val < 0.05 & t < 0) %>% pull(symbol) 

## run gprofiler
library(gProfileR)
gp_over <- gprofiler(g_name_ov, organism = "hsapiens")
gp_under <- gprofiler(g_name_un, organism = "hsapiens")

## Save in file
over_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO/overexpressed_gp/overexpressed", sample, 'txt', sep = ".")
under_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO/underexpressed_gp/underexpressed", sample, 'txt', sep = ".")

write.table(gp_over, over_out, sep='\t', quote = F, row.names = F, col.names = T)
write.table(gp_under,under_out, sep='\t', quote = F, row.names = F, col.names = T)