options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)

f = '~/GoogleDrive/gastric_cancer_samples/Tables/DEP_blastp/table.DEP_blastp.N39T40.txt'
# f = args[1]
sample = strsplit(strsplit(f, 'table.DEP_blastp.', fixed = T)[1][[2]], '\\.')[[1]][2]

## Load data
m = as.data.frame(read.delim(f))

# Sort the proteins that include both DEP_up and DEP_down peptide
library(dplyr)
library(rlist)
total_protein <- m %>% pull(symbol) %>% unique()
DEP_up_protein <- m %>% filter(adj.P.Val < 0.05 & t > 0) %>% pull(symbol) %>% unique()
DEP_down_protein <- m %>% filter(adj.P.Val < 0.05 & t < 0) %>% pull(symbol) %>% unique()
total_protein %>% list.filter((total_protein %in% DEP_up_protein &
                                 total_protein %in% DEP_down_protein) == 'TRUE')

