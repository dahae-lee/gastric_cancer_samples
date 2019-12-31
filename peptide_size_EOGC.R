# Proteomics data_peptide size
options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)

# load data
# f = '~/GoogleDrive/EOGC_Nitration/summary_files_201912/N13T236_Global_nitration.xlsx'

f = args[1]
sample = strsplit(strsplit(f, 'summary_files_201912/', fixed = T)[[1]][2], '_')[[1]][1]

## Load data
e = as.data.frame(read_excel(f, sheet='Peptide View'))
samples = c('norm1', 'tumor1', 'norm2', 'tumor2')
colnames(e)[5:8] <- samples

## Remove peptides with missing expression
e = e[!(e$norm1==0 | e$norm2==0 | e$tumor1==0 | e$tumor2==0),]

# pull the number of rows as data frame
t <- as.data.frame(nrow(e))
t$name <- sample
names(t) <- c("peptide_size", "sample")
t_out = paste('~/GoogleDrive/gastric_cancer_samples/Tables/peptide_size/
              table.peptide_size', sample, 'txt', sep='.')
write.table(t, t_out, sep='\t', quote = F, row.names = F, col.names = T)
