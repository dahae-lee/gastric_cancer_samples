# Proteomics data_peptide size
options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)

# f = '~/GoogleDrive/EOGC_Nitration/summary_files_201912/N13T236_Global_nitration.xlsx'

# Split the sample ID
f = args[1]
sample = strsplit(strsplit(f, 'summary_files_201912/', fixed = T)[[1]][2], '_')[[1]][1]

## Load data
e = as.data.frame(read_excel(f, sheet='Peptide View'))
samples = c('norm1', 'tumor1', 'norm2', 'tumor2')
colnames(e)[5:8] <- samples

## Parse nitrated peptides
e$isY45 <- ifelse(grepl('Y[45]', e$Sequence, fixed = T), 1, 0)
e$isC29 <- ifelse(grepl('C[29]', e$Sequence, fixed = T), 1, 0)
e$isNitrated <- ifelse(e$isC29==1 | e$isY45 == 1, 1, 0)

## Remove peptides with missing expression
e = e[!(e$norm1==0 | e$norm2==0 | e$tumor1==0 | e$tumor2==0),]

#sort out only nitrated ones
library(tidyverse)
e1 <- e %>% filter(isNitrated == '1')

# pull the number of rows as data frame
t <- as.data.frame(nrow(e1))
t$name <- sample
names(t) <- c("peptide_size", "sample")
t_out = paste('~/GoogleDrive/gastric_cancer_samples/Tables/nitrated_peptide_size/
              table.nitrated_peptide_size', sample, 'txt', sep='.')
write.table(t, t_out, sep='\t', quote = F, row.names = F, col.names = T)

