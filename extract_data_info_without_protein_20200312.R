# Proteomics data_peptide size
options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)
date = "20200417"

# load data
f = '~/GoogleDrive/EOGC_Nitration/no_protein/N155T156_Nitration_in_Global_embed.xlsx'

# f = args[1]
sample = strsplit(strsplit(f, 'no_protein/', fixed = T)[[1]][2], '_')[[1]][1]

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

#sort out only nitrated ones
library(tidyverse)
e1 <- e %>% filter(isNitrated == '1')

# pull the number of rows as data frame
t <- data.frame("Date" = "201912", "sample_id" = sample, "peptide_size" = nrow(e), 
                "nitrated_peptide_size" = nrow(e1),
                "protein_size" = "0")

t_out = paste('~/GoogleDrive/gastric_cancer_samples/Tables/data_QC.202000417/201912_data_info', sample, date, 'txt', sep='.')
write.table(t, t_out, sep='\t', quote = F, row.names = F, col.names = F)


