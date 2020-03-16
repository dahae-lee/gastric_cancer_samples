# Proteomics data_peptide size
options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)
date = "20200311"

# load data
# f = '~/GoogleDrive/EOGC_Nitration/summary_files_20200309/N111T112_Global_nitration.xlsx'

f = args[1]
sample = strsplit(strsplit(f, 'summary_files_201912/', fixed = T)[[1]][2], '_')[[1]][1]

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

## Load protein data
d <- as.data.frame(read_excel(f, sheet='Protein View'))
colnames(d)[10:13] <- samples

## Remove proteins with missing expression
d = d[!(d$norm1==0 | d$norm2==0 | d$tumor1==0 | d$tumor2==0),]

# pull the number of rows as data frame
t <- data.frame("Date" = "201912", "sample_id" = sample, "peptide_size" = nrow(e), 
                "nitrated_peptide_size" = nrow(e1),
                "protein_size" = nrow(d))

t_out = paste('~/GoogleDrive/gastric_cancer_samples/Tables/data_QC.20200312/201912_data_info', sample, date, 'txt', sep='.')
write.table(t, t_out, sep='\t', quote = F, row.names = F, col.names = F)
