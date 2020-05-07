# Proteomics data_total protein expression level
options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)
date = "20200421"

# load data
# f = '/Users/dahaelee/GoogleDrive/EOGC_Nitration/summary_files_20200309/N111T112_Global_nitration.xlsx'

f = args[1]
sample = strsplit(strsplit(f, '/', fixed = T)[[1]][7], '_')[[1]][1]
Date = strsplit(strsplit(f, 'summary_files_')[[1]][2], '/')[[1]][1]


## Load protein data
d <- as.data.frame(read_excel(f, sheet='Protein View'))
samples = c('norm1', 'tumor1', 'norm2', 'tumor2')
colnames(d)[10:13] <- samples

## Sum the protein expression data
d1 <- d[,10:13]
total_protein_expression <- apply(d1,2,sum)

## Write table
t <- data.frame("date" = Date, "sample_id" = sample)
t1 <- t(as.data.frame(total_protein_expression))
t2 <- cbind.data.frame(t,t1)
rownames(t2) <- NULL

t2_out = paste('~/GoogleDrive/gastric_cancer_samples/Tables/total_protein_expression.20200421/total_protein_expression', Date, sample, date, 'txt', sep='.')
write.table(t2, t2_out, sep='\t', quote = F, row.names = F, col.names = F)
