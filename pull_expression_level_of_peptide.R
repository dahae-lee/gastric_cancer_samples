options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)

#f = '~/GoogleDrive/EOGC_Nitration/summary_files/N155T156_Nitration_in_Global.xlsx'
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

e$Sequence2 <- gsub("[^a-zA-Z]", "", e$Sequence)
e$norm1 <- log2(e$norm1 + 1)
e$norm2 <- log2(e$norm2 + 1)
e$tumor1 <- log2(e$tumor1 + 1)
e$tumor2 <- log2(e$tumor2 + 1)
e$pid <- paste('p', seq(1,nrow(e)), sep='')

## Remove peptides with missing expression
e = e[!(e$norm1==0 | e$norm2==0 | e$tumor1==0 | e$tumor2==0),]

# select the expression levels only and add the sample number in new column
e <- e[,5:8]
e$sample_id <- sample
e_out = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/peptide_expression/table.peptide_expression', sample, 'txt', sep='.')
write.table(e, e_out, sep='\t', quote = F, row.names = F, col.names = T)
