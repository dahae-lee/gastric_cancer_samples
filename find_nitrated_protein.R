# Proteomics data_peptide size
options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)
library(tidyverse)
library(dplyr)
date = "20200603"

# load data
# f = '~/GoogleDrive/EOGC_Nitration/summary_files_20200309/N123T124_Global_nitration.xlsx'

f = args[1]
sample = strsplit(strsplit(f, 'summary_files_20200309/', fixed = T)[[1]][2], '_')[[1]][1]

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

#e$Sequence에서 알파벳만 빼서 e$Sequence2에 넣어 준다.
e$Sequence2 <- gsub("[^a-zA-Z]", "", e$Sequence)
e$norm1 <- log2(e$norm1 + 1)
e$norm2 <- log2(e$norm2 + 1)
e$tumor1 <- log2(e$tumor1 + 1)
e$tumor2 <- log2(e$tumor2 + 1)
e$pid <- paste('p', seq(1,nrow(e)), sep='')

#sort out only nitrated ones
library(tidyverse)
# e1 <- e %>% filter(isNitrated == '1')

## Create an input fasta
writeLines(text = paste('>', e$pid, '\n', e$Sequence2, sep=''), con = 'test.fa')

## Run blastp
blastp = "blastp"
blast_db = "/Users/dahaelee/GoogleDrive/resources/uniprot/Uniprot_homo_sapiens_20170724_reviewed"
input = "test.fa"
evalue = 1e-6
format = 6
nThreads = 2
f_blastp = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/blastp.20200603/table.blastp', sample, date, 'txt.gz', sep='.')

if (!file.exists(f_blastp)){
  system2(command = blastp, 
          args = c("-db", blast_db, 
                   "-query", input, 
                   "-outfmt", format, 
                   "-evalue", evalue, 
                   "-num_threads", nThreads,
                   "| gzip >", f_blastp))
}


## Sort the proteins with only 100 pident
res = read.table(gzfile(f_blastp), header = F, col.names = c('pid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'))
res = res[res$pident >= 100,]
res = na.omit(res)
m = merge(e[,c('pid', 'isNitrated')], res[,c('pid', 'sseqid')], by='pid')
m$sseqid <- as.character(res$sseqid)
m$uniprot_ids <- as.character(do.call(rbind.data.frame, strsplit(res$sseqid, '|', fixed = T))[[2]])

## Make the nitrated protein table
a <- m %>% group_by(uniprot_ids) %>% summarise('Number_all_peptide' = n())
b <- m %>% filter(isNitrated == 1) %>% group_by(uniprot_ids) %>% summarise('Number_nitrated_peptide' = n())
d <- merge(a, b, by='uniprot_ids', all.x = T)
d[is.na(d)] <- 0
d$Proportion_of_Nitrated_peptide <- d$Number_nitrated_peptide/d$Number_all_peptide
d$Sample_ID <- sample
d <- d[,c(5,1,2,3,4)]

d_out = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/Protein_with_nitrated_peptide_by_sample.20200603/table.Nitrated_Peptide', sample, date, 'txt', sep='.')
write.table(d, d_out, sep='\t', quote = F, row.names = F, col.names = F)
