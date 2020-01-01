options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)

## Merge to a super table
res = read.table(f_blastp, header = F, col.names = c('pid', 'sseqid', 'pident', 
                                                     'length', 'mismatch', 'gapopen', 
                                                     'qstart', 'qend', 'sstart', 'send', 
                                                     'evalue', 'bitscore'))
res = res[res$pident >= 100,]
m = merge(m, res[,c('pid', 'sseqid', 'pident', 'sstart', 'send')], by='pid', all.x=T)
m$sseqid <- as.character(m$sseqid)
m$uniprot_ids <- as.character(do.call(rbind.data.frame, strsplit(m$sseqid, '|', fixed = T))[[2]])

## Add gene information
hgnc = read.delim('~/GoogleDrive/resources/hgnc/hgnc_complete_set.txt')
m = merge(m, hgnc[,c('uniprot_ids', 'ensembl_gene_id', 'symbol', 'hgnc_id')], by='uniprot_ids')

# Save File
# f_out = paste('Tables/DEP_blastp/table.DEP_blastp', sample, 'txt', sep='.')
# write.table(m, f_out, sep='\t', quote = F, row.names = F, col.names = T)

# Abstract the information of protein size, DEP up&down protein size
library(dplyr)
total_protein <- m %>% count(symbol) %>% nrow()
DEP_up_peptide <- m %>% filter(adj.P.Val < 0.05 & t > 0) %>% count(Sequence) %>% nrow()
DEP_down_peptide <- m %>% filter(adj.P.Val < 0.05 & t < 0) %>% count(Sequence) %>% nrow()
DEP_up_protein <- m %>% filter(adj.P.Val < 0.05 & t > 0) %>% count(symbol) %>% nrow()
DEP_down_protein <- m %>% filter(adj.P.Val < 0.05 & t < 0) %>% count(symbol) %>% nrow()

# Make a table that includes the size information

y <- data.frame("sample_id" = sample, "Total_Protein" = total_protein, 
                "DEP_Peptide_Up" = DEP_up_peptide, "DEP_Peptide_Down" = DEP_down_peptide, 
                "DEP_Protein_Up" = DEP_up_protein, "DEP_Protein_Down" = DEP_down_protein)

## Save to file
y_out = paste('Tables/sample_info/table.sample_info', sample, 'txt', sep='.')
write.table(y, y_out, sep='\t', quote = F, row.names = F, col.names = T)

