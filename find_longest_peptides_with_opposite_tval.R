options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)

# f = '~/GoogleDrive/gastric_cancer_samples/Tables/DEP_blastp/table.DEP_blastp.N39T40.txt'
f = args[1]
sample = strsplit(f, '.', fixed = T)[[1]][3]

## Load data
m = as.data.frame(read.delim(f))

## filter the peptides with longest length
library(tidyverse)
m$peptide_length <- m$Sequence2 %>% nchar()
t <- m %>% group_by(symbol) %>% top_n(1, peptide_length) %>% select(symbol,peptide_length,t)

proteins <- t %>% pull(symbol) %>% unique()

## use the function to know whether there are some longest peptides with opposite t value
l <- lapply(proteins, function(x){
  t1 <- t %>% filter(symbol == x)
  return(c(x, all(t1$t < 0) | all(t1$t > 0))) 
})

l <- as.data.frame(l)
l<- t(l)
rownames(l) <- NULL
colnames(l) <- c("protein", "opposite_t")
l <- as.data.frame(l)
# table(l$opposite_t)

## Save to file
f <- l %>% filter(opposite_t == FALSE)

m1 <- m %>% filter(symbol %in% f$protein)
m1_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/longest_and_opposite_t/protein_longest_and_opposite_t", sample, 'txt', sep = ".")
write.table(m1, m1_out, sep='\t', quote = F, row.names = F, col.names = T)
