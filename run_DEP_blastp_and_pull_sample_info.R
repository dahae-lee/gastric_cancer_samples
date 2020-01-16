options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)
library(tidyverse)

date = "20200116"

option_test = F
if (option_test){
  f = '~/GoogleDrive/EOGC_Nitration/summary_files_201912/N13T236_Global_nitration.xlsx' 
} else {
  f = args[1]  
}

sample = strsplit(strsplit(f, 'summary_files_201912/', fixed = T)[[1]][2], '_')[[1]][1]

## Load data
e = as.data.frame(read_excel(f, sheet='Peptide View'))
samples = c('norm1', 'tumor1', 'norm2', 'tumor2')
colnames(e)[5:8] <- samples

## Parse nitrated peptides
e$isY45 <- ifelse(grepl('Y[45]', e$Sequence, fixed = T), 1, 0)
e$isC29 <- ifelse(grepl('C[29]', e$Sequence, fixed = T), 1, 0)
#둘 중 하나라도 1이면 1이 되는 column을 만든다.
e$isNitrated <- ifelse(e$isC29==1 | e$isY45 == 1, 1, 0)

#e$Sequence에서 알파벳만 빼서 e$Sequence2에 넣어 준다.
e$Sequence2 <- gsub("[^a-zA-Z]", "", e$Sequence)
e$norm1 <- log2(e$norm1 + 1)
e$norm2 <- log2(e$norm2 + 1)
e$tumor1 <- log2(e$tumor1 + 1)
e$tumor2 <- log2(e$tumor2 + 1)
e$pid <- paste('p', seq(1,nrow(e)), sep='')
#paste함수는 문자열을 합쳐 주는데, 이 때 sep=''은 그 문자열 사이에 무언가를 넣어주라는 
# 의미이다. 여기에서는 예를 들어 `paste(1,2,3,4,sep='-')`는 `1-2-3-4`로 출력된다. 여기에서는 
# 아무것도 넣지 않았는데 이는 그냥 합친다는 뜻이 아니라 공백없이 붙이라는 의미가 된다. 
# sep=''를 넣지 않은 경우 p와 seq(1,nrow(e)) 사이에 띄어쓰기가 된 것처럼 공백이 생긴다.

## Remove peptides with missing expression
e = e[!(e$norm1==0 | e$norm2==0 | e$tumor1==0 | e$tumor2==0),]

## ------------------- vv
## DEP
# if (!requireNamespace("BiocManager", quietly = TRUE)) + install.packages("BiocManager")
# BiocManager::install("edgeR")
library(edgeR)


# d = e[nchar(e$Sequence2)/e$`Molecular Weight` >= 0.008, c(samples, 'pid')]
# d = e[nchar(e$Sequence2) >= 20,c(samples, 'pid')]

d = e[,c(samples, 'pid')]
#d의 rowname에 pid column을 넣어 주고 pid column은 지운다.
rownames(d) <- d$pid
d$pid <- NULL

#아까 지정했던 samples object에서 1또는2를 제거해준다.
group <- gsub('1|2', '', samples)
batch <- gsub('norm|tumor', '', samples)

# Remove batch effect
library(sva)
d = limma::removeBatchEffect(d, batch = batch)

svafit <- sva(d, model.matrix( ~ group ), B = 10)

#개체들 사이의 유사성/비유사성을 측정하여 2차원 또는 3차원 공간상에 점으로  분석 방법. 개체들간의 근접성(proximity)을 시각화하여 데이터 속에 잠재해 있는 패턴이나 구조를 찾아내표현하는는 통계 기법. 개체들간의 거리 계산은 유클리드 거리 행렬을 사용한다. 상대적 거리의 정확도를 높이기 위해 적합한 정도를 스트레스 값(stress value)으로 나타낸다.
plotMDS(d, col = as.numeric(group))

if (svafit$n.sv != 0){
  # SVA is available
  mm <- model.matrix(~ group + svafit$sv)  
  
  # Re-scaled expression matrix
  d2 = limma::removeBatchEffect(d, batch = svafit$sv, design = model.matrix(~ group))
  d2 = as.data.frame(d2) %>% mutate(pid = rownames(d2))
} else {
  # SVA is not available
  mm <- model.matrix(~ group)  
  d2 = d
  d2 = as.data.frame(d2) %>% mutate(pid = rownames(d2))
}

f_voom = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Figures/voom_variance', date, '/plot.voom_variance', sample, date, 'pdf', sep='.')
pdf(f_voom, width = 5, height = 5)
#voom은 RNA sequence data를 linear modeling으로 바꾸어준다.
y <- voom(d, mm, plot = T)
dev.off()

fit <- lmFit(y, mm)
fit <- eBayes(fit)

qvals <- p.adjust(fit$p.value[,'grouptumor'], method = 'BH')
res <- tibble(t = fit$coefficients[,'grouptumor'], 
              pval = fit$p.value[,'grouptumor'],
              adj.P.Val = qvals,
              pid = rownames(d),
              is_significant = adj.P.Val <= 0.05) 

table(res$adj.P.Val <= 0.05)

d3 = merge(res, d2, by='pid')

m = merge(e %>% select(-one_of(samples)), 
          d3, by='pid')

## ------------------- vv
## blastp

## Create an input fasta
writeLines(text = paste('>', m$pid, '\n', m$Sequence2, sep=''), con = 'test.fa')

## Run blastp
blastp = "blastp"
blast_db = "/Users/dahaelee/GoogleDrive/resources/uniprot/Uniprot_homo_sapiens_20170724_reviewed"
input = "test.fa"
evalue = 1e-6
format = 6
nThreads = 2
f_blastp = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/blastp', date, '/table.blastp', sample, date, 'txt.gz', sep='.')

if (!file.exists(f_blastp)){
  system2(command = blastp, 
          args = c("-db", blast_db, 
                   "-query", input, 
                   "-outfmt", format, 
                   "-evalue", evalue, 
                   "-num_threads", nThreads,
                   "| gzip >", f_blastp))
}


## Merge to a super table
res = read.table(gzfile(f_blastp), header = F, col.names = c('pid', 'sseqid', 'pident', 
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
f_out = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/DEP_blastp', date, '/table.DEP_blastp', sample, date, 'txt', sep='.')
write.table(m, f_out, sep='\t', quote = F, row.names = F, col.names = T)

# Abstract the information of protein size, DEP up&down protein size
library(dplyr)
total_protein <- m %>% pull(symbol) %>% unique() %>% length
DEP_up_peptide <- m %>% filter(adj.P.Val < 0.05 & t > 0) %>% count(Sequence) %>% nrow()
DEP_down_peptide <- m %>% filter(adj.P.Val < 0.05 & t < 0) %>% count(Sequence) %>% nrow()
DEP_up_protein_Nitrated <- m %>% filter(adj.P.Val < 0.05 & t > 0 & isNitrated == '1') %>% pull(symbol) %>% unique() %>% length
DEP_down_protein_Nitrated <- m %>% filter(adj.P.Val < 0.05 & t < 0 & isNitrated == '1') %>% pull(symbol) %>% unique() %>% length
DEP_up_protein_NonNitrated <- m %>% filter(adj.P.Val < 0.05 & t > 0 & isNitrated == '0') %>% pull(symbol) %>% unique() %>% length
DEP_down_protein_NonNitrated <- m %>% filter(adj.P.Val < 0.05 & t < 0 & isNitrated == '0') %>% pull(symbol) %>% unique() %>% length
NonDEP_nitrated_Protein <- m %>% filter(adj.P.Val >= 0.05 & isNitrated == 1)  %>% pull(symbol) %>% unique() %>% length
NonDEP_NonNitrated_Protein <- m %>% filter(adj.P.Val >= 0.05 & isNitrated == 0)  %>% pull(symbol) %>% unique() %>% length

# Make a table that includes the size information

y <- data.frame("sample_id" = sample, "Total_Protein" = total_protein, 
                "DEP_Peptide_Up" = DEP_up_peptide, "DEP_Peptide_Down" = DEP_down_peptide, 
                "DEP_up_Nitrated_Protein" = DEP_up_protein_Nitrated, 
                "DEP_down_Nitrated_Protein" = DEP_down_protein_Nitrated,
                "DEP_up_NonNitrated_Protein" = DEP_up_protein_NonNitrated, 
                "DEP_down_NonNitrated_Protein" = DEP_down_protein_NonNitrated,
                "NonDEP_Nitrated_Protein" = NonDEP_nitrated_Protein, 
                "NonDEP_NonNitrated_Protein" = NonDEP_NonNitrated_Protein)

## Save to file
y_out = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/sample_info', date, '/table.sample_info', sample, date, 'txt', sep='.')
write.table(y, y_out, sep='\t', quote = F, row.names = F, col.names = T)

