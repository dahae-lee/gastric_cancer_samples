options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)

f = '~/GoogleDrive/EOGC_Nitration/summary_files_201912/N5357T5374_Global_nitration.xlsx'
# f = args[1]
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

## ------------------- vv
## DEP
library(edgeR)
d = e[,c(samples, 'pid')]
rownames(d) <- d$pid

group <- gsub('1|2', '', samples)
plotMDS(d, col = as.numeric(group))
mm <- model.matrix(~ 0 + group)

f_voom = paste('Figures/plot.voom_variance', sample, 'pdf', sep='.')
pdf(f_voom, width = 5, height = 5)
y <- voom(d, mm, plot = T)
dev.off()

fit <- lmFit(y, mm)
contr <- makeContrasts(groupnorm - grouptumor, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
res <- topTable(tmp, sort.by = "P", n = Inf)
table(res$adj.P.Val < 0.05)
hist(res$adj.P.Val, breaks = 1000)
res$pid <- rownames(res)

res = res[res$adj.P.Val <= 0.05,]
e = e[e$pid %in% res$pid,]
m = merge(e, res, by='pid')

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
f_blastp = paste('Tables/blastp/table.blastp', sample, 'txt.gz', sep='.')

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
res = read.table(f_blastp, header = F, col.names = c('pid', 'sseqid', 'pident', 
                                                     'length', 'mismatch', 'gapopen', 
                                                     'qstart', 'qend', 'sstart', 'send', 
                                                     'evalue', 'bitscore'))
res = res[res$pident >= 100,]
m = merge(m, res[,c('pid', 'sseqid', 'pident', 'sstart', 'send')], by='pid', all.x=T)
m$uniprot_ids <- as.character(do.call(rbind.data.frame, strsplit(m$sseqid, '|', fixed = T))[[2]])

## Add gene information
hgnc = read.delim('~/Dropbox/Resources/hgnc/hgnc_complete_set.txt')
m = merge(m, hgnc[,c('uniprot_ids', 'ensembl_gene_id', 'symbol', 'hgnc_id')], by='uniprot_ids')


## Save to file
f_out = paste('Tables/table.DEP_blastp', sample, 'txt', sep='.')
write.table(m, f_out, sep='\t', quote = F, row.names = F, col.names = T)


