options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)

# f = '~/GoogleDrive/EOGC_Nitration/summary_files_201912/N39T40_Global_nitration.xlsx'
f = args[1]
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
d = e[,c(samples, 'pid')]
#d의 rowname에 pid column을 넣어 주고 pid column은 지운다.
rownames(d) <- d$pid
d$pid <- NULL

#아까 지정했던 samples object에서 1또는2를 제거해준다.
group <- gsub('1|2', '', samples)
#개체들 사이의 유사성/비유사성을 측정하여 2차원 또는 3차원 공간상에 점으로  분석 방법. 개체들간의 근접성(proximity)을 시각화하여 데이터 속에 잠재해 있는 패턴이나 구조를 찾아내표현하는는 통계 기법. 개체들간의 거리 계산은 유클리드 거리 행렬을 사용한다. 상대적 거리의 정확도를 높이기 위해 적합한 정도를 스트레스 값(stress value)으로 나타낸다.
plotMDS(d, col = as.numeric(group))
mm <- model.matrix(~ 0 + group)

f_voom = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Figures/voom_variance/plot.voom_variance', sample, 'pdf', sep='.')
pdf(f_voom, width = 5, height = 5)
#voom은 RNA sequence data를 linear modeling으로 바꾸어준다.
y <- voom(d, mm, plot = T)
dev.off()

fit <- lmFit(y, mm)
#makeContrasts 함수는 매개 변수 집합(a set of parameter)간의 대비(contrast)를 숫자 행렬로 표현한다.
contr <- makeContrasts(groupnorm - grouptumor, levels = colnames(coef(fit)))
#contrasts.fit는 microarray data에 적합한 linear 모델이 주어지면 주어진 set of contrasts에 대한 추정 계수 및 표준 오류를 계산한다.
tmp <- contrasts.fit(fit, contr)
#eBayes 함수는 microarray linear model fit이 주어지면 표준 오차를 공통 값으로 경험적으로 Bayes 중재함으로써 중간 t- 통계량, F- 통계량, 미분 표현의 로그 차수를 계산한다.
tmp <- eBayes(tmp)
#topTable 함수는 a linear model fit에서 a table of the top-ranked genes를 extract한다. n = Inf라는 건 n이 무한대까지, 즉 모든 data를 고른다는 의미인 것 같다.
res <- topTable(tmp, sort.by = "P", n = Inf)
#p value가 0.05 미만인 것들과 그렇지 않은 것들의 숫자를 고른다.
table(res$adj.P.Val < 0.05)
hist(res$adj.P.Val, breaks = 1000)
#res에서 pid라는 column을 만들어서 res의 rowname에 있는 정보들을 하나의 column으로 빼 준다.
res$pid <- rownames(res)

#res에서 p value가 0.05 미만인 것들만 추린다. 
# 다만, 여기에서는 전체 protein과 pepetide 갯수를 세어 주기 위해 이 과정은 생략한다.
# res = res[res$adj.P.Val <= 0.05,]
#e의 pid column이 res에 포함되어 있는 것들, 즉 p value가 0.05 이하인 항목들만 골라낸다.
e = e[e$pid %in% res$pid,]
#e랑 m을 옆으로(가로로) pid에 따라 합쳐버린다.
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
f_blastp = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/blastp/table.blastp', sample, 'txt.gz', sep='.')

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
f_out = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/DEP_blastp/table.DEP_blastp', sample, 'txt', sep='.')
write.table(m, f_out, sep='\t', quote = F, row.names = F, col.names = T)

# Abstract the information of protein size, DEP up&down protein size
library(dplyr)
total_protein <- m %>% pull(symbol) %>% unique() %>% length
DEP_up_peptide <- m %>% filter(adj.P.Val < 0.05 & t > 0) %>% count(Sequence) %>% nrow()
DEP_down_peptide <- m %>% filter(adj.P.Val < 0.05 & t < 0) %>% count(Sequence) %>% nrow()
DEP_up_protein <- m %>% filter(adj.P.Val < 0.05 & t > 0) %>% pull(symbol) %>% unique() %>% length
DEP_down_protein <- m %>% filter(adj.P.Val < 0.05 & t < 0) %>% pull(symbol) %>% unique() %>% length
NonDEP_nitrated_Protein <- m %>% filter(adj.P.Val >= 0.05 & isNitrated == 1)  %>% pull(symbol) %>% unique() %>% length
NonDEP_NonNitrated_Protein <- m %>% filter(adj.P.Val >= 0.05 & isNitrated == 0)  %>% pull(symbol) %>% unique() %>% length

# Make a table that includes the size information

y <- data.frame("sample_id" = sample, "Total_Protein" = total_protein, 
                "DEP_Peptide_Up" = DEP_up_peptide, "DEP_Peptide_Down" = DEP_down_peptide, 
                "DEP_Protein_Up" = DEP_up_protein, "DEP_Protein_Down" = DEP_down_protein,
                "NonDEP_Nitrated_Protein" = NonDEP_nitrated_Protein, 
                "NonDEP_NonNitrated_Protein" = NonDEP_NonNitrated_Protein)

## Save to file
y_out = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/sample_info/table.sample_info', sample, 'txt', sep='.')
write.table(y, y_out, sep='\t', quote = F, row.names = F, col.names = T)

