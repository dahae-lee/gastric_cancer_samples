options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)

# f = '~/GoogleDrive/EOGC_Nitration/summary_files_201912/N5357T5374_Global_nitration.xlsx'
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

#change data frame into the tidy form in order to make plot
library(tidyverse)
n1 <- e[,c(1,5)] %>% mutate(celltype = 'norm1') %>% rename('expression' = 'norm1')
t1 <- e[,c(2,5)] %>% mutate(celltype = 'tumor1') %>% rename('expression' = 'tumor1')
n2 <- e[,c(3,5)] %>% mutate(celltype = 'norm2') %>% rename('expression' = 'norm2')
t2 <- e[,c(4,5)] %>% mutate(celltype = 'tumor2') %>% rename('expression' = 'tumor2')
a <- rbind.data.frame(n1, t1, n2, t2)

#draw plot

pa <- a %>% 
  ggplot(aes(x = celltype, y = expression, fill = celltype))+
  geom_boxplot(position=position_dodge(1))+
  scale_fill_manual(breaks = c("norm1", "norm2", "tumor1", "tumor2"), 
                    values=c("grey", "grey", "firebrick", "firebrick"))+
  xlab("Cell Type")+
  ylab("Peptide epxpression")+
  ggtitle(sample)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pa

#save plot
pa_out <- paste('~/GoogleDrive/gastric_cancer_samples/Figures/peptide_expression/
                plot_peptide_expression_level', sample, 'pdf', sep = ".")
ggsave(pa_out, pa)
