options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)

# f = '~/GoogleDrive/gastric_cancer_samples/Tables/DEP_blastp/table.DEP_blastp.N39T40.txt'
f = args[1]
sample = strsplit(f, '.', fixed = T)[[1]][3]

## Load data
m = as.data.frame(read.delim(f))

# Sort the proteins that include both DEP_up and DEP_down peptide
library(dplyr)
total_protein <- m %>% pull(symbol) %>% unique()
DEP_up_protein <- m %>% filter(adj.P.Val < 0.05 & t > 0) %>% pull(symbol) %>% unique()
DEP_down_protein <- m %>% filter(adj.P.Val < 0.05 & t < 0) %>% pull(symbol) %>% unique()

opposite_protein <- total_protein[(total_protein %in% DEP_up_protein &
                 total_protein %in% DEP_down_protein) == 'TRUE']
e <- m[m$symbol %in% opposite_protein,]

# Save the protein with opposite peptides proportion
total_protein <- m %>% pull(symbol) %>% unique() %>% length
opposite_protein_num <- e %>% pull(symbol) %>% unique() %>% length
opposite_proportion <- opposite_protein_num / total_protein

y <- data.frame("Sample_ID" = sample, "Opposite_Protein_Proportion" = opposite_proportion)
y_out = paste('~/GoogleDrive/gastric_cancer_samples/Tables/opposite_protein/table.opposite_protein', sample, 'txt', sep='.')
write.table(y, y_out, sep='\t', quote = F, row.names = F, col.names = T)

# make plot between peptide length and t value
library(ggplot2)
e$peptide_length <- e$Sequence2 %>% nchar()
p <- e %>% ggplot(aes(x = peptide_length, y = abs(t)))+
  geom_point()+
  xlab("Peptide Length")+
  ylab("t value(abs)")+
  ggtitle(sample)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save plot
p_out <- paste('~/GoogleDrive/gastric_cancer_samples/Figures/t_by_peptide_length/
                plot_t_by_peptide_length', sample, 'pdf', sep = ".")
ggsave(p_out, p, width = 6, height = 4)

