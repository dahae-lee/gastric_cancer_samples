options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
library(readxl)

f = '~/GoogleDrive/gastric_cancer_samples/Tables/peptide_expression/table.peptide_expression.N5357T5374.txt'
# f = args[1]
sample = strsplit(strsplit(f, 'expression.', fixed = T)[1][[2]], '.')[[1]][2]

