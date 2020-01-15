options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)

# f = '~/GoogleDrive/gastric_cancer_samples/Tables/GO/overexpressed_gp/overexpressed.N13T236.txt'
f = args[1]
sample = strsplit(f, '.', fixed = T)[[1]][2]

## Load data
m = as.data.frame(read.delim(f))

## Write table about GO size
GO_size <- nrow(m)

# y <- data.frame('Sample_ID' = sample, 'GO_overexpressed' = GO_size)
y <- data.frame('Sample_ID' = sample, 'GO_underexpressed' = GO_size)

## Save to file
# y_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO_size/overexpressed/GO_size_overexpressed", sample, 'txt', sep = ".")
y_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO_size/underexpressed/GO_size_underexpressed", sample, 'txt', sep = ".")
write.table(y, y_out, sep='\t', quote = F, row.names = F, col.names = T)
