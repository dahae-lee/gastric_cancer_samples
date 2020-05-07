library(ggplot2)
library(tidyverse)
date = '20200422'

## Load data
d <- read.delim('~/GoogleDrive/gastric_cancer_samples/Tables/table.DEG_and_GO_size_by_sample.20200421.txt')

## Draw DEG plot
p <- mosaicplot(~)

##Save
p_out <- paste('~/GoogleDrive/gastric_cancer_samples/Figures/Figure.4C', date,'pdf', sep = ".")
ggsave(p_out, p, width = 6, height = 4)