library(tidyverse)
library(ggplot2)
date = '20200423'

## Load data
dir = "~/GoogleDrive/gastric_cancer_samples/Tables/DESeq2_GO.200421/up"
sample_raw <- list.files(dir)
sampleFiles <- file.path(dir, sample_raw)
tmp <- lapply(seq(1,length(sampleFiles)), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

a <- rbind_list(tmp)
colnames(a) <- as.list(a[1,])
a <- a[-1,]
a1 <- a %>% select(term.name) %>% add_column('GO' = 'Up')

dir_d = "~/GoogleDrive/gastric_cancer_samples/Tables/DESeq2_GO.200421/down"
sample_raw_d <- list.files(dir_d)
sampleFiles_d <- file.path(dir_d, sample_raw_d)
tmp_d <- lapply(seq(1,length(sampleFiles_d)), function(i){
  X <- read.delim(sampleFiles_d[i], header = F, sep = '\t')
  X})

b <- rbind_list(tmp_d)
colnames(b) <- as.list(b[1,])
b <- b[-1,]
b1 <- b %>% select(term.name) %>% add_column('GO' = 'Down')

## Merge data
c <- rbind.data.frame(a1,b1)

a2 <- a1 %>% group_by(term.name) %>% summarise('freq' = n()) %>% arrange(desc(freq))
a2 <- a2[2:21,]
b2 <- b1 %>% group_by(term.name) %>% summarise('freq' = n()) %>% arrange(desc(freq))
b2 <- b2[2:21,]

## Draw plot
pa <- a2 %>% ggplot(aes(x=reorder(str_wrap(term.name,45),freq), y=freq))+
  geom_bar(stat = 'identity', fill = '#CD2836')+
  coord_flip()+
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size=7, face = "bold.italic"))+
  xlab("")+
  ylab("")
pa

pb <- b2 %>% ggplot(aes(x=reorder(str_wrap(term.name,45),freq), y=freq))+
  geom_bar(stat = 'identity', fill = '#8985E6')+
  coord_flip()+
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size=7, face = "bold.italic"))+
  xlab("")+
  ylab("")
pb

pa_out <- paste('~/GoogleDrive/gastric_cancer_samples/Figures/Fig.2C/up', date,'pdf', sep = ".")
ggsave(pa_out, pa, width = 6, height = 6)

pb_out <- paste('~/GoogleDrive/gastric_cancer_samples/Figures/Fig.2C/down', date,'pdf', sep = ".")
ggsave(pb_out, pb, width = 6, height = 6)

