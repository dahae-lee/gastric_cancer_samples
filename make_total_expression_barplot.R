# 데이터의 상위 폴더를 dir에 지정
dir = "~/GoogleDrive/gastric_cancer_samples/Tables/total_protein_expression.20200312"
date = '20200313'

# 폴더에 존재하는 하위 폴더(=샘플명) 리스트 불러오기
# (하위 폴더가 없으면 생략)
sample_raw <- list.files(dir)

# 샘플 파일 주소 리스트 만들기(file.path = paste(..., sep = '/'))
sampleFiles <- file.path(dir, sample_raw)

# 파일 주소에 해당하는 txt파일 읽어들이기
tmp <- lapply(seq(1,63), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

# 매트릭스 하나로 bind
library(tidyverse)
a <- rbind_list(tmp)
a$batch <- ifelse(a$V1 == 201912, '1', '2')
a <- a %>% rename('Date' = 'V1', 'Sample_ID' = 'V2',
                  'norm1' = 'V3', 'tumor1' = 'V4', 
                  'norm2' = 'V5', 'tumor2' = 'V6')

# Save to file
a_out = paste('~/GoogleDrive/gastric_cancer_samples/Tables/total_protein_expression', date, 'txt', sep='.')
write.table(a, a_out, sep='\t', row.names = F, col.names = T)

# Change data frame into the tidy form in order to make plot
n1 <- a[,c(1,2,3,7)] %>% mutate(celltype = 'norm1') %>% 
  rename('expression' = 'norm1')
n2 <- a[,c(1,2,5,7)] %>% mutate(celltype = 'norm2') %>% 
  rename('expression' = 'norm2')
t1 <- a[,c(1,2,4,7)] %>% mutate(celltype = 'tumor1') %>% 
  rename('expression' = 'tumor1')
t2 <- a[,c(1,2,6,7)] %>% mutate(celltype = 'tumor2') %>% 
  rename('expression' = 'tumor2')
a = rbind.data.frame(n1, n2, t1, t2)

# Remove duplicated samples
c <- as.data.frame(table(a$Sample_ID))
dup <- c %>% filter(Freq != 4) %>% pull(Var1)
a1 <- a %>% filter(!(Sample_ID %in% dup & Date == 202003))

# make plots
pa <- a1 %>%
  ggplot(aes(x = reorder(Sample_ID,Date), y = log(expression), fill = batch))+
  geom_bar(stat = "identity")+
  scale_colour_viridis_d(alpha = 1, begin = 0, end = 1,
                         direction = 1, option = "D", aesthetics = "fill")+
  facet_grid(celltype~.)+
  xlab("Sample ID")+
  ylab("Total Protein Expression Size")+
  ggtitle("Total Protein Expression Size by Samples")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank())
pa

ggsave("Figures/plot.total_protein_size_by_samples_203012.pdf", pa, width = 8, height = 4)
