# 데이터의 상위 폴더를 dir에 지정
dir = "~/GoogleDrive/gastric_cancer_samples/Tables/peptide_expression"

# 폴더에 존재하는 하위 폴더(=샘플명) 리스트 불러오기
# (하위 폴더가 없으면 생략)
sample_raw <- list.files(dir)

# 샘플 파일 주소 리스트 만들기(file.path = paste(..., sep = '/'))
sampleFiles <- file.path(dir, sample_raw)

# 파일 주소에 해당하는 txt파일 읽어들이기
tmp <- lapply(seq(1,30), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

# 매트릭스 하나로 bind
library(tidyverse)
a <- rbind_list(tmp)
a <- a[-1,]
colnames(a) <- c('norm1', 'tumor1', 'norm2',	'tumor2', 'sample_id')

#change data frame into the tidy form in order to make plot
n1 <- a[,c(1,5)] %>% mutate(celltype = 'norm1') %>% rename('expression' = 'norm1')
t1 <- a[,c(2,5)] %>% mutate(celltype = 'tumor1') %>% rename('expression' = 'tumor1')
n2 <- a[,c(3,5)] %>% mutate(celltype = 'norm2') %>% rename('expression' = 'norm2')
t2 <- a[,c(4,5)] %>% mutate(celltype = 'tumor2') %>% rename('expression' = 'tumor2')
a <- rbind.data.frame(n1, t1, n2, t2)

# make plots
pa <- a %>% 
  ggplot(aes(x = sample_id, y = expression, fill = celltype))+
  geom_boxplot(position=position_dodge(1))+
  scale_colour_viridis_d(alpha = 1, begin = 0, end = 1,
                         direction = 1, option = "D", aesthetics = "fill")
pa

ggsave("Figures/plot.peptide_expression_by_samples_191231.pdf", pa, width = 30, height = 10)

