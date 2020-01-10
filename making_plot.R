# 데이터의 상위 폴더를 dir에 지정
dir = "~/GoogleDrive/gastric_cancer_samples/Tables/peptide_size"
dir_n = "~/GoogleDrive/gastric_cancer_samples/Tables/nitrated_peptide_size"

# 폴더에 존재하는 하위 폴더(=샘플명) 리스트 불러오기
# (하위 폴더가 없으면 생략)
sample_raw <- list.files(dir)
sample_raw_n <- list.files(dir_n)

# 샘플 파일 주소 리스트 만들기(file.path = paste(..., sep = '/'))
sampleFiles <- file.path(dir, sample_raw)
sampleFiles_n <- file.path(dir_n, sample_raw_n)

# 파일 주소에 해당하는 txt파일 읽어들이기
tmp <- lapply(seq(1,30), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

tmp_n <- lapply(seq(1,30), function(i){
  X <- read.delim(sampleFiles_n[i], header = F, sep = '\t')
  X})

# 매트릭스 하나로 bind
library(tidyverse)
a <- rbind_list(tmp)
i <- seq(1,59,by=2)
a <- a[-i,]
a <- a %>% rename('size' = 'V1', 'sample_id' = 'V2')
a <- a %>% as.numeric(as.factor(size))
a$size <- as.numeric(a$size)

b <- rbind_list(tmp_n)
i <- seq(1,59,by=2)
b <- b[-i,]
b <- b %>% rename('size' = 'V1', 'sample_id' = 'V2')
b$size <- as.numeric(b$size)

# Save to file
write.table(a, 'Tables/Peptide_size_by_samples', col.names = T)
write.table(b, 'Tables/Nitrated_Peptide_size_by_samples', col.names = T)

# make plots
pa <- a %>%
  ggplot(aes(x = sample_id, y = size, fill = sample_id))+
  geom_bar(stat = "identity")+
  scale_colour_viridis_d(alpha = 1, begin = 0, end = 1,
                         direction = 1, option = "D", aesthetics = "fill")+
  xlab("Sample ID")+
  ylab("Peptide Size")+
  ggtitle("Peptide Size by Samples")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
  axis.text.x = element_blank())
pa

pb <- b %>% 
  ggplot(aes(x = sample_id, y = size, fill = sample_id))+
  geom_bar(stat = "identity")+
  scale_colour_viridis_d(alpha = 1, begin = 0, end = 1,
                         direction = 1, option = "D", aesthetics = "fill")+
  xlab("Sample ID")+
  ylab("Nitrated Peptide Size")+
  ggtitle("Nitrated Peptide Size by Samples")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank())
pb

ggsave("Figures/plot.peptide_size_by_samples_191230.pdf", pa, width = 8, height = 4)
ggsave("Figures/plot.nitrated_peptide_size_by_samples_191230.pdf", pb, width = 8, height = 4)
