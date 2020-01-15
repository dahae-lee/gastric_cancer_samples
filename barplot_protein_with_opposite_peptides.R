# 데이터의 상위 폴더를 dir에 지정
dir = "~/GoogleDrive/gastric_cancer_samples/Tables/opposite_protein"

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
i <- seq(1,59,by=2)
a <- a[-i,]
a <- a %>% rename('sample_id' = 'V1', 'opposite_peptide_prop' = 'V2')
a$opposite_peptide_prop <- as.numeric(a$opposite_peptide_prop)


# make plot
pa <- a %>%
  ggplot(aes(x = sample_id, y = opposite_peptide_prop, fill = sample_id))+
  geom_bar(stat = "identity")+
  scale_colour_viridis_d(alpha = 1, begin = 0, end = 1,
                         direction = 1, option = "D", aesthetics = "fill")+
  xlab("Sample ID")+
  ylab("Proportion")+
  ggtitle("Proportion of Proteins with Opposite Peptides by Sample")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank())
pa

ggsave("Figures/plot.proportion_of_proteins_with_opposite_peptides.pdf", pa, width = 8, height = 4)

