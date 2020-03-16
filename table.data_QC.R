# 데이터의 상위 폴더를 dir에 지정
dir = "~/GoogleDrive/gastric_cancer_samples/Tables/data_QC.20200312/"
date = '20200313'

# 폴더에 존재하는 하위 폴더(=샘플명) 리스트 불러오기
# (하위 폴더가 없으면 생략)
sample_raw <- list.files(dir)

# 샘플 파일 주소 리스트 만들기(file.path = paste(..., sep = '/'))
sampleFiles <- file.path(dir, sample_raw)

# 파일 주소에 해당하는 txt파일 읽어들이기
tmp <- lapply(seq(1,94), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

# 매트릭스 하나로 bind
library(tidyverse)
a <- rbind_list(tmp)
a$batch <- ifelse(a$V1 == 201912, 1, 2)
a <- a %>% rename('Date' = 'V1', 'Sample_id' = 'V2',
                  'Peptide_size' = 'V3', 'Nitrated_Peptide_size' = 'V4', 
                  'Protein_size' = 'V5')

# Save to file
a_out = paste('~/GoogleDrive/gastric_cancer_samples/Tables/table.data_OC', date, 'txt', sep='.')
write.table(a, a_out, sep='\t', row.names = F, col.names = T)
