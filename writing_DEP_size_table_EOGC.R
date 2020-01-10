# 데이터의 상위 폴더를 dir에 지정
dir = "~/GoogleDrive/gastric_cancer_samples/Tables/sample_info"

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
a <- a %>% rename('Sample_ID' = 'V1', 'Total_Protein' = 'V2', 'DEP_Peptide_Up' = 'V3', 
                  'DEP_Peptide_Down' = 'V4', 'DEP_Protein_Up' = 'V5', 
                  'DEP_Protein_Down' = 'V6', 'NonDEP_Nitrated_Protein' = 'V7',
                  'NonDEP_NonNitrated_Protein' = 'V8')

# Save to file

a_out = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/DEP_sample_info', 'txt', sep='.')
write.table(a, a_out, sep='\t', quote = F, row.names = F, col.names = T)
