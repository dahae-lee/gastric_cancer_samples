# 데이터의 상위 폴더를 dir에 지정
dir = "/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/Protein_with_nitrated_peptide_by_sample.20200603/"

# 폴더에 존재하는 하위 폴더(=샘플명) 리스트 불러오기
# (하위 폴더가 없으면 생략)
sample_raw <- list.files(dir)

# 샘플 파일 주소 리스트 만들기(file.path = paste(..., sep = '/'))
sampleFiles <- file.path(dir, sample_raw)

# 파일 주소에 해당하는 txt파일 읽어들이기
tmp <- lapply(seq(1,65), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

# 매트릭스 하나로 bind
library(tidyverse)
a <- rbind_list(tmp)
a <- a %>% rename('Sample_ID' = 'V1', 
                  'Uniprot_ID' = 'V2',
                  'Number_whole_peptide' = 'V3',
                  'Number_nitrated_peptide' = 'V4',
                  'Proportion_nitrated_peptide' = 'V5')

# Save to file
a_out = paste('/Users/dahaelee/GoogleDrive/gastric_cancer_samples/Tables/Nitrated_sample_info.200605', 'txt', sep='.')
write.table(a, a_out, sep='\t', quote = F, row.names = F, col.names = T)
