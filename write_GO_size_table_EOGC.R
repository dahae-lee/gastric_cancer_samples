# 데이터의 상위 폴더를 dir에 지정
dir = "~/GoogleDrive/gastric_cancer_samples/Tables/GO_size/overexpressed"
dir_u = "~/GoogleDrive/gastric_cancer_samples/Tables/GO_size/underexpressed"

# 폴더에 존재하는 하위 폴더(=샘플명) 리스트 불러오기
# (하위 폴더가 없으면 생략)
sample_raw <- list.files(dir)
sample_raw_u <- list.files(dir_u)

# 샘플 파일 주소 리스트 만들기(file.path = paste(..., sep = '/'))
sampleFiles <- file.path(dir, sample_raw)
sampleFiles_u <- file.path(dir_u, sample_raw_u)

# 파일 주소에 해당하는 txt파일 읽어들이기
tmp <- lapply(seq(1,30), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

tmp_u <- lapply(seq(1,30), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

# 매트릭스 하나로 bind
library(tidyverse)
a <- rbind_list(tmp)
i <- seq(1,59,by=2)
a <- a[-i,]
a <- a %>% rename('Sample_ID' = 'V1', 'GO_overexpressed' = 'V2')

b <- rbind_list(tmp_u)
i <- seq(1,59,by=2)
b <- b[-i,]
b <- b %>% rename('Sample_ID' = 'V1', 'GO_underexpressed' = 'V2')

c <- merge(a,b,by = 'Sample_ID')

c_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO_size_table", 'txt', sep = ".")
write.table(c, c_out, sep='\t', quote = F, row.names = F, col.names = T)
