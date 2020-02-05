date = '20200117'

# 데이터의 상위 폴더를 dir에 지정
dir_u_nit = "~/GoogleDrive/gastric_cancer_samples/Tables/GO_size.20200116./up_nitrated"
dir_d_nit = "~/GoogleDrive/gastric_cancer_samples/Tables/GO_size.20200116./down_nitrated"
dir_u_non = "~/GoogleDrive/gastric_cancer_samples/Tables/GO_size.20200116./up_nonnitrated"
dir_d_non = "~/GoogleDrive/gastric_cancer_samples/Tables/GO_size.20200116./down_nonnitrated"

# 폴더에 존재하는 하위 폴더(=샘플명) 리스트 불러오기
# (하위 폴더가 없으면 생략)
sample_raw_u_nit <- list.files(dir_u_nit)
sample_raw_d_nit <- list.files(dir_d_nit)
sample_raw_u_non <- list.files(dir_u_non)
sample_raw_d_non <- list.files(dir_d_non)


# 샘플 파일 주소 리스트 만들기(file.path = paste(..., sep = '/'))
sampleFiles_u_nit <- file.path(dir_u_nit, sample_raw_u_nit)
sampleFiles_d_nit <- file.path(dir_d_nit, sample_raw_d_nit)
sampleFiles_u_non <- file.path(dir_u_non, sample_raw_u_non)
sampleFiles_d_non <- file.path(dir_d_non, sample_raw_d_non)

# 파일 주소에 해당하는 txt파일 읽어들이기
tmp_u_nit <- lapply(seq(1,30), function(i){
  X <- read.delim(sampleFiles_u_nit[i], header = F, sep = '\t')
  X})

tmp_d_nit <- lapply(seq(1,30), function(i){
  X <- read.delim(sampleFiles_d_nit[i], header = F, sep = '\t')
  X})

tmp_u_non <- lapply(seq(1,30), function(i){
  X <- read.delim(sampleFiles_u_non[i], header = F, sep = '\t')
  X})

tmp_d_non <- lapply(seq(1,30), function(i){
  X <- read.delim(sampleFiles_d_non[i], header = F, sep = '\t')
  X})

# 매트릭스 하나로 bind
library(tidyverse)
a <- rbind_list(tmp_u_nit)
i <- seq(1,59,by=2)
a <- a[-i,]
a <- a %>% rename('Sample_ID' = 'V1', 'GO_Nitrated_Overexpressed' = 'V2')

b <- rbind_list(tmp_d_nit)
i <- seq(1,59,by=2)
b <- b[-i,]
b <- b %>% rename('Sample_ID' = 'V1', 'GO_Nitrated_Underexpressed' = 'V2')

c <- rbind_list(tmp_u_non)
i <- seq(1,59,by=2)
c <- c[-i,]
c <- c %>% rename('Sample_ID' = 'V1', 'GO_NonNitrated_Overexpressed' = 'V2')

d <- rbind_list(tmp_d_non)
i <- seq(1,59,by=2)
d <- d[-i,]
d <- d %>% rename('Sample_ID' = 'V1', 'GO_NonNitrated_Underexpressed' = 'V2')

e <- merge(a,b, by = 'Sample_ID')
e <- merge(e,c, by = 'Sample_ID')
e <- merge(e,d, by = 'Sample_ID')

e_out <- paste("~/GoogleDrive/gastric_cancer_samples/Tables/GO_size_table", date, 'txt', sep = ".")
write.table(e, e_out, sep='\t', quote = F, row.names = F, col.names = T)
