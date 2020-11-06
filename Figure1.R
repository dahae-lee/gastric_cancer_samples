library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(forcats)
date = '20200613'

## Figure 1B
options(stringsAsFactors = F)
rm(list=ls())

outtag = 20200427
d <- read_excel('~/GoogleDrive/gastric_cancer_samples/80eogcH_edit.xlsx')


## Q1 : Gender
Meta1 <- d %>% group_by(Gender) %>% summarise(n = n())
Meta1 <- data.frame("Plot" = c("Gender"), "Variables"=c("Male", "Female"),
                    "Count"=c(36, 44))

## Q2 : Age
Meta2 <- as.data.frame(0)
colnames(Meta2) <- "Age"
Meta2[c(1),] <- length(which(d$`Age at the time of surgery`>= 20 & d$`Age at the time of surgery` < 30))
Meta2[c(2),] <- length(which(d$`Age at the time of surgery`>= 30 & d$`Age at the time of surgery` < 40))
Meta2[c(3),] <- length(which(d$`Age at the time of surgery`>= 40))

Meta2 <- data.frame("Plot" = c("Age"), "Variables"= c("20~29", "30~39", ">= 40"),
                    "Count"=c(7, 41, 32))

## Q3 : Location
Meta3 <- d %>% group_by(`Location of primary tumor in the stomach`) %>% summarise(Count = n())
Meta3 <- Meta3[c(4,3,1,2),]
colnames(Meta3) <- c("Variables", "Count")
Meta3 <- transform(Meta3, Variables = factor(Variables, 
                                             levels = c("Proximal 1/3", "Middle 1/3", "Distal 1/3", "Entire")))
Meta3 <- mutate(Meta3, Plot = c("Location"))
Meta3 <- Meta3[,c(3,1,2)]

## Q4 : Stage
Meta4 <- d %>% group_by(Stage) %>% summarise(Count = n())
colnames(Meta4) <- c("Variables", "Count")
Meta4 <- transform(Meta4, Variables = factor(Variables, 
                                             levels = c("IA", "IB", "IIA", "IIB", "IIIA", "IIIB", "IIIC", "IV")))
Meta4 <- mutate(Meta4, Plot = c("Stage"))
Meta4 <- Meta4[,c(3,1,2)]

## Draw figures

Meta <- rbind(Meta1, Meta2, Meta3, Meta4)
Meta$Plot <- factor(Meta$Plot, levels = c("Stage", "Location", "Age", "Gender"))
Meta$Variables <- fct_rev(factor(Meta$Variables, levels = c("Male", "Female", "20~29", "30~39", ">= 40", 
                                                            "Proximal 1/3", "Middle 1/3", "Distal 1/3", "Entire", 
                                                            "IA", "IB", "IIA", "IIB", "IIIA", "IIIB", "IIIC", "IV")))
source('~/GoogleDrive/gastric_cancer_samples/gastric_cancer_samples/cohort_colors.R')

p <- ggplot(Meta, aes(Plot, Count, fill = Variables))+
  geom_bar(stat='identity', width = 0.9)+ 
  labs(x = "", y = "", shape = "Cohort")+ 
  scale_fill_manual(values = cohort_colors)+
  coord_flip(expand = F)+ 
  theme(axis.text.y = element_text(angle=0, face="bold", size=16))+
  theme(axis.text.x = element_text(size=12))+
  theme(legend.position = "none")+
  theme(panel.background = element_blank())

p1 <- Meta %>%
  filter(Plot %in% c("Gender", "Age")) %>%
  ggplot(aes(Plot, Count, fill = Variables))+
  geom_bar(stat='identity')+ 
  scale_fill_manual(values = cohort_colors, name = "Cohort", 
                    breaks = c("Male", "Female", "20~29", "30~39", ">= 40"))+
  theme(legend.title = element_text(color = "black", size = 15))+ 
  theme(legend.text = element_text(color = "black", size = 13))

p2 <- Meta %>%
  filter(Plot %in% c("Location")) %>%
  ggplot(aes(Plot, Count, fill = Variables)) + geom_bar(stat='identity')+ 
  scale_fill_manual(values = cohort_colors, name = "", 
                    breaks = c("Proximal 1/3", "Middle 1/3", "Distal 1/3", "Entire"))+
  theme(legend.title = element_text(color = "black", size = 10))+ 
  theme(legend.text = element_text(color = "black", size = 13))  

p3 <- Meta %>%
  filter(Plot %in% c("Stage")) %>%
  ggplot(aes(Plot, Count, fill = Variables)) + geom_bar(stat='identity')+ 
  scale_fill_manual(values = cohort_colors, name = "", 
                    breaks = c("IA", "IB", "IIA", "IIB", "IIIA", "IIIB", "IIIC", "IV"))+
  guides(fill = guide_legend(nrow=4))+
  theme(legend.title = element_text(color = "black", size = 10))+ 
  theme(legend.text = element_text(color = "black", size = 13))  

## Merge figures :: library(cowplot)
cow <- plot_grid(
  p, 
  plot_grid(get_legend(p1),
            get_legend(p2),
            get_legend(p3),
            nrow = 1, rel_widths = c(1,1,1)),
  nrow = 2, rel_heights = c(8,2))

## Figure 1C
# Load data
d = read.delim('~/GoogleDrive/gastric_cancer_samples/Tables/table.data_QC.20200421.txt')

# Remove duplicated samples
c <- as.data.frame(table(d$Sample_id))
dup <- c %>% filter(Freq != 1) %>% pull(Var1) %>% as.character()
d1 <- d %>% filter(Sample_id %in% dup)
d2 <- d %>% filter(!(Sample_id %in% dup & Date == 201912))

## Draw Peptide size plot
pa <- d2 %>%
  ggplot(aes(x = Sample_id, y = Peptide_size, fill = Sample_id))+
  geom_bar(stat = "identity", show.legend = FALSE)+
  scale_colour_viridis_d(alpha = 1, begin = 0, end = 1,
                         direction = 1, option = "D", aesthetics = "fill")+
  xlab("Sample ID")+
  ylab("Peptide Size")+
  ggtitle("Peptide Size by Samples")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank())

## Draw Nitrated Peptide size plot
pb <- d2 %>% 
  ggplot(aes(x = Sample_id, y = Nitrated_Peptide_size, fill = Sample_id))+
  geom_bar(stat = "identity", show.legend = FALSE)+
  scale_colour_viridis_d(alpha = 1, begin = 0, end = 1,
                         direction = 1, option = "D", aesthetics = "fill")+
  xlab("Sample ID")+
  ylab("Nitrated Peptide Size")+
  ggtitle("Nitrated Peptide Size by Samples")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank())

## Draw Protein size plot
# 데이터의 상위 폴더를 dir에 지정
dir = "~/GoogleDrive/gastric_cancer_samples/Tables/total_protein_expression.20200421"
date = '20200421'

# 폴더에 존재하는 하위 폴더(=샘플명) 리스트 불러오기
# (하위 폴더가 없으면 생략)
sample_raw <- list.files(dir)

# 샘플 파일 주소 리스트 만들기(file.path = paste(..., sep = '/'))
sampleFiles <- file.path(dir, sample_raw)

# 파일 주소에 해당하는 txt파일 읽어들이기
tmp <- lapply(seq(1,length(sampleFiles)), function(i){
  X <- read.delim(sampleFiles[i], header = F, sep = '\t')
  X})

# 매트릭스 하나로 bind
library(tidyverse)
a <- rbind_list(tmp)
a$batch <- ifelse(a$V1 == 201912, '1', '2')
a <- a %>% rename('Date' = 'V1', 'Sample_ID' = 'V2',
                  'norm1' = 'V3', 'tumor1' = 'V4', 
                  'norm2' = 'V5', 'tumor2' = 'V6')


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
a1 <- a %>% filter(!(Sample_ID %in% dup & Date == 20200309))

# make plots
pc <- a1 %>%
  ggplot(aes(x = reorder(Sample_ID,Date), y = log(expression), fill = batch))+
  geom_bar(stat = "identity")+
  scale_colour_viridis_d(alpha = 1, begin = 0, end = 1,
                         direction = 1, option = "D", aesthetics = "fill")+
  facet_grid(celltype~.)+
  xlab("Sample ID")+
  ylab("Total Protein Expression Size")+
  ggtitle("Total Protein Expression Size by Samples")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank())

## Merge plots
c <- plot_grid(pa, pb, pc, ncol=1)
p <- plot_grid(cow, c, ncol = 2, nrow = 1, labels = c('B', 'C'))
p <- ggdraw()+
  draw_plot(pa, x = 0, y = 2, width = 1, height = 1)+
  draw_plot(pb, x = 0, y = 1, width = 1, height = 1)+
  draw_plot(pc, x = 0, y = 0, width = 1, height = 1)+
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0), y = c(2, 1, 0.5))

# with A
p <- ggdraw() +
  draw_plot(cow, x = 0, y = 0.05, width = .4, height = 0.5) +
  draw_plot(c, x = .5, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0.5), y = c(1, 0.55, 1))
# without A
p <- ggdraw() +
  draw_plot(cow, x = 0, y = .3, width = .4, height = 0.7) +
  draw_plot(c, x = .5, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("B", "C"), size = 15,
                  x = c(0, 0.5), y = c(1, 1))
p

ggsave('Figures/Fig.1.200615.pdf', p, width = 10, height = 8)

