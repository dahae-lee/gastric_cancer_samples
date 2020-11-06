options(stringsAsFactors = F)
rm(list=ls())

library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(forcats)

outtag = 20200427
d <- read_excel('~/Analysis/gastric/Documents/80eogcH_edit.xlsx')


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
source('~/Dropbox/Resources/Colors/cohort_colors.R')

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

ggsave2(filename = paste('Figures/plot.cohort', outtag, 'png', sep='.'), 
        plot = cow, width = 11, height = 12)
