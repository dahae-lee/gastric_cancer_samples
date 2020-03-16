options(stringsAsFactors = F)
rm(list=ls())
library(WGCNA)
library(preprocessCore)
library(tidyverse, quietly = T)
library(RColorBrewer)
library(fields)
library(ggplot2)
source('~/Dropbox/Scripts/rnaseq/function_wgcna.R')

## Load data: normalized counts
option_merge = F

if (option_merge){
  fs = list.files('Inputs', pattern='norm', full.names = T)
  
  sid = strsplit(fs[1], '.', fixed = T)[[1]][3]
  d = read.delim(fs[1]) 
  colnames(d)[2:5] <- paste(colnames(d)[2:5], sid, sep='_')
  
  for (f in fs[2:length(fs)]){
    print (f)
    sid = strsplit(f, '.', fixed = T)[[1]][3]
    d1 = read.delim(f)
    colnames(d1)[2:5] <- paste(colnames(d1)[2:5], sid, sep='_')
    
    d = merge(d, d1, by='gene_id')
  }
  write.table(d, 'table.merged_expr_norm.20200311.txt', sep='\t', quote = F, row.names = F, col.names = T)
} else {
  d = read.delim('table.merged_expr_norm.20200311.txt')
}

## Prepare inputs
e = as.matrix(d %>% select(-gene_id))
e = normalize.quantiles(as.matrix(e)) # quantile normalization
rownames(e) <- d$gene_id 

gsg = goodSamplesGenes(t(e), verbose = 5)
table(gsg$goodGenes) # just check whether any bad gene
table(gsg$goodSamples)  # just check whether any bad sample


outtag = 'all_20200311'

## Determine soft power cutoff for WGCNA
option_determineSoftPowerWGCNA = T
if (option_determineSoftPowerWGCNA){
  sft0 <- determineSoftPowerWGCNA(
    data1 = e, 
    outtag = outtag,
    sft1=1, sft2=30, sft_by=1
  )
  
  sft <- sft0$fitIndices
  sft_cut = min(sft[sft$SFT.R.sq >= 0.8 & sft$mean.k. <= 50,]$Power)
  print (paste('The sft cutoff for this analysis is', sft_cut, sep=' '))
}

## Set a run parameter for WGCNA
type <- 'signed'
type2 <- ifelse(type == 'signed hybrid', 'hybrid', type)
minModuleSize <- 20
deepSplit <- 4

outtag = paste(outtag, 
               'minSize', minModuleSize,
               'type',  type2,
               'ds', deepSplit,
               sep='_')

f_oriMat = paste('Tables/table.wgcna_original_module_genes', outtag, 'txt', sep='.')
f_ME = paste('Tables/table.ME', outtag, 'txt', sep='.')
f_kME = paste('Tables/table.kME', outtag, 'txt', sep='.')

## Adjacency matrix
TOMsim <- TOMsimilarityFromExpr(t(e), 
                                corType = "pearson", 
                                networkType = type, 
                                power = sft_cut, 
                                TOMType = "signed", 
                                nThreads = 4)

dissTOMA1 <- 1 - TOMsim
geneTree = hclust(as.dist(dissTOMA1), method="average")

tree = cutreeHybrid(dendro = geneTree,
                    minClusterSize= minModuleSize,
                    pamStage=FALSE,
                    cutHeight = 0.999,
                    deepSplit= deepSplit,
                    distM= dissTOMA1 )

merged = mergeCloseModules(exprData= t(e),
                           colors = tree$labels,
                           cutHeight=0.10,
                           verbose = 3)

## Merge information into the data frame
res = setNames(as.data.frame(cbind(tree$labels, merged$colors)), c('Raw', 'Merged'))
res$gid <- rownames(e)

res$Merged_color = labels2colors(res$Merged)
sort(table(res$Merged_color))
length(unique(res$Merged_color))

# write.table(res, f_oriMat, sep='\t', quote=F, col.names = T, row.names = F)
moddat = merged$newMEs
rownames(moddat) <- colnames(e)
# write.table(moddat, f_ME, sep='\t', quote=F, col.names = NA, row.names = T)

## Plot WGCNA modules
f_plot = paste('Figures/plot.wgcna_original_module_genes', outtag, 'png', sep='.')
png(f_plot, width=1600, height=400)
plotDendroAndColors(geneTree,
                    res$Merged_color,
                    "TCGA-LUAD-TP",
                    addGuide = F,
                    colorHeightMax = 0.9,
                    dendroLabels = FALSE)
dev.off()
