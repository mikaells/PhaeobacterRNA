
#install the package if you dont have it
# BiocManager::install("DESeq2")
# install.packages("vegan") #for multivariate stuff
# install.packages("stringr") #for text manipulation
# install.packages("pheatmap") #for pretty heatmaps
# install.packages("readxl") #for pretty heatmaps

library(stringr)
library(DESeq2)
library(vegan)
library(pheatmap)
library(readxl)

#setup pretty plotting parameters
par(mar=c(2.5,2.5,2,0.1), mgp=c(1.4,0.4,0),font.lab=2,mfrow=c(1,1))

ReadData = function(statsPath="Data/stats.txt", annotPath="Data/annot.txt", featurePath="Data/features.txt") {
  
  #Read in the full data set
  stats=read.table(statsPath)
  annot=read.table(annotPath)
  features=read.table(featurePath)
  #KEGG=data.frame(read_excel("out.emapper.annotations.xlsx", skip=2))
  
  
  features=features[,naturalsort::naturalorder( colnames(features))]
  
  #make meta-data to keep track of samples, which we need for deseq2
  meta=data.frame(read_excel("Data/S26_transcriptome_sample_descriptions.xlsx"))
  meta$combi=factor(paste(meta$Sample.type,meta$Time,sep="_"), ordered = F, levels=c( "WT_24h", "WT_72h","dtdaB_24h", "dtdaB_72h"))
  meta$Sample.type=factor(meta$Sample.type)
  meta$Time=factor(meta$Time)
  rownames(meta) = make.unique(as.character(meta$combi))
  
  colnames(features)=meta$combi
  
  return(list(annot,features,meta))
  
}
