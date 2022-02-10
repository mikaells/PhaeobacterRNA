rm(list=ls()); gc()

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

#####
#0) read in the files and set up the meta

#Read in the full data set
stats=read.table("stats.txt")
annot=read.table("annot.txt")
features=read.table("features.txt")
KEGG=data.frame(read_excel("out.emapper.annotations.xlsx", skip=2))



features=features[,naturalsort::naturalorder( colnames(features))]

#make meta-data to keep track of samples, which we need for deseq2
meta=data.frame(read_excel("S26_transcriptome_sample_descriptions.xlsx"))
meta$combi=factor(paste(meta$Sample.type,meta$Time,sep="_"), ordered = F, levels=c( "WT_24h", "WT_72h","dtdaB_24h", "dtdaB_72h"))
meta$Sample.type=factor(meta$Sample.type)
meta$Time=factor(meta$Time)
rownames(meta) = make.unique(as.character(meta$combi))

colnames(features)=meta$combi

#####
#1) make deseq object 
#####

#put the gene names in first column for compatibility with DESeq2
features_DS=data.frame(rownames(features),features)


#turn into a DESeq2 object
dds_dat <- DESeqDataSetFromMatrix(countData=features_DS, 
                                  colData=meta, 
                                  design=~combi, tidy = TRUE)

#Run the DESeq analysis
dds <- DESeq(dds_dat)


#####
#2) Have a look at the DESeq object
#####

## Check the size factors
sizeFactors(dds)

## Total number of raw counts per sample
plot(colSums(counts(dds)))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Dispersion plots like in the paper
plotDispEsts(dds)

#####
#3) Multivariate stuff
#####

#first a heatmap
#note how we can fish out the counts from the dds object with counts()
#we are overwriting the labels of the columns, because the names are so long
pheatmap(log10(counts(dds, normalized=T)+1), labels_row = "")


#make PCA using transformed reads
#Transformation is basically a log of normalized counts, i think
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="combi")

#####
#4) individual genes
#####

#the TDA cluster

plotCounts(dds, gene="LysR family transcriptional regulator_3",  intgroup="combi", pch=16)
plotCounts(dds, gene="LysR family transcriptional regulator_2",  intgroup="combi", pch=16) #A
plotCounts(dds, gene="glutathione S-transferase family protein_1",  intgroup="combi", pch=16) #B
plotCounts(dds, gene="prephenate dehydratase_1",  intgroup="combi", pch=16) #C
plotCounts(dds, gene="acyl-CoA thioesterase_2",  intgroup="combi", pch=16) #D
plotCounts(dds, gene="prephenate dehydratase_1",  intgroup="combi", pch=16) #E




#####
#5) Differential expression
#####

#lets compare "WT_72h" with"dtdaB_72h"
#note notation for the contrast argument is first the variable and next the groups to compare  
res_WT72_dtdaB72 <- results(object = dds,contrast = c("combi", "WT_72h","dtdaB_72h"))

#look at the summary
summary(res_WT72_dtdaB72)

#order by log2FoldChange rather than alphabetic, then the largest wil be on top
#ordering by p-value might also make sense
#recall that res_WT72_dtdaB72 behaves like a matrix and if you can re-order the rows like this
res_WT72_dtdaB72_sig=res_WT72_dtdaB72[order(res_WT72_dtdaB72$padj),]


#plot all significant and large-ish differences
sigLarge=subset(res_WT72_dtdaB72_sig, padj<.05 & abs(log2FoldChange)>1)



#Make a nice overview of architecture of significant genes
#organizing  annotation by the significant genes
sigCoord0=data.frame(annot[match(rownames(sigLarge),annot$GeneID),], sigLarge)
sigCoord=sigCoord0[order(sigCoord0$Chr,sigCoord0$Start),]

#calculating the distances between genes
sigCoord$gap=-1

for(i in 1:(NROW(sigCoord)-1)) {
  sigCoord$gap[i]=sigCoord$Start[i+1] - sigCoord$End[i]
  if(sigCoord$Chr[i+1]!=sigCoord$Chr[i]) {
    sigCoord$gap[i]=0
  }
  
}

#Printing all significant genes

rownames(sigLarge)
par(mfrow=c(5,5))
for(i in rownames(sigLarge)) {
  if(grepl("hypothetical",i)) next
  plotCounts(dds, gene=i, intgroup="combi" )
}


#looking for phages

phageAnnot1=subset(subset(annot, Chr=="CP080275"), (Start>626399-10) & (End <657712+10) )

phageRes=res_WT72_dtdaB72[rownames(res_WT72_dtdaB72) %in% phageAnnot1$GeneID,] 

par(mfrow=c(5,5))
for(i in rownames(phageRes)) {
  plotCounts(dds, gene=i, intgroup="combi")
}


phageAnnot2=subset(subset(annot, Chr=="CP080275"), (Start>1889854) & (End <1905516+10) )

phageRes2=res_WT72_dtdaB72[rownames(res_WT72_dtdaB72) %in% phageAnnot2$GeneID,] 

par(mfrow=c(5,5))
for(i in rownames(phageRes2)) {
  plotCounts(dds, gene=i, intgroup="combi")
}



metal=res_WT72_dtdaB72_sig[grep("metal|Zn|zinc|ferro|iron|copper",rownames(res_WT72_dtdaB72_sig),ignore.case = T),]
metalSigLarge=subset(metal, padj<.05 & abs(log2FoldChange)>.1)

par(mfrow=c(3,4))
for(i in rownames(metalSigLarge)) {
  plotCounts(dds, gene=i, intgroup="combi" )
}

par(mfrow=c(1,1))



#####
#6) Volcano plot
#####

par(mfrow=c(1,1))
#plot all points, log2fold vs -log10(pvalues)
with(res_WT72_dtdaB72, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot" ))
#make all significant points blue
with(subset(res_WT72_dtdaB72, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
#make all signifcant AND at least 2-fold differential red
with(subset(res_WT72_dtdaB72, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

######
#Print genome context
######
j=sort(unique(annot$Chr))[2]
for(j in sort(unique(annot$Chr))) {
  
  annotChr=subset(x = annot, Chr==j )
  
  res_WT72_dtdaB72_Chr=res_WT72_dtdaB72[rownames(res_WT72_dtdaB72) %in% annotChr$GeneID,]
  
  #plot(x=annotChr$Start, -log10(res_WT72_dtdaB72_Chr@listData$padj), main=j, cex=0.2, pch=16, col=ifelse(res_WT72_dtdaB72_Chr@listData$padj<0.01,2,1))
  plot(x=annotChr$Start, abs(res_WT72_dtdaB72_Chr@listData$log2FoldChange),main=j, ylab="log2-change", cex=0.5, pch=16, col=ifelse(res_WT72_dtdaB72_Chr@listData$padj<0.01,2,1))
  
  metChr=annotChr[grepl("metal|Zn|zinc|ferro|iron|copper", annotChr$GeneID),]
  
  metChr_res=res_WT72_dtdaB72_Chr[metChr$GeneID,]
  
  points(x=metChr$Start, -log10(metChr_res@listData$padj), cex=0.8, pch=23, bg=4)
  
  
  if(j=="CP080275") {
    lines(x=c(626399, 657712), y=c(-.1,-.1), col=2, lwd=3)
    lines(x=c(1889854, 1905516), y=c(-.1,-.1), col=2, lwd=3)
  }
  
  
  
  if(j=="CP080276") { 
    res_WT72_dtdaB72_Chr[c((grep("glutathione S-transferase family protein_1", rownames(res_WT72_dtdaB72_Chr))-5):(grep("glutathione S-transferase family protein_1", rownames(res_WT72_dtdaB72_Chr))+1)),]
    TDA=annotChr[c((grep("glutathione S-transferase family protein_1", annotChr$GeneID)-5):(grep("glutathione S-transferase family protein_1", annotChr$GeneID)+1)),]
    lines(x=c(min(TDA$Start),max(TDA$End)),y=c(0,0),col=3, lwd=3)
    lines(x=c(min(TDA$Start),max(TDA$End)),y=c(max(res_WT72_dtdaB72_Chr$log2FoldChange)-1,max(res_WT72_dtdaB72_Chr$log2FoldChange)-1),col=3, lwd=3)
  }
}  

