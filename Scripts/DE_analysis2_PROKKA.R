rm(list=ls()); gc()

writeFig=F
writeTable=F
doAllPlots=T

source("Scripts/functions.R")

if(writeFig){
  pdf(file = "Figures//allFigs.pdf", onefile = T, width = 16,height = 12)
  par(mar=c(2.5,2.5,2,0.1), mgp=c(1.4,0.4,0),font.lab=2,mfrow=c(1,1))
}
#####
#0) read in the files and set up the meta

All=ReadData(statsPath = "PROKKA_data/stats.txt", annotPath = "PROKKA_data/annot.txt", featurePath = "PROKKA_data/features.txt")

annot=All[[1]]
features=All[[2]]
meta=All[[3]]

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
#3) Multivariate stuff
#####

plot(0,0, col=0,xaxt='n',yaxt='n',frame.plot=F, xlab="",ylab="", ann=F)
text(0,0,"Multivariate", cex=3, xaxt='n')

#first a heatmap
pheatmap(log10(counts(dds, normalized=T)+1), labels_row = "")

#make PCA using transformed reads
#Transformation is basically a log of normalized counts, i think
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="combi")

#Finding co-occuring genes
COR=cor(t((counts(dds, normalized=T))),method = "spear")
pheatmap(COR[1:200,1:200], labels_col = NA,labels_row = NA, cutree_rows = 5, cutree_cols = 5)

#####
#4) individual genes
#####

#the TDA cluster
plot(0,0, col=0,xaxt='n',yaxt='n',frame.plot=F, xlab="",ylab="", ann=F)
text(0,0,"TDA-Cluster", cex=3, xaxt='n')

plotCounts(dds, gene="LysR family transcriptional regulator_3",     intgroup="combi", pch=16)
plotCounts(dds, gene="LysR family transcriptional regulator_2",     intgroup="combi", pch=16) #A
plotCounts(dds, gene="glutathione S-transferase family protein_1",  intgroup="combi", pch=16) #B
plotCounts(dds, gene="prephenate dehydratase_1",                    intgroup="combi", pch=16) #C
plotCounts(dds, gene="acyl-CoA thioesterase_2",                     intgroup="combi", pch=16) #D
plotCounts(dds, gene="prephenate dehydratase_1",                    intgroup="combi", pch=16) #E


#####
#5) Differential expression
#####

#lets compare "WT_72h" with"dtdaB_72h"
#note notation for the contrast argument is first the variable and next the groups to compare  
res_WT72_dtdaB72 <- results(object = dds,contrast = c("combi", "WT_72h","dtdaB_72h"))

summary(res_WT72_dtdaB72)

#Make overview of all genes, their coordinates and if they are differential
AllCoord0=data.frame(annot[match(rownames(res_WT72_dtdaB72),annot$GeneID),], res_WT72_dtdaB72)
AllCoord=AllCoord0[order(AllCoord0$Chr,AllCoord0$Start),]

#calculating the distances between genes
AllCoord$gap=-1

for(i in 1:(NROW(AllCoord)-1)) {
  AllCoord$gap[i]=AllCoord$Start[i+1] - AllCoord$End[i]
  if(AllCoord$Chr[i+1]!=AllCoord$Chr[i]) {
    AllCoord$gap[i]=0
  }
}

#write to file
if(writeTable) {
  write.csv(x = AllCoord, file = "PROKKA_summaries//allGenes.csv", row.names = F)
}

#Make a nice overview of architecture of significant genes
#organizing  annotation by the significant genes
SigCoord=subset(AllCoord, padj<0.05 & abs(log2FoldChange)>1)

#write to file
if(writeTable) {
  write.csv(x = sigCoord, file = "PROKKA_summaries/sigGenes.csv", row.names = F)
}

######
#looking for phages
#####

#Phage1
phage1=subset(subset(AllCoord, Chr=="1"), (Start>626399-10) & (End <657712+10) )

par(mfrow=c(5,5))
plot(0,0, col=0,xaxt='n',yaxt='n',frame.plot=F, xlab="",ylab="", ann=F)
text(0,0,"Phage1", cex=3, xaxt='n')
for(i in phage1$GeneID) {
  if(doAllPlots) plotCounts(dds, gene=i, intgroup="combi", col=c(rep(1,9),2,1,1))
}

#Phage2
phage2=subset(subset(AllCoord, Chr=="1"), (Start>1889854) & (End <1905516+10) )

par(mfrow=c(5,5))
plot(0,0, col=0,xaxt='n',yaxt='n',frame.plot=F, xlab="",ylab="", ann=F)
text(0,0,"Phage2", cex=3, xaxt='n')
for(i in phage2$GeneID) {
  if(doAllPlots) plotCounts(dds, gene=i, intgroup="combi")
}

####
#Look at metal genes
####

MetalCoord=AllCoord[grep("Zn|zinc|ferro|iron|copper", AllCoord$GeneID,ignore.case = T),]

MetalCoordSig=subset(MetalCoord, padj<.05 & abs(log2FoldChange)>.1)
#write to file
if(writeTable) {
  write.csv(metalCoord,"Summaries/metals.csv" )
}

#plot"
par(mfrow=c(5,5))
plot(0,0, col=0,xaxt='n',yaxt='n',frame.plot=F, xlab="",ylab="", ann=F)
text(0,0,"Metals", cex=3, xaxt='n')
for(i in MetalCoordSig[order(MetalCoordSig$padj),]$GeneID) {
  if(doAllPlots) plotCounts(dds, gene=i, intgroup="combi" ,col=c(rep(1,9),2,1,1))
}

par(mfrow=c(1,1))

#####
#6) Volcano plot
#####


#plot all points, log2fold vs -log10(pvalues)
with(AllCoord, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot" , xlim=c(-8,8)))
arrows(x0 = 7,x1 = 8,y0 = 40,y1 = 40,length = .1)
text(7.5,50, "tdaB, lfc>16")
#make all significant points blue
with(subset(AllCoord, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
#make all signifcant AND at least 2-fold differential red
with(subset(AllCoord, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#with(subset(AllCoord, padj<.1E-20 & abs(log2FoldChange)>5 ), text(log2FoldChange, -log10(pvalue),GeneID,))
with(phage2, points(log2FoldChange, -log10(pvalue), pch=21, bg="green"))
with(MetalCoordSig, points(log2FoldChange, -log10(pvalue), pch=21, bg="grey"))

legend("topleft",legend = c( "p<0.05","lfc>2","GTA","Metal"), pt.bg = c("blue", "red", "green", "grey"), pch=21)

######
#Print genome context
######

plot(0,0, col=0,xaxt='n',yaxt='n',frame.plot=F, xlab="",ylab="", ann=F)
text(0,0,"Gene architecture", cex=3, xaxt='n')

j=sort(unique(AllCoord$Chr))[1]
for(j in sort(unique(AllCoord$Chr))) {
  
  ChrCoord=subset(AllCoord, Chr==j)
  
  with(ChrCoord,plot(x=Start, abs(log2FoldChange),main=j, ylab="log2-change", pch=16,col=0))
  abline(v = seq(from=0, to=max(ChrCoord$End),length.out=20), col="grey90")
  
  if(j=="CP080275") {
    lines(x=c(626399, 657712), y=c(-.1,-.1), col=5, lwd=5)
    lines(x=c(1889854, 1905516), y=c(-.1,-.1), col=5, lwd=5)
    rect(xleft = 626399, ybottom = 0,xright = 657712,ytop = max(ChrCoord$log2FoldChange), col = "grey90", lty = 0)
    rect(xleft = 1889854, ybottom = 0,xright = 1905516,ytop = max(ChrCoord$log2FoldChange), col = "grey90", lty = 0)
  }
  
  with(ChrCoord,points(x=Start, abs(log2FoldChange),main=j, ylab="log2-change", pch=16, cex=ifelse(padj<0.03,.8,.1), col=ifelse(padj<0.01,2,1)))
  with(subset(MetalCoordSig, Chr==j), points(x=Start, abs(log2FoldChange), cex=1, pch=21, bg=4))
}  

if(writeFig){
  dev.off() 
}

