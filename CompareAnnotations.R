rm(list=ls())

library(Rsubread)

safA=flattenGTF(
  GTFfile =  "PROKKA_03102021_correct_locus_tag.gff", 
  GTF.featureType = "CDS",
  GTF.attrType = "product",
  # the option specifying the merging algorithm
  method = "merge")


safA2=safA[order( safA$Chr, safA$Start),]
safA2$Length=safA2$End - safA2$Start



safB=flattenGTF(
  GTFfile =  "Phaeobacter_S26_PGAP_220208.gff", 
  GTF.featureType = "CDS",
  GTF.attrType = "product",
  # the option specifying the merging algorithm
  method = "merge")

safB2=safB[order(safB$Chr,safB$Start),]
safB2$Length=safB2$End - safB2$Start


length(which(safA2$GeneID=="hypothetical protein"))
length(which(safB2$GeneID=="hypothetical protein"))


SafA_chr1=subset(safA2, Chr=="1")
SafB_chr1=subset(safB2, Chr=="CP080275")
both_chr1=merge(x = SafA_chr1 ,y =  SafB_chr1, by.x="Start", by.y="Start", all=T)

