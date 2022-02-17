rm(list=ls())

library(Rsubread)

safA=flattenGTF(
  GTFfile =  "Data/PROKKA_03102021_correct_locus_tag.gff", 
  GTF.featureType = "CDS",
  GTF.attrType = "locus_tag",#"product",
  # the option specifying the merging algorithm
  method = "merge")


safA2=safA[order( safA$Chr, safA$Start),]
safA2$Length=safA2$End - safA2$Start



safB=flattenGTF(
  GTFfile =  "Data/Phaeobacter_S26_PGAP_220208.gff", 
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

both_chr_all=both_chr1
for(i in 2:5) {
  SafA_chr=subset(safA2, Chr==unique(safA2$Chr)[i])
  SafB_chr=subset(safB2, Chr==unique(safB2$Chr)[i])
  
  both_chr=merge(x = SafA_chr ,y =  SafB_chr, by.x="Start", by.y="Start", all=T)
  
  both_chr_all=rbind(both_chr_all, both_chr )
}

write.csv(x = both_chr_all, file = "BothAnnotations.csv", row.names = F)


