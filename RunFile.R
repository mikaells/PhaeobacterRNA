#clean out existing stuff to start fresh
rm(list = ls())

#0. SETUP THE LIBRARIES

#Uncomment if you dont have the packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("Rsubread")

#read in the libraries we need
library("Rsubread")


#####
#1. SETUP THE SEQUENCING FILES 
#####

#get the directory of the raw files
fastq_directory <- "inputs/"
#get the adresses of the files we need
#remember that all samples have two files
reads1 <- list.files(path = file.path(fastq_directory), pattern = "*1.fq.gz$", full.names = TRUE)
reads2 <- list.files(path = file.path(fastq_directory), pattern = "*2.fq.gz$", full.names = TRUE)

#do we have all 12 files? Should return 12 for both
cat(paste("reads1 =",length(reads1),"\nreads2 =", length(reads2)))


#####
#2. MAKE THE INDEX TO ALIGN AGAINST
#####

#first, we make a folder/directory for the index
dir.create("rsub_indx")

#Then we build the actual index from the whole genome file
#and we will put it in the indx/-folder we just made
buildindex(basename="rsub_indx/PROKKA_03102021", #the name of our index
           reference="PGAP_annot/Phaeobacter_S26_PGAP_220208.fa", #what file to make it from?
           memory=16000, #max memory is surely enough for a smallish bacterial genome
           gappedIndex = F) #?

#####
#3. MAKE THE MAP BETWEEN POSITION AND GENE
#####

# now we map the reference
#The .gff3 file has all that info, so we will turn it into a SAF-file and use that
saf=flattenGTF(
  GTFfile =  "PGAP_annot/Phaeobacter_S26_PGAP_220208.gff", 
  GTF.featureType = "CDS",
  GTF.attrType = "product",
  # the option specifying the merging algorithm
  method = "merge")

gff=readLines("PGAP_annot/Phaeobacter_S26_PGAP_220208.gff")[-c(1:7)]

gff2=stringr::str_split_fixed(gff,"\t", 10)

gff3=gff2[grep(gff2[,3],pattern = "CDS"  ),  ] 

CDS=stringr::str_split_fixed(gff3[,9]  , ";", 10)
length(grep("hypothetic",saf$GeneID ))


saf2=saf[order(saf$Start),]

dubliNames=names(which (table(saf2$GeneID  ) >1))

for(j in dubliNames){ 
  counter=1
  for(i in 1:NROW(saf2)) {
    
    if(saf2$GeneID[i] ==j ){
      newName=paste(saf2$GeneID[i] , counter,sep="_") 
      
      saf2$GeneID[i] = newName
      #print(newName)
      counter=counter+1
    }  
  }
}

saf2[ grep("phage", saf2$GeneID ),] 

#####
#4. ALIGN THE SEQUENCES TO INDEX
#####

#we now align our sequencing files to the reference
#e.g. each read pair will be aligned to the genome
#note that the alignments are saved as external files rather than as R-objects
#****
#LOOK AT THE nthreads option
#if your computer only has 4 threads, its gonna die with 12. Set it at your threads minus 1.
#****
#will be fast on a linux (perhaps on a mac too) and not on a windows (something weird with the 
#index-load on some machines)
#it is working though, note bam-files appearing in the subset-folder

#make a folder for the bam-files
dir.create("bam_out")

#then the command for alignment
alignInfo=align(index = "rsub_indx/PROKKA_03102021",   #the index
                readfile1 = reads1,            #read pairs 1
                readfile2 = reads2,            #read pairs 2
                type = "rna",                  #we are doing rna
                input_format = "gzFASTQ",        #files are fastq-files, could also be fastq.gz
                output_format = "BAM",         #output will be BAM files
                PE_orientation = "fr",         #reads are [f]orward and [r]everse
                nthreads = 20,                  #number of threads
                output_file = gsub("_1.fq.gz",".bam",gsub("inputs","bam_out",reads1)) 
)


#####
#5.  COUNT HOW THE SEQUENCES MAPPED TO THE GENES
#####

#lets get the list of BAM files
bam.files <- list.files("bam_out/", pattern = "bam$", full.names = T)

#Run the actual counting.
fc=featureCounts(files = bam.files,      #the bam-files we are counting 
                 annot.ext = saf2,nthreads=20,   #the external object we are annotating with
                 isPairedEnd = T)   #and our files are obviously paired


stats=fc$stat

features=fc$counts

annot=fc$annotation 

write.table(stats,"stats.txt")
write.table(features,"features.txt")
write.table(annot,"annot.txt")


#features[order(features, decreasing = T)]   
