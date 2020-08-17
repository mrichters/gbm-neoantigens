# 9.20.19

library(sciClone)
library(readxl)
library(tidyverse)
library(ggsci)
library(plyr); library(dplyr)
library(GGally)
library(UpSetR)
library(gridExtra)
library(gtable)
library(scales)

# browser() command will stop the code at specific points

#set wd
setwd("~/Desktop/GBM/GBM065/sciclone")

#read in folders with 2 required file types
vaf_files <- unlist(list.files("vaf_files"))
#cn_files <- unlist(list.files("sciclone_seg_files"))

# create vector with all sample / tumor names 
sample_names <- unlist(lapply(vaf_files, function(x) unlist(strsplit(x, split = ".tsv"))[1]))
tumor_names <- "GBM065" #unique(unlist(lapply(sample_names, function(x) unlist(strsplit(x, split = "-"))[2])))

# loop over sample / tumor names and run sciclone
for ( tumor in tumor_names ) { #2:length(tumor_names)] ) {
  samples <- sample_names[grep(tumor, sample_names)]
  if ( length(samples) == 2 ) {
    # define files for input into sciclone
    v1 = read.table(paste("vaf_files/",samples[1], ".tsv", sep=""), header=T)
    v2 = read.table(paste("vaf_files/",samples[2], ".tsv", sep=""), header=T)
    cn1 = read.table(paste("sciclone_seg_files/",samples[1], ".mod.seg", sep=""), header=T)
    cn2 = read.table(paste("sciclone_seg_files/",samples[2], ".mod.seg", sep=""), header=T)
    # clustering algorithm:
    sc = sciClone(vafs=list(v1,v2), copyNumberCalls=list(cn1,cn2), sampleNames=samples, minimumDepth = 50, cnCallsAreLog2=TRUE)
  } else if ( length(samples) == 3) {
    # define files for input into sciclone
    v1 = read.table(paste("vaf_files/",samples[1], ".tsv", sep=""), header=T)
    v2 = read.table(paste("vaf_files/",samples[2], ".tsv", sep=""), header=T)
    v3 = read.table(paste("vaf_files/",samples[3], ".tsv", sep=""), header=T)
    cn1 = read.table(paste("sciclone_seg_files/",samples[1], ".mod.seg", sep=""), header=T)
    cn2 = read.table(paste("sciclone_seg_files/",samples[2], ".mod.seg", sep=""), header=T)
    cn3 = read.table(paste("sciclone_seg_files/",samples[3], ".mod.seg", sep=""), header=T)
    # clustering algorithm:
    sc = sciClone(vafs=list(v1,v2,v3), copyNumberCalls=list(cn1,cn2,cn3), sampleNames=samples, minimumDepth = 50, cnCallsAreLog2=TRUE)
  } else if ( length(samples) == 4) {
    v1 = read.table(paste("vaf_files/",samples[1], ".tsv", sep=""), header=T)
    v2 = read.table(paste("vaf_files/",samples[2], ".tsv", sep=""), header=T)
    v3 = read.table(paste("vaf_files/",samples[3], ".tsv", sep=""), header=T)
    v4 = read.table(paste("vaf_files/",samples[4], ".tsv", sep=""), header=T)
    cn1 = read.table(paste("sciclone_seg_files/",samples[1], ".mod.seg", sep=""), header=T)
    cn2 = read.table(paste("sciclone_seg_files/",samples[2], ".mod.seg", sep=""), header=T)
    cn3 = read.table(paste("sciclone_seg_files/",samples[3], ".mod.seg", sep=""), header=T)
    cn4 = read.table(paste("sciclone_seg_files/",samples[4], ".mod.seg", sep=""), header=T)
    # clustering algorithm:
    sc = sciClone(vafs=list(v1,v2,v3,v4), copyNumberCalls=list(cn1,cn2,cn3,cn4), sampleNames=samples, minimumDepth = 50, cnCallsAreLog2=TRUE)
  } else if ( length(samples) == 5) {
    v1 = read.table(paste("vaf_files/",samples[1], ".tsv", sep=""), header=T)
    v2 = read.table(paste("vaf_files/",samples[2], ".tsv", sep=""), header=T)
    v3 = read.table(paste("vaf_files/",samples[3], ".tsv", sep=""), header=T)
    v4 = read.table(paste("vaf_files/",samples[4], ".tsv", sep=""), header=T)
    v5 = read.table(paste("vaf_files/",samples[5], ".tsv", sep=""), header=T)
    #cn1 = read.table(paste("sciclone_seg_files/",samples[1], ".mod.seg", sep=""), header=T)
    #cn2 = read.table(paste("sciclone_seg_files/",samples[2], ".mod.seg", sep=""), header=T)
    #cn3 = read.table(paste("sciclone_seg_files/",samples[3], ".mod.seg", sep=""), header=T)
    #cn4 = read.table(paste("sciclone_seg_files/",samples[4], ".mod.seg", sep=""), header=T)
    #cn5 = read.table(paste("sciclone_seg_files/",samples[5], ".mod.seg", sep=""), header=T)
    # clustering algorithm:
    #sc = sciClone(vafs=list(v1,v2,v3,v4,v5), copyNumberCalls=list(cn1,cn2,cn3,cn4,cn5), sampleNames=samples, minimumDepth = 50, cnCallsAreLog2=TRUE)
    sc = sciClone(vafs=list(v1,v2,v3,v4,v5), sampleNames=samples, minimumDepth = 10)
  } else if ( length(samples) == 6) {
    # define files for input into sciclone
    v1 = read.table(paste("vaf_files/",samples[1], ".tsv", sep=""), header=T)
    v2 = read.table(paste("vaf_files/",samples[2], ".tsv", sep=""), header=T)
    v3 = read.table(paste("vaf_files/",samples[3], ".tsv", sep=""), header=T)
    v4 = read.table(paste("vaf_files/",samples[4], ".tsv", sep=""), header=T)
    v5 = read.table(paste("vaf_files/",samples[5], ".tsv", sep=""), header=T)
    v6 = read.table(paste("vaf_files/",samples[6], ".tsv", sep=""), header=T)
    cn1 = read.table(paste("sciclone_seg_files/",samples[1], ".mod.seg", sep=""), header=T)
    cn2 = read.table(paste("sciclone_seg_files/",samples[2], ".mod.seg", sep=""), header=T)
    cn3 = read.table(paste("sciclone_seg_files/",samples[3], ".mod.seg", sep=""), header=T)
    cn4 = read.table(paste("sciclone_seg_files/",samples[4], ".mod.seg", sep=""), header=T)
    cn5 = read.table(paste("sciclone_seg_files/",samples[5], ".mod.seg", sep=""), header=T)
    cn6 = read.table(paste("sciclone_seg_files/",samples[6], ".mod.seg", sep=""), header=T)
    # clustering algorithm:
    sc = sciClone(vafs=list(v1,v2,v3,v4,v5,v6), sampleNames=samples, minimumDepth = 50)
  }
  #create output
  writeClusterTable(sc, paste(tumor, "_clusters.tsv", sep=""))
  sc.plot1d(sc, paste(tumor, "_clusters.1d.pdf", sep=""))
  sc.plot2d(sc, paste(tumor, "_clusters.2d.pdf", sep=""))
}


 ############# PARSE THE OUTPUT TO DETERMINE TUMOR PURITY ################

# 11.1.19

#for ( tumor in tumor_names ) {

#  final_cluster_file <- read.table("~/Desktop/GBM/sciclone/BrMET010/clusters.tsv", header = TRUE)
#  founder_clone <- filter(final_cluster_file, cluster == 2)
#  sample1_purity <- median(founder_clone$BrMET010.1.vaf)*2
#  sample2_purity <- median(founder_clone$BrMET010.2.vaf)*2
#  sample3_purity <- median(founder_clone$BrMET010.3.vaf)*2

#}

