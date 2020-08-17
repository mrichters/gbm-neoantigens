# 10.24.19

# load packages
library(TPES)

# set wd
setwd("~/Desktop/GBM/tpes")

# example command
#TPES_purity(ID, SEGfile, SNVsReadCountsFile, ploidy, RMB = 0.47, maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)

# need 3 files to run tpes:
# 1. MAF - created from readcount files - sample, chr, start, end, ref.count, alt.count
# 2. ploidy - status per sample? how do i find this?
# 3. seg file - raw output from cnvkit

# read in folders with 3 required file types
seg_files <- unlist(lapply(list.files("seg_files"), function(x) paste("seg_files", x, sep = '/')))
maf_files <- unlist(lapply(list.files("maf_files"), function(x) paste("maf_files", x, sep = '/')))
ploidy_files <- unlist(lapply(list.files("ploidy_files"), function(x) paste("ploidy_files", x, sep = '/')))

# create a vector with all sample names
sample_names <- unlist(lapply(list.files("seg_files"), function(x) unlist(strsplit(x, split = ".seg"))[1]))

# loop over sample names and run tpes
results_df <- data.frame()
for ( sample in sample_names ) {
  
  # grep files from list vectors
  seg <- as.data.frame(read.table(seg_files[grep(paste(sample, ".seg", sep = ""), seg_files)], header=T))
  maf <- as.data.frame(read.table(maf_files[grep(paste(sample, ".maf.tsv", sep = ""), maf_files)], header=T))
  ploidy <- as.data.frame(read.table(ploidy_files[grep(paste(sample, ".ploidy.tsv", sep = ""), ploidy_files)], header=T))
  
  # run purity algorithm 
  purity_value <- TPES_purity(ID = sample, SEGfile = seg, SNVsReadCountsFile = maf, ploidy = ploidy, RMB = 0.47, maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)
  results_df <- rbind(results_df, purity_value)

}

rownames(results_df) <- NULL
write.table(results_df, "TPES_results.tsv", sep = '\t', quote = F, row.names = F)
