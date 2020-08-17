# 10.24.19

# ABSOLUTE package downloaded as tar.gz file
#install.packages("numDeriv")
#install.packages("~/Documents/Rpackages/ABSOLUTE_1.0.6.tar.gz", repos = NULL, type ='source')

# load packages
library(numDeriv)
library(ABSOLUTE)

# load input data
home_dir <- "~/Desktop/GBM/ABSOLUTE"
input_data_folder <- "~/Desktop/GBM/ABSOLUTE/seg_files_absolute"
seg_files <- list.files(input_data_folder)
file <- seg_files[1]
# run abolute
for ( file in seg_files ) {
  # example command
  #RunAbsolute(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, platform, sample.name, results.dir, max.as.seg.count, max.non.clonal, max.neg.genome, copy_num_type, maf.fn=NULL, min.mut.af=NULL, output.fn.base=NULL, verbose=FALSE)
  sample_name <- unlist(strsplit(file, '.seg'))
  RunAbsolute(seg.dat.fn=paste(input_data_folder, file, sep='/'), sigma.p=0, max.sigma.h=0.02, min.ploidy=0.95, max.ploidy=10, primary.disease="Lung adenocarcinoma", platform = "Illumina_WES", sample.name = sample_name, results.dir = "~/Desktop/GBM/ABSOLUTE", max.as.seg.count=1500, max.non.clonal=0, max.neg.genome=0, copy_num_type = "total")

CreateReviewObject(obj.name = sample_name, absolute.files = "~/Desktop/GBM/ABSOLUTE/H_TC-BrMET008-008-Tumor-1-DNA.ABSOLUTE.Rdata", indv.results.dir = "~/Desktop/GBM/ABSOLUTE/test", copy_num_type = "total")
    
}