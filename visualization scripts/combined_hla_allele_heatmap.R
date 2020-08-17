# make one full prediction heatmap per hla allele

# read in txt file
my_data <- read.table('hla-a.txt', sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE)
# convert data frame to a matrix
my_matrix <- as.matrix(my_data[ ,c(17:21)])
# create row labels
#rownames(my_matrix) <- c('Normal-DNA',paste('Tumor-DNA', 1:4),paste('Tumor-RNA', 1:4))
#rownames(my_matrix) <- c(paste('BrMET008', 1:9), paste('BrMET009', 1:7), paste('GBM030', 1:5), paste('GBM032', 1:8), paste('GBM047', 1:6), paste('GBM051', 1:7), paste('GBM052', 1:6))
rownames(my_matrix) <- c('BrMET008 Normal',paste('BrMET008 tDNA', 1:4),paste('BrMET008 tRNA', 1:4), 'BrMET009 Normal',paste('BrMET009 tDNA', 1:3),paste('BrMET009 tRNA', 1:3), 'GBM030 Normal',paste('GBM030 tDNA', 1:2),paste('GBM030 tRNA', 1:2), 'GBM032 Normal',paste('GBM032 tDNA', 1:4),paste('GBM032 tRNA', 1:3), 'GBM047 Normal',paste('GBM047 tDNA', 1:3),paste('GBM047 tRNA', 2:3), 'GBM051 Normal',paste('GBM051 tDNA', 1:3),paste('GBM051 tRNA', 1:3),'GBM052 Normal',paste('GBM052 tDNA', 1:3),paste('GBM052 tRNA', 2:3))
# Now we make the heatmap 
Heatmap(my_matrix, row_names_side = "left", name = "HLA alleles", column_title = "DRB1", na_col = "grey")

