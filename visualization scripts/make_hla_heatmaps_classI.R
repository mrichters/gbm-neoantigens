# make heatmap with hla types

# read in txt file
my_data <- read.table('classI_hla-a.txt', sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE)

# brmet008
# convert data frame to a matrix
my_matrix <- as.matrix(my_data[c(1:9),c(2:11)])
# create row labels
rownames(my_matrix) <- c('Normal-DNA',paste('Tumor-DNA', 1:4),paste('Tumor-RNA', 1:4))
# save as jpeg
jpeg('HLA_typing_heatmaps/BrMET008B.jpg')
# Now we make our first heatmap 
Heatmap(my_matrix, row_names_side = "left", name = "HLA alleles", column_title = "BrMET008")
# this is necessary
dev.off()

# brmet009
my_matrix <- as.matrix(my_data[c(10:16),c(2:11)])
rownames(my_matrix) <- c('Normal-DNA',paste('Tumor-DNA', 1:3),paste('Tumor-RNA', 1:3))
jpeg('HLA_typing_heatmaps/BrMET009B.jpg')
Heatmap(my_matrix, row_names_side = "left", name = "HLA alleles", column_title = "BrMET009") 
dev.off()

#gbm 030
my_matrix <- as.matrix(my_data[c(17:21),c(2:11)])
rownames(my_matrix) <- c('Normal-DNA',paste('Tumor-DNA', 1:2),paste('Tumor-RNA', 1:2))
jpeg('HLA_typing_heatmaps/GBM030B.jpg')
Heatmap(my_matrix, row_names_side = "left", name = "HLA alleles", column_title = "GBM030")
dev.off()

#gbm 032
my_matrix <- as.matrix(my_data[c(22:29),c(2:11)])
rownames(my_matrix) <- c('Normal-DNA',paste('Tumor-DNA', 1:4),paste('Tumor-RNA', 1:3))
jpeg('HLA_typing_heatmaps/GBM032B.jpg')
Heatmap(my_matrix, row_names_side = "left", name = "HLA alleles", column_title = "GBM032")
dev.off()

#gbm 047
my_matrix <- as.matrix(my_data[c(30:35),c(2:11)])
rownames(my_matrix) <- c('Normal-DNA',paste('Tumor-DNA', 1:3),paste('Tumor-RNA', 2:3))
jpeg('HLA_typing_heatmaps/GBM047B.jpg')
Heatmap(my_matrix, row_names_side = "left", name = "HLA alleles", column_title = "GBM047")
dev.off()

#gbm 051
my_matrix <- as.matrix(my_data[c(36:42),c(2:11)])
rownames(my_matrix) <- c('Normal-DNA',paste('Tumor-DNA', 1:3),paste('Tumor-RNA', 1:3))
jpeg('HLA_typing_heatmaps/GBM051A.jpg')
Heatmap(my_matrix, row_names_side = "left", name = "HLA alleles", column_title = "GBM051")
dev.off()

#gbm 052
my_matrix <- as.matrix(my_data[c(43:48),c(2:11)])
rownames(my_matrix) <- c('Normal-DNA',paste('Tumor-DNA', 1:3),paste('Tumor-RNA', 2:3))
jpeg('HLA_typing_heatmaps/GBM052A.jpg')
Heatmap(my_matrix, row_names_side = "left", name = "HLA alleles", column_title = "GBM052")
dev.off()