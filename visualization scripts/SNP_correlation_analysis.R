library(stringr)
library(corrplot)
library("Hmisc")
library(readxl)
library(ggcorrplot)

file <- read.table("../correlation_p_values_new_cutoff_mr.tsv", sep="\t", header=FALSE)
file.dna <- read.table("../cor_matrix_DNA.tsv", sep="\t", header=FALSE)
file.rna <- read.table("../cor_matrix_RNA.tsv", sep="\t", header=FALSE)
file.dna.rna <- read.table("../cor_matrix_DNA_RNA.tsv", sep="\t", header=FALSE)
colnames(file.rna) <- c("row", "column", "cor", "p")

cor_data <- as.matrix(spread(file.dna[,-c(4)], key = column, value = cor))
rownames(cor_data) <- cor_data[,1]
cor_data <- cor_data[,-c(1)]

cor_data2<- as.matrix(spread(file.rna[,-c(4)], key = row, value = cor))
rownames(cor_data2) <- cor_data2[,1]
cor_data2 <- cor_data2[,-c(1)]

#then manually imported the matrices into excel and created a correlation matrix (basically filled out other half)
#now import excel back into Rstudio

DNA_matrix <- as.data.frame(read_excel("../SNP_correlation_sample_swap.xlsx", sheet ="DNA"))
RNA_matrix <- as.data.frame(read_excel("../SNP_correlation_sample_swap.xlsx", sheet ="RNA"))

rownames(DNA_matrix) <- DNA_matrix[,1]
DNA_matrix <- DNA_matrix[,-c(1)]
D <- data.matrix(DNA_matrix)
cor_D <- cor(D)
pcor_D <- cor_pmat(D)
ggcorrplot(cor_D, outline.col = "white")
ggcorrplot(cor_D, hc.order = TRUE, outline.col = "white")

rownames(RNA_matrix) <- RNA_matrix[,1]
RNA_matrix <- RNA_matrix[,-c(1)]
R <- data.matrix(RNA_matrix)
cor_R <- cor(R)
pcor_R <- cor_pmat(R)
ggcorrplot(cor_R, outline.col = "white")
ggcorrplot(cor_R, hc.order = TRUE, outline.col = "white", lab=TRUE)

# grep worked - and you dont need exact matches
#rbind(file[grep("TWCK.BrMET018.BrMET018_Tumor", file$row),], file[grep("TWCK.BrMET018.BrMET018_Tumor", file$column),])

