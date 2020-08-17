library(readxl)
library(tidyverse)

germline <- read_excel("~/Desktop/GBM/GBM065/germline_DNA/GBM065_germline_variants_annotated_vaf_mmrgenes.xlsx")

germline <- filter(germline, Consequence %in% c("missense_variant", "missense_variant&splice_region_variant", "inframe_deletion"))

germline <- as.data.frame(germline[,c(1:5,9,10,15,16)])

colnames(germline) <- c("Chr", "Pos", "Ref", "Alt", "GT", "Consequence", "Gene Name", "Protein Position", "Amino Acids")

mutation_col <- c()
for (row in 1:nrow(germline)) {
    name <- germline[row, "Gene Name"]
    position <- germline[row, "Protein Position"]
    aa <- unlist(strsplit(germline[row, "Amino Acids"], '/'))
    mutation <- paste(aa[1], position, aa[2], sep = "")
    mutation_col <- c(mutation_col, mutation)
}

germline["Mutation"] = mutation_col

germline[ ,c("Protein Position", "Amino Acids")]

germline <- as.data.frame(germline[,c(1:7,10)])

germline <- arrange(germline, `Gene Name`)

write_tsv(germline, "~/Desktop/gbm065_germline.tsv")
