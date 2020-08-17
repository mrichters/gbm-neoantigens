################################################################################
################### Read in the data and load libraries ########################
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(viridis)
library(ggdendro)
library(gridExtra)
library(reshape2)
library(gtable)
library(data.table)
library(readxl)
library(dplyr)
library(mgsub)

setwd("~/Desktop/GBM/deconstructSigs/")

# read in data  
read_data <- function(batch) {
  final_df <- c()
  variant_path <- paste("~/Desktop/GBM/gbm_google_drive/", batch, "_annotated_variants.xlsx", sep="")
  samples <- excel_sheets(path = variant_path)
  if ( batch == "Batch3") {
    tumors <- unique(sapply(mgsub(samples, c("19-", "_"), c("", "-")), function(v) return(unlist(strsplit(v, '-'))[1])))
  } else {
    tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[2])))
  }
  for (tumor in tumors) {
    tumor_df <- c()
    for ( name in samples[grep(tumor, samples)] ) {
      myData <- read_excel(variant_path, sheet = name)
      myData <- myData[,c("CHROM", "POS", "REF", "ALT")]
      myData <- myData[,c("CHROM", "POS", "REF", "ALT")]
      colnames(myData) <- c("Chromosome", "Stop", "Reference", "Variant")
      tumor_df <- rbind(tumor_df, myData)
    }
    tumor_df <- tumor_df[order(tumor_df$Chromosome, tumor_df$Stop),]
    unique_tumor_df <- distinct(tumor_df)
    unique_tumor_df$sample <- tumor
    final_df <- rbind(final_df, unique_tumor_df)
  }
  return(final_df)
}

batch1 <- read_data("Batch1")
batch2 <- read_data("Batch2")
batch3 <- read_data("Batch3")

variants <- rbindlist(list(batch1, batch2, batch3))

################################################################################
#################### Filter the data to only snvs and format ###################

# Remove anything not a SNP and limit to just standard chromosomes
keep <- c("A", "C", "T", "G")
variants <- variants[Reference %chin% keep & Variant %chin% keep]
chromosomes <- paste0("chr", c(1:22, "X", "Y"))
variants <- variants[Chromosome %chin% chromosomes,]

# format the input
variants <- variants[,c("sample", "Chromosome", "Stop", "Reference", "Variant")]
colnames(variants) <- c("Sample", "chr", "pos", "ref", "alt")

# remove samples with less than 50 mutatations
sample_count <- as.data.frame(table(variants$Sample))
keep_sample <- as.character(filter(sample_count, Freq > 50)$Var1)
variants <- variants[Sample %chin% keep_sample]
variants$pos <- as.numeric(variants$pos)
################################################################################
############ convert to signature matrix #######################################
sigs.input <- mut.to.sigs.input(mut.ref = variants, sample.id = "Sample", chr = "chr",
                                pos = "pos", ref = "ref", alt = "alt",
                                bsg=BSgenome.Hsapiens.UCSC.hg38)


sig_data <- list()
a <- function(x){
  output.sigs <- whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.cosmic,
                                 sample.id = x, contexts.needed = TRUE, tri.counts.method = "exome2genome")
  pdf(file=paste0(x, "sigplot.s.pdf"), height=8, width=11)
  plotSignatures(output.sigs)
  dev.off()
  pdf(file=paste0(x, "pieplot.s.pdf"), height=8, width=8)
  makePie(output.sigs)
  dev.off()
  sig_data[[x]] <- output.sigs$weights
}
sig_data <- do.call(rbind, lapply(rownames(sigs.input), a))
sig_data$sample <- rownames(sig_data)
sig_data <- reshape2::melt(data=sig_data, id.vars=c("sample"))

# add back the unknown signatures
sig_data_unknown <- aggregate(data=sig_data, value ~ sample, sum)
sig_data_unknown$value <- 1-sig_data_unknown$value
sig_data_unknown$variable <- "Unknown"
sig_data <- rbind(sig_data, sig_data_unknown)

# write it out
write.table(sig_data, file="signature_weights.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# do clustering on the data 
sig_data_2 <- reshape2::dcast(sig_data, sample~variable, value.var="value")
rownames(sig_data_2) <- sig_data_2$sample
sig_data_2 <- sig_data_2[,!colnames(sig_data_2) %in% "sample"]
dist_matrix <- dist(sig_data_2)
fit <- hclust(dist_matrix)

# make the value column a proportion across the entire cohort
sig_data$value <- sig_data$value/length(unique(sig_data$sample))
################################################################################
############## create plot #####################################################

# Create a barchart
pdf(file="deconstructSig.barchart.pdf", height=8, width=12)
ggplot(sig_data, aes(x=variable, y=value, fill=sample)) + geom_bar(stat="identity") +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), plot.title = element_text(hjust = 0.5)) + ggtitle("GBM Signature Landscape") +
  xlab("Cosmic Signature") + ylab("Proportion")
dev.off()

################################################################################
########### create a heatmap ###################################################

# create a dendogram
dhc <- as.dendrogram(fit)
ddata <- dendro_data(dhc, type="triangle")
dendrogram <- ggplot(segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + 
  scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand=c(0.015, 0.015))+ theme_dendro()

# Create a heatmap
variable_order <- aggregate(data=sig_data, value ~ variable, sum)
variable_order <- as.character(variable_order[order(-variable_order$value),]$variable)
sig_data$variable <- factor(sig_data$variable, levels=variable_order)

sig_data$sample <- factor(sig_data$sample, levels=as.character(ddata[["labels"]]$label))

heatmap <- ggplot(sig_data, aes(x=variable, y=sample, fill=value)) + geom_tile() + scale_fill_viridis("Proportion", option="plasma", trans="sqrt") +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("Signature") + ylab("Sample") + 
  theme(legend.position="top", legend.text = element_text(angle=45, hjust=1))

# arrange plot data
dendrogram <- ggplotGrob(dendrogram)
heatmap <- ggplotGrob(heatmap)
dendrogram <- gtable_add_rows(dendrogram, heights=unit(1, "null"), pos=5)
dendrogram <- gtable_add_rows(dendrogram, heights=unit(1, "null"), pos=5)
maxheight <- grid::unit.pmax(dendrogram$heights, heatmap$heights)
dendrogram$heights <- as.list(maxheight)
heatmap$heights <- as.list(maxheight)


pdf(file="deconstructSig.heatmap.2.pdf", height=12, width=14)
grid.arrange(dendrogram, heatmap, widths=c(.2,.8))
dev.off()



