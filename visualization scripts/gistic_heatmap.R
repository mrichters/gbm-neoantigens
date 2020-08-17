library(tidyverse)
library(wesanderson)
library(stringr)
library(colorRamps)
library(scico)

get_gbm_order <- function(df) {
    recurrent_samples <- sort(as.character(unique(df[grep("Re", df$Samples), ]$Samples)))
    not_recurrent_samples <- sort(as.character(unique(filter(df, !Samples %in% recurrent_samples)$Samples)))
    sample_order <- c(not_recurrent_samples, recurrent_samples)
    #order_list <- list(sample_order, tumor_order, tumor_type)
    return(sample_order)
}

df <- as.data.frame(read_tsv("~/Desktop/GBM/cnvkit/07.2020/brmet_lung/all_lesions.conf_75.txt"))
df <- filter(df, ! grepl("CN values", `Unique Name`))
df <- df[,c(1,2,10:27)]


df$`Unique Name` <- sapply(df$`Unique Name`, function(x) return(substr(x, 1, 3)))
df.c <- unite(df, CNV, `Unique Name`, Descriptor, sep=":")
df.g <- gather(df.c, key = Samples, value = `CNV Threshold`, -CNV)
gbm_order <- get_gbm_order(df.g)
df.g$Samples <- factor(df.g$Samples, levels = gbm_order)
df.g <- df.g[order(df.g$Samples), ]
# amp/del column
df.g$CNV <- factor(df.g$CNV, levels = unique(str_sort(df.g$CNV, numeric=T, decreasing=T)))
df.g <- df.g[order(df.g$CNV), ]
df.g$Type <- sapply(df.g$CNV, function(x) return(substr(x, 1, 3)))
#df.g$Status <- "BrMET008-1"

for (n in 1:nrow(df.g)) {
    if (grepl("Del", df.g$CNV[n])) {
        df.g$`CNV Threshold`[n] <- df.g$`CNV Threshold`[n] * -1
    }
}

pal = wes_palette("Zissou1", 5, type = "continuous")

ggplot(df.g, aes(x=Samples, y=CNV, fill=`CNV Threshold`)) + geom_tile() + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.title.x.bottom=element_blank(), axis.ticks.y=element_blank(), legend.position = "right", text = element_text(size=13), plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio=0.75) + scale_fill_scico(palette = "vikO", limit = c(-2,2)) + ggtitle("BrMET NSCLC")
#+ scale_fill_gradientn(colours = pal,limits=c(-2, 2), breaks=seq(-2,2,by=1))
ggplot(df.g, aes(x=Status, y=CNV, fill=Type)) + geom_tile() + coord_fixed(ratio=1) + theme(axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", axis.text.x=element_text(angle=65, hjust=1), text = element_text(size=14))
 
