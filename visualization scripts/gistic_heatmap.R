# 11.5.20

# import packages
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(colorRamps)
library(patchwork)

# functions
get_gbm_order <- function(vector) {
    re_samples <- unique(vector[grepl("Re", vector)])
    not_re_samples <- unique(vector[!grepl("Re", vector)])
    sample_order <- c(not_re_samples, re_samples)
    return(sample_order)
}

cnv_proportion <- function(cnv_id, df) {
    cnv_df <- filter(df, CNV == cnv_id)
    cnv_df[cnv_df == 2] <- 1
    cnv_samples <- filter(gather(cnv_df, key = "Sample", value = "CNV"), CNV == 1)$Sample
    cnv_tumors <- unique(sapply(cnv_samples, function(x) return(unlist(strsplit(x, '-'))[1])))
    all_tumors <- unique(sapply(colnames(cnv_df[,-1]), function(x) return(unlist(strsplit(x, '-'))[1])))
    proportion <- c(length(cnv_tumors), length(all_tumors))
    return(proportion)
}

format_df <- function(file_name) {
    # get CNV presence / absence data from gistic output - all_lesions.conf_75.txt
    cnv_df <- as.data.frame(read_tsv(file_name))
    cnv_df <- filter(cnv_df, ! grepl("CN values", `Unique Name`))
    cnv_df <- cnv_df[,c(1,2,11:ncol(cnv_df)-1)]
    # change unique name to Amp or Del
    cnv_df$Type <- sapply(cnv_df$`Unique Name`, function(x) return(substr(x, 1, 3)))
    # gather CNV df to get unique CNVs
    cnv_df <- gather(cnv_df, key = Samples, value = `CNV Threshold`, -c(Type, `Unique Name`, Descriptor))
    # get sample order, convert to factor and order df
    cnv_df$Samples_format <- sapply(cnv_df$Samples, function(x) return(gsub('-', '  ', x)))
    gbm_order <- get_gbm_order(cnv_df$Samples_format)
    cnv_df$Samples_format <- factor(cnv_df$Samples_format, levels = gbm_order)
    cnv_df <- cnv_df[order(cnv_df$Samples_format), ]
    # amp/del column, factor and order
    cnv_df$Descriptor <- factor(cnv_df$Descriptor, levels = unique(str_sort(cnv_df$Descriptor, numeric=T, decreasing=T)))
    cnv_df <- cnv_df[order(cnv_df$Descriptor), ]
    # format type column
    cnv_df$Type <- factor(cnv_df$Type, levels = c("Amp", "Del"))
    cnv_df <- cnv_df[order(cnv_df$Type), ]
    # make deletions negative for the heatmap
    for (n in 1:nrow(cnv_df)) {
        if (grepl("Del", cnv_df$Type[n])) {
            cnv_df$`CNV Threshold`[n] <- cnv_df$`CNV Threshold`[n] * -1
        }
    }
    return(cnv_df)
}

plot_heatmap <- function(df, line_vector, title, labels) {
    ggplot(df, aes(x=Samples_format, y=Descriptor, fill=`CNV Threshold`)) + geom_tile() + geom_vline(xintercept=line_vector) + theme_bw() + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.title.x.bottom=element_blank(), axis.ticks.y=element_blank(), legend.position = "right", text = element_text(size=8), legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) + scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limit = c(-2,2)) + ggtitle(title) + scale_x_discrete(labels = labels)
}

lung_df <- format_df("~/Desktop/GBM/cnvkit/07.2020/brmet_lung/all_lesions.conf_75.txt") 
breast_df <- format_df("~/Desktop/GBM/cnvkit/07.2020/brmet_breast/all_lesions.conf_75.txt")
gbm_df <- format_df("~/Desktop/GBM/cnvkit/07.2020/gbm_0.1/all_lesions.conf_75.txt")

lung_plot <- plot_heatmap(lung_df, c(4.5,7.5,10.5,13.5,16.5), "BrMET NSCLC", c("BrMET008  1" = "1", "BrMET008  3" = "3", "BrMET008  4" = "4", "BrMET009  1" = "1", "BrMET009  3" = "3", "BrMET010  1" = "1", "BrMET010  3" = "3", "BrMET019  1" = "1", "BrMET019  3" = "3", "BrMET025  1" = "1", "BrMET025  3" = "3"))
breast_plot <- plot_heatmap(breast_df, c(3.5,6.5,9.5), "BrMET BRCA", c("BrMET018  1" = "1", "BrMET018  3" = "3", "BrMET023  1" = "1", "BrMET023  3" = "3", "BrMET024  1" = "1", "BrMET024  3" = "3", "BrMET027  1" = "1", "BrMET027  3" = "3"))
gbm_plot <- plot_heatmap(gbm_df, c(2.5,6.5,9.5, 12.5, 15.5, 18.5, 21.5, 24.5, 27.5, 30.5, 33.5, 36.5, 39.5, 42.5, 45.5, 48.5, 52.5, 55.5), "GBM", c("GBM030  2" = "2","GBM032  1" = "1", "GBM032  3" = "3", "GBM032  4" = "4", "GBM051  1" = "1", "GBM051  3" = "3", "GBM052  1" = "1", "GBM052  3" = "3", "GBM055  1" = "1", "GBM055  3" = "3", "GBM056  1" = "1", "GBM056  3" = "3", "GBM059  1" = "1", "GBM059  3" = "3", "GBM062  1" = "1", "GBM062  3" = "3", "GBM063  1" = "1", "GBM063  3" = "3", "GBM064  1" = "1", "GBM064  3" = "3", "GBM069  1" = "1", "GBM069  3" = "3", "GBM070  1" = "1", "GBM070  3" = "3", "GBM074  1" = "1", "GBM074  3" = "3", "GBM079  1" = "1", "GBM079  3" = "3", "GBM083  1" = "1", "GBM083  3" = "3", "GBM018.Re  1" = "1", "GBM018.Re  3" = "3", "GBM031.Re  1" = "1", "GBM031.Re  3" = "3", "GBM031.Re  4" = "4", "GBM047.Re  1" = "1", "GBM047.Re  3" = "3", "GBM065.Re  1" = "1", "GBM065.Re  4" = "4", "GBM065.Re  5" = "5"))

cnv_heatmap <- gbm_plot + lung_plot + breast_plot + plot_layout(ncol = 3, widths = c(4,1.1,0.8), guides = 'collect')

## get proportion of samples with CNVs ##
# gbm
EGFR <- cnv_proportion("Amp:7p11.2")
CDKN2A <- cnv_proportion("Del:9p21.3")
# breast
ERBB2 <- cnv_proportion("Amp:17q12")
# lung
KRAS <- cnv_proportion("Amp:12p12.1")

ggplot(df.g, aes(x=Samples, y=CNV, fill=`CNV Threshold`)) + geom_tile() + geom_vline(xintercept=c(4.5,7.5,10.5,13.5,16.5,17.5)) + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.title.x.bottom=element_blank(), axis.ticks.y=element_blank(), legend.position = "right", text = element_text(size=13), plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio=0.75) + scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limit = c(-2,2)) + ggtitle("BrMET NSCLC")
# coolwarm()
# brewer.pal(RdBu)

ggplot(df.g, aes(x=Status, y=CNV, fill=Type)) + geom_tile() + coord_fixed(ratio=1) + theme(axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", axis.text.x=element_text(angle=65, hjust=1), text = element_text(size=14))
 
# CNV heterogeneity estimates
# to create df for variant clonality plots (clonal, subclonal shared, private)
clonality_plot <- function(file_name) {
    df <- as.data.frame(read_tsv(file_name))
    df <- filter(df, !grepl("CN", `Unique Name`))
    tumors <- unique(sapply(names(df)[grepl("BrMET|GBM", names(df))], function(x) unlist(strsplit(x, '-'))[1]))
    clonal_df <- c()
    for (t in tumors) {
        df.spread <- df[,grepl(t, names(df))]
        df.spread[df.spread == 2] = 1
        df.spread$row_sums <- rowSums(df.spread)
        max_clonal_num <- max(df.spread$row_sums)
        df.spread <- filter(df.spread, row_sums != 0)
        for ( line in 1:nrow(df.spread) ) {
            count <- df.spread[line, "row_sums"]
            if ( max_clonal_num == 2 ) {
                if ( count == 1 ) {
                    df.spread[line, "Coverage"] <- "Subclonal Private"
                } else if ( count == 2 ) {
                    df.spread[line, "Coverage"] <- "Clonal"
                }
            } else {
                if ( count == 1 ) {
                    df.spread[line, "Coverage"] <- "Subclonal Private"
                } else if ( count %in% 2:(max_clonal_num-1)) {
                    df.spread[line, "Coverage"] <- "Subclonal Shared"
                } else if ( count == max_clonal_num ) {
                    df.spread[line, "Coverage"] <- "Clonal"
                }
            }
        }
        final_spread <- data.frame(Tumor = t, as.data.frame(table(df.spread$Coverage)))
        colnames(final_spread) <- c("Tumor", "Variant Type", "Variant Count")
        clonal_df <- rbind(clonal_df, final_spread)
    }
    
    clonal_df$`Variant Type` <- factor(clonal_df$`Variant Type`, levels = c("Subclonal Private", "Subclonal Shared", "Clonal"))
    clonal_df$Tumor <- factor(clonal_df$Tumor, levels = get_gbm_order(clonal_df$Tumor))
    return(clonal_df)
}  
    # create clonal
plot_cnv_clonal <- function(df, title) {
    ggplot(df, aes(x=Tumor, y=`Variant Count`, fill=`Variant Type`)) + geom_bar(stat="identity") + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=8), legend.position = "right", plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) + ylab("CNV Count") + scale_fill_manual(values = c("forestgreen", "royalblue3", "firebrick3")) + ggtitle(title) + coord_cartesian(expand = F)
}

lung_df2 <- clonality_plot("~/Desktop/GBM/cnvkit/07.2020/brmet_lung/all_lesions.conf_75.txt")
breast_df2 <- clonality_plot("~/Desktop/GBM/cnvkit/07.2020/brmet_breast/all_lesions.conf_75.txt")
gbm_df2 <- clonality_plot("~/Desktop/GBM/cnvkit/07.2020/gbm_0.1/all_lesions.conf_75.txt")

l <- plot_cnv_clonal(lung_df2, "BrMET NSCLC")
b <- plot_cnv_clonal(breast_df2, "BrMET BRCA")
g <- plot_cnv_clonal(gbm_df2, "GBM")

cnv_clonality <- g + l + b + plot_layout(widths = c(4,1.2,0.8), guides = 'collect')
