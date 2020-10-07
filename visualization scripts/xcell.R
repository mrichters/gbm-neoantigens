# 02.26.20
# analyze xCell output
# Usage: Rscript xcell.R

# load packages
library(tidyverse)
library(readxl)
library(textshape)
library(ggsci)
library(ComplexHeatmap)
library(BBmisc)
library(wesanderson)
library(grid)
library(gridExtra)
library(stringi)
library(mgsub)

#functions
create_initial_df <- function(input_file, transpose=TRUE) {
    if ( transpose == T ) {
        # create rownames from samples
        df <- column_to_rownames(input_file, loc = 1)
        # transpose the data
        df_t <- as.data.frame(t(df))
        # undo rownames and make Samples a column again
        df_t_c <- rownames_to_column(df_t, var = "Samples")
        # change samples from "BrMET008.1" "BrMET008-1"
        df_t_c$Samples <- sapply(df_t_c$Samples, function(x) return(stri_replace_last_fixed(x, '.', '-')))
        # rename df_t_c to results
        results <- df_t_c
    } else {
        # read in file
        results <- read_tsv(input_file, header=T)
    }
    # add Tumor column by removing sample number from each line
    #results$Tumor <- sapply(results$Samples, function(x) return(substr(x, 1, nchar(x)-2)))
    # add Tumor Type column
    #results$`Tumor Type` <- sapply(results$Tumor, function(x) return(str_extract(x, "[a-zA-Z]+" )))
    return(results)
}

get_gbm_order <- function(vector) {
    re_samples <- unique(vector[grepl("Re", vector)])
    not_re_samples <- unique(vector[!grepl("Re", vector)])
    sample_order <- c(not_re_samples, re_samples)
    tumor_order <- unique(sapply(sample_order, function(x) return(substr(x, 1, nchar(x)-2))))
    order_list <- list(sample_order, tumor_order)
    return(sample_order)
}

convert_to_tumor <- function(df) {
    output_df <- data.frame()
    for ( tumor in unique(df$Tumor) ) {
        df_subset <- filter(df, Tumor == tumor)
        df_subset_spread <- spread(df_subset, Samples, Score)
        score_means <- rowMeans(df_subset_spread[4:ncol(df_subset_spread)])
        new_line <- data.frame(df_subset_spread[1:3], "Score" = score_means)
        output_df <- rbind(output_df, new_line)
    }
    return(output_df)
}

# factor_and_order <- function(df, col, order) {
#     df[col] <- factor(df[col], levels = order)
#     df.ordered <- df[order(df[col]), ]
#     return(df.ordered)
# }

format_names <- function(name) {
    if ( grepl("Re", name) ) {
        if ( grepl("Primary", name) ) {
            new_name <- "Re.GBM065-primary"
            #tumor <- "GBM065"
        } else {
            new_name <-  paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
            #tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
        }
    } else if ( grepl("^19-0", name) ) {
        new_name <- mgsub(name, c("19-\\d*_", "Tumor", "_"), c("", "", "-"))
        #tumor <- unlist(strsplit(new_name, '-'))[1]
    } else if ( grepl("Re", name) == FALSE ) {
        new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
        #tumor <- unlist(str_split(name, '-'))[2]
    }
    return(new_name)
}

# setwd
setwd("~/Desktop/GBM/immune_infiltration/xcell/")

# load files
results <- read_tsv("xCell_kallisto_abundance_matrix_strandfix_09-21-20.txt")

# create initial df
xcell_df <- create_initial_df(results)

# get correct gbm sample/tumor order
#gbm_order <- get_gbm_order(xcell_df)
# factor and order
#xcell_df$Samples <- factor(xcell_df$Samples, levels = gbm_order)
#xcell_df <- xcell_df[order(xcell_df$Samples), ]
#xcell_df$Samples <- as.character(xcell_df$Samples)

#add tumor and tumor type
xcell_df$Tumor <- sapply(xcell_df$Samples, function(x) return(substr(x, 1, nchar(x)-2)))
xcell_df$Tumor_Type <- sapply(xcell_df$Tumor, function(x) return(unlist(stri_extract_all_regex(x, "[a-zA-Z.]+"))))

#cell_levels <- names(sort(apply(df_celltypes[,2:ncol(df_celltypes)], 2, mean),decreasing = TRUE))

#xcell_df.tidy <- gather(xcell_df, key = Cell_Type, value = xCell_Score, -Samples, -Tumor, -Tumor_Type)

column_totals <- apply(select_if(xcell_df, is.numeric), 2, sum)
nonzero_cols <- names(Filter(function(x) x > 0, column_totals))
xcell_df_subset <- xcell_df[,c("Samples", "Tumor", "Tumor_Type", nonzero_cols)]


xcell_norm <- normalize(xcell_df_subset, method = "standardize", range = c(0,1), margin=1)

#xcell_gather <- gather(xcell_norm, Cell_Type, Score, -Samples, -Tumor, -Tumor_Type)
#xcell_order <- factor_and_order(xcell_gather, xcell_gather$Samples, sample_order)

# create matrix for heatmap (Complex Heatmap)
xcell_transpose <- as.data.frame(t(column_to_rownames(xcell_df_subset[,-c(2:3)], loc=1)))
xcell_norm <- normalize(xcell_transpose, method = "range", range = c(0,1), margin=1)
xcell_matrix <- as.matrix(xcell_norm)
# create Heatmap!
Heatmap(xcell_matrix, name = "Norm. xCell Score", row_order = rownames(xcell_transpose), column_order = get_gbm_order(names(xcell_transpose)),width = unit(30, "cm"), height = unit(24, "cm"))

#### T CELLS ONLY SAMPLES ####
# vector of columns I want to include
Tcell_cols <- c("Samples", "Tumor", colnames(xcell_df)[grep("CD", colnames(xcell_df))], "Tumor Type")
# subset original df
Tcell_results <- xcell_df[Tcell_cols]
# these values are sums of all CD4 and CD8 cell types in xcell
Tcell_results$`CD4+ Cells` <- rowSums(Tcell_results[,3:7])
Tcell_results$`CD8+ Cells` <- rowSums(Tcell_results[,8:11])
# subset the data to new columns
Tcell_results_subset <- Tcell_results[c("Samples", "CD8+ Cells", "CD4+ Cells", "Tumor", "Tumor Type")]
# normalize 2 cell type values for placement in heatmap
Tcell_norm <- normalize(Tcell_results_subset, margin=2, range = c(0,1), method="range")
# gather cell types
Tcell_norm_gather <- gather(Tcell_norm, Cell_Type, Score, `CD4+ Cells`, `CD8+ Cells`)
# factor and order - ready to plot
Tcell_norm_gather$Samples <- factor(Tcell_norm_gather$Samples, levels = sample_order)
Tcell_norm_gather.ordered <- Tcell_norm_gather[order(Tcell_norm_gather$Samples), ]

#### T CELLS ONLY TUMOR ####
# now do same thing but with tumor average of scores instead of samples
Tcell_tumor <- convert_to_tumor(Tcell_norm_gather)
# factor and order - ready to plot
Tcell_tumor.ordered <- factor_and_order(Tcell_tumor, Tcell_tumor$Tumor, tumor_order)

#### DC CELLS ONLY SAMPLES ####
# plot all DCs in heatmap
# subset df for all types of T cells
DC_cols <- c("Samples", "Tumor", "Tumor Type", colnames(xcell_df)[grep("DC", colnames(xcell_df))])
DC_results <- xcell_df[DC_cols]
# these values are sums of all DC cell types in xcell
DC_results$`DC Cells` <- rowSums(DC_results[,4:8])
# normalize
DC_norm <- normalize(DC_results, margin=2, range = c(0,1), method="range")
# gather cell types
DC_norm_gather <- gather(DC_norm, Cell_Type, Score, -Samples, -Tumor, -`Tumor Type`)
# factor and order - ready to plot
DC_norm_gather$Samples <- factor(DC_results_norm_gather$Samples, levels = sample_order)
DC_norm_gather.ordered <- DC_results_norm_gather[order(DC_results_norm_gather$Samples), ]
# make cell types a factor with specified levels
dc_levels <- rev(c("DC", "aDC", "cDC", "iDC", "pDC", "DC Cells"))
DC_norm_gather.ordered$Cell_Type <- factor(DC_norm_gather.ordered$Cell_Type, levels = dc_levels)

#### DC CELLS ONLY TUMOR ####
DC_tumor <- convert_to_tumor(DC_norm_gather)
# factor and order - ready to plot
DC_tumor$Tumor <- factor(DC_tumor$Tumor, levels = tumor_order)
DC_tumor.ordered <- DC_tumor[order(DC_tumor$Tumor), ]
# make cell types a factor with specified levels
dc_levels <- rev(c("DC", "aDC", "cDC", "iDC", "pDC", "DC Cells"))
DC_tumor.ordered$Cell_Type <- factor(DC_tumor.ordered$Cell_Type, levels = dc_levels)

# palettes
my_pal <- c(pal_npg("nrc", alpha = 0.8)(10)[3:4], pal_simpsons("springfield", alpha = 1)(10))
pal = wes_palette("Zissou1", 100, type = "continuous")

#### BAR / BOX / VIOLIN ####
# bar
ggplot(xcell_df, aes(x=Samples, y=ImmuneScore)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# stacked bar
ggplot(Tcell_norm_gather.ordered, aes(x=Samples, y=Score, fill=Cell_Type)) + geom_bar(stat="identity", width=0.95)+ theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") + guides(fill=guide_legend(title="Cell Type")) #+ scale_fill_manual(values = my_pal) 

# split violin plot
ggplot(Tcell_norm_gather.ordered, aes(x=`Tumor Type`, y=Score, fill=Cell_Type)) + geom_violin(scale="width", trim=F) + scale_fill_manual(values = my_pal) + theme_bw() + theme(legend.position = "bottom", legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_rect(color = "black", fill=NA, size=1), text = element_text(size=12)) + ylab("Normalized xCell score") #+ ylim(0,0.3)

# box & whisker / violin plot
xcell_df$Tumor <- factor(xcell_df$Tumor, levels = tumor_order)
xcell_df.ordered <- xcell_df[order(xcell_df$Tumor), ]
ggplot(xcell_df.ordered, aes(x=Tumor, y=ImmuneScore)) + geom_boxplot(color="black", fill=NA) + geom_point(color="grey40",size=2) + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), legend.position = "none" , axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_rect(color = "black", fill=NA, size=1), text = element_text(size=16))

#### HEATMAPS ####

# heatmap with ggtile - using coord_fixed to fix dimensions!
# BrMET
p1 <- ggplot(filter(Tcell_norm_gather.ordered, `Tumor Type` == "BrMET"), aes(x=Samples, y=Cell_Type, fill=Score)) + geom_tile() + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.title.x.bottom=element_blank(), legend.position = "none", axis.ticks.y=element_blank(), text = element_text(size=14), plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio=1) + scale_fill_gradientn(colours = pal) + ggtitle("BrMET")

# GBM
p2 <- ggplot(xcell_gather, aes(x=Samples, y=Cell_Type, fill=Score)) + geom_tile() + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text = element_text(size=14), legend.position = "bottom", legend.key.width=unit(1.3,"cm"), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio=1) + scale_fill_gradientn(colours = pal) #+ ggtitle("Primary/Recurrent GBM")

# combine the two plots
grid_plot <- grid.arrange(arrangeGrob(p1, p2, heights=c(1.15,2)))
# save the combined plot
pdf("~/Desktop/GBM/immune_infiltration/xcell/gbm_xcell/heatmap_tcells_sample.pdf", width = 12, height = 5)
plot(grid_plot)
dev.off()

# plot with tumor 
ggplot(Tcell_norm_tumor.ordered, aes(x=Tumor,  y=Cell_Type, fill=Score)) + geom_tile() + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text = element_text(size=14), legend.position = "bottom", legend.key.width=unit(1,"cm"), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio=1) + scale_fill_gradientn(colours = pal)

# plot every cell type in large heatmap
ggplot(results_t_c.gather.ordered, aes(x=Samples,  y=Cell_Type, fill=Score)) + geom_tile() + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text = element_text(size=14), legend.position = "bottom", legend.key.width=unit(1,"cm"), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio=1) + scale_fill_gradientn(colours = pal)

# plot all DCs in heatmap
ggplot(xcell_order, aes(x=Samples,  y=Cell_Type, fill=Score)) + geom_tile() + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text = element_text(size=12), legend.position = "bottom", legend.key.width=unit(1,"cm"), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio=1) + scale_fill_gradientn(colours = pal) #, limits=c(0, 1), breaks=seq(0,1,by=0.25))
DC_norm_gather.ordered

# plot with tumor 
ggplot(DC_tumor.ordered, aes(x=Tumor,  y=Cell_Type, fill=Score)) + geom_tile() + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text = element_text(size=14), legend.position = "bottom", legend.key.width=unit(1,"cm"), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio=1) + scale_fill_gradientn(colours = pal, limits=c(0, 1), breaks=seq(0,1,by=0.25))
