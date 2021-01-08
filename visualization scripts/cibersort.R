# import packages
library(readxl)
library(tidyverse)
library(stringr)
library(ggsci)
library(viridis)
library(RColorBrewer)

# set wd
setwd("~/Desktop/GBM/immune_infiltration/cibersort/")
input_file <- "CIBERSORT.Output.xlsx"

# format_names <- function(name) {
#     if ( grepl("Re", name) ) {
#         if ( grepl("Primary", name) ) {
#             new_name <- "Re.GBM065-primary"
#             #tumor <- "GBM065"
#         } else {
#         new_name <-  paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
#         #tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
#         }
#     } else if ( grepl("^19-0", name) ) {
#         new_name <- mgsub(name, c("19-\\d*_", "Tumor", "_"), c("", "", "-"))
#         #tumor <- unlist(strsplit(new_name, '-'))[1]
#     } else if ( grepl("Re", name) == FALSE ) {
#         new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
#         #tumor <- unlist(str_split(name, '-'))[2]
#     }
#     return(new_name)
# }

get_gbm_order <- function(vector) {
    re_samples <- unique(vector[grepl("Re", vector)])
    not_re_samples <- unique(vector[!grepl("Re", vector)])
    sample_order <- c(not_re_samples, re_samples)
    tumor_order <- unique(sapply(sample_order, function(x) return(substr(x, 1, nchar(x)-2))))
    order_list <- list(sample_order, tumor_order)
    return(sample_order)
}

# open cibersort results
df <- as.data.frame(read_excel(input_file))
# format sample names
#df$`Input Sample` <- sapply(df$`Input Sample`, format_names)
# subset df to reflect cell types
df_celltypes <- subset(df, select = -c(`P-value`, `Pearson Correlation`, RMSE))
# levels - ordering by mean cell fraction
cell_levels <- names(sort(apply(df_celltypes[,2:ncol(df_celltypes)], 2, mean),decreasing = TRUE))
# reorder df by levels
df_celltypes.reorder <- df_celltypes[c("Input Sample", cell_levels)]
# tidy the df - gather
df_celltypes.tidy <- gather(df_celltypes.reorder, key = `Cell Type`, value = Proportion, -`Input Sample`)
# change the order of the legend based on proportion
df_celltypes.tidy$`Cell Type` <- factor(df_celltypes.tidy$`Cell Type`, levels = unique(df_celltypes.tidy$`Cell Type`))
# create levels and change order of df
df_celltypes.tidy$`Input Sample` <- factor(df_celltypes.tidy$`Input Sample`, levels = get_gbm_order(unique(df_celltypes.tidy$`Input Sample`)))
df_celltypes.tidy <- df_celltypes.tidy[order(df_celltypes.tidy$`Input Sample`), ]
###

# change column names
colnames(df_celltypes.tidy) <- c("Samples", "Cell Type", "Proportion")

# create a color palette 
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
color_palette <- heat.colors(8)
my_pal <- brewer.pal(8, "Set1") #pal_npg("nrc")(8) #, pal_npg("nrc")(10), pal_simpsons("springfield")(5))

# make a stacked bar plot showing proportions of cell types per sample
# really messy because 22 cell types
ggplot(df_celltypes.tidy, aes(x=Samples, y=Proportion, fill=`Mod_cell_type`)) + geom_bar(stat="identity", width=0.95)+ theme(legend.position="top", axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=6)) + guides(fill=guide_legend(title="Cell Type")) + scale_fill_manual(values = viridis(8)) + coord_cartesian(expand = F)

#### COMBINE CELL TYPES FOR CLEARER PLOT ####
# make these into lists
# combined cell types
CD8_cells <- c("CD8+ T cells", "T cells CD8")
CD4_cells <- c("CD4+ T cells", "T cells CD4 memory activated", "T cells CD4 naive", "T cells CD4 memory resting", "T cells follicular helper")
B_cells <- c("B cells", "B cells naive", "B cells memory")
NK_cells <- c("NK cells", "NK cells resting", "NK cells activated")
DC_cells <- c("DC cells", "Dendritic cells activated", "Dendritic cells resting")
Macrophages <- c("Macrophages", "Macrophages M0", "Macrophages M1", "Macrophages M2")
Monocytes <- c("Monocytes")
Other <- c("Other", "Neutrophils", "Eosinophils", "Mast cells activated", "Mast cells resting", "Plasma cells", "T cells gamma delta", "T cells regulatory (Tregs)")
# combined cell types vector 
cell_types <- list(CD8_cells = CD8_cells, CD4_cells = CD4_cells, B_cells = B_cells, NK_cells = NK_cells, DC_cells = DC_cells, Macrophages = Macrophages, Monocytes = Monocytes, Other = Other)

CellTYPE <- list(CD8_cells, CD4_cells, B_cells, NK_cells, DC_cells, Macrophages, Monocytes, Other)

df_celltypes.tidy$`Mod Cell Type` <- "type"
for (i in 1:nrow(df_celltypes.tidy)) {
    for (j in Cell.types) {
        print(df_celltypes.tidy[1,"Cell Type"])
        if (filter(i, `Cell Type` %in% j)) {
            df_celltypes.tidy[i, "Mod. Cell Type"] <- j
        }
    }
}

# access cell type per line
find_types <- function(v) {
    for (j in 1:length(cell_types)) {
        #print(cell_types[[j]][1])
        if (v %in% cell_types[[j]]) {
            return(cell_types[[j]][1])
        }
    }
}
test <- sapply(df_celltypes.tidy$`Cell Type`, find_types)
# see if cell type is in loop of general cell type vectors
# if yes replace row entry with general cell type

find_types("Neutrophils")
# create total_df for plotting combined cell types
#tumors <- unique(sapply(df_celltypes.tidy$Samples, function(x) return(substr(x, 1, nchar(x)-2))))
samples <- unique(df_celltypes.tidy$`Input Sample`)
total_df <- data.frame()
for ( item in samples ) {
    sample_df <- data.frame()
    for ( num in 1:length(cell_types) ) {
        cell_types_df <- filter(df_celltypes.tidy, Samples == item)
        new_line <- data.frame(item, cell_types[[num]][1], sum(filter(cell_types_df, Samples == item & `Cell Type` %in% unlist(cell_types[[num]]))$Proportion))
        colnames(new_line) <- c("Samples", "Cell Type", "Proportion")
        sample_df <- rbind(sample_df, new_line)
    }
    print(sum(sample_df$Proportion))
    total_df <- rbind(total_df, sample_df)
}

#total_df.s <- spread(total_df, key=`Cell Type`, value=Proportion)
# levels - ordering by mean cell fraction
#cell_levels <- names(sort(apply(total_df.s[,2:ncol(total_df.s)], 2, mean),decreasing = TRUE))
# change the order of the legend based on preference
total_df$`Cell Type` <- factor(total_df$`Cell Type`, levels = c("CD8+ T cells", "CD4+ T cells", "Monocytes", "Macrophages", "B cells", "DC cells","NK cells", "Other"))
# change the order of the samples
total_df$Samples <- factor(total_df$Samples, levels = get_gbm_order(unique(total_df$Samples)))

# adjusted color palette
my_pal <- c(pal_nejm("default", alpha = 1)(7), "grey60")
# make a stacked bar plot showing proportions of cell types per sample
ggplot(total_df, aes(x=Samples, y=Proportion, fill=`Cell Type`)) + geom_bar(stat="identity", width=0.95)+ theme(legend.position="top", axis.text.x = element_text(angle = 65, hjust = 1), text = element_text(size=12)) + guides(fill=guide_legend(title="Cell Type")) + scale_fill_manual(values = my_pal) 

# to create Heatmap
total_df.s <- column_to_rownames(spread(total_df, key=`Cell Type`, value=Proportion), loc=1)
total_df.t <- as.matrix(t(total_df.s))
# plot!
Heatmap(total_df.t, name = "Relative cell abundance", row_order = c("CD8+ T cells", "CD4+ T cells", "Monocytes", "Macrophages", "B cells", "DC cells","NK cells", "Other"), column_order = get_gbm_order(unique(total_df$Samples)), width = unit(28, "cm"), height = unit(4, "cm"))
# column_order = gsub('-', '\\.', gbm_order)


