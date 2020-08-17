# 8.09.19
## TUMOR PURITY + DENSITY/COVERAGE PLOTS ##

# import packages
library(readxl)
library(tidyverse)
library(ggsci)
library(plyr); library(dplyr)
library(GGally)
library(UpSetR)
library(gridExtra)
library(gtable) 
library(scales)
library(openxlsx)

# set wd
setwd("~/Desktop/GBM_plots/CURRENT_PLOTS/tumor_purity/")

## TUMOR PURITY DENSITY PLOTS ##
# strategy for estimating tumor purity
# use only the following variants:
#1. 50x or greater coverage
#2. called by atleast 2 callers
#3. get rid of everything but chr1-22
#then take the median from all remaining variants

#gbm_batches <- c("Batch1", "Batch2", "Batch3")
gbm_batches <- c("Batch1", "Batch2", "Batch3", "TR")

tumor_purity_df <- data.frame()
variant_info_df <- data.frame()
for ( batch in gbm_batches ) {
  variant_path <- paste("~/Desktop/GBM/gbm_google_drive/", batch, "_annotated_variants.xlsx", sep="")
  samples <- excel_sheets(path = variant_path)
  tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[1])))
  for ( name in samples ) {
    # read in correct batch variant excel file
    df <- as.data.frame(read_excel(variant_path, sheet = name))
    # create standardized shorter sample names
     if ( grepl("Re", name) == FALSE && batch %in% c("Batch1", "Batch2") ) {
       new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
       tumor <- unlist(str_split(name, '-'))[2]
     } else if (batch == "Batch3") {
       new_name <- name
       tumor <- unlist(strsplit(name, '-'))[1]
     } else if (batch == "TR") {
       new_name <- paste(name, "TR", sep="_")
       tumor <- paste(unlist(str_split(name, '-'))[1], "TR", sep="_")
     } else {
       new_name <- paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
       tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")  
     } 
    df.subset <- df %>% 
      select(CHROM, POS, TUMOR.AD, TUMOR.DP, `Tumor VAF (%)`, set) 
    df.subset$`Tumor VAF (%)` <- as.numeric(df.subset$`Tumor VAF (%)`)
    df.subset$TUMOR.DP <- as.numeric(df.subset$TUMOR.DP)
    df.separate <- df.subset %>%
      separate(TUMOR.AD, into = c("REF.READS", "VAR.READS"), sep = ',')
    # count how many callers called the variant (column = set)
    set_nums <- c()
    for ( item in df.separate$set ) {
    set_nums <- c(set_nums, length(unlist(strsplit(item, "-"))))
    }
    df.final <- cbind(Samples = new_name, Tumors = tumor, df.separate, `set count` = set_nums)
    filtered_samples_df <- filter(df.final, CHROM != "chrX" & CHROM != "chrY" & TUMOR.DP >= 50 & `set count` > 1)
    variant_info_df <- rbind(variant_info_df, filtered_samples_df)
    tumor_purity <- median(filtered_samples_df$`Tumor VAF (%)`)*2
    purity_line <- data.frame(Sample = new_name, Tumor = tumor, `Total Variants` = nrow(df.final), `Filtered Variants` = nrow(filtered_samples_df), `Tumor Purity` = tumor_purity)
    tumor_purity_df <- rbind(tumor_purity_df, purity_line)
  }
}
# write df to excel file
write.xlsx(tumor_purity_df, "tumor_purity_median_VAF_estimates_all.xlsx")

# creating manual color palette for plotting
my_pal <- c(pal_npg("nrc", alpha = 1)(10), pal_simpsons("springfield", alpha = 1)(10))

# PLOTTING bar plot with tumor purity values
tumor_purity_df <- filter(tumor_purity_df, Filtered.Variants >= 10)
tumor_purity_df$Tumor <- factor(tumor_purity_df$Tumor, levels = unique(sort(as.character(tumor_purity_df$Tumor))))
tumor_purity_df <- tumor_purity_df[order(tumor_purity_df$Tumor),]
#order(levels(tumor_purity_df$Sample)

#tumor_purity_df$Sample[order()]

#df_celltypes.tidy <- df_celltypes.tidy[order(df_celltypes.tidy$sort.Samples.), ]

ggplot(tumor_purity_df, aes(x = Tumor, y = Tumor.Purity)) + geom_boxplot() + geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1))



# PLOTTING combined density/coverage plots
for ( tumor in unique(variant_info_df$Tumors) ) {
  #select_tumors <- c("BrMET025", "BrMET025tr", "BrMET028", "BrMET028tr")
  #if (tumor %in% select_tumors) {
    # subset a single tumor
  plotting_df <- filter(variant_info_df, Tumors == tumor)
  
  # calculating median VAF (same as above)
  median_values <- ddply(plotting_df, "Samples", summarise, value=median(`Tumor VAF (%)`))
  
  # plotting jitter w/ x=VAF, y=TUMOR.DP
  coverage_plot <- ggplot(plotting_df, aes(`Tumor VAF (%)`, TUMOR.DP, color = Samples)) + geom_point() + xlim(0,100) + ylab("Coverage") + scale_color_manual(values = my_pal) + theme(legend.position = "bottom")
  
  # plotting density plots w/ median vline x=VAF, y=density
  density_plot <- ggplot(plotting_df, aes(`Tumor VAF (%)`, fill = Samples, colour = Samples)) + geom_density(alpha=0.4) + xlim(0,100) + scale_fill_manual(values = my_pal) + scale_color_manual(values = my_pal) + geom_vline(data=median_values, aes(xintercept = value, color = Samples), linetype="dashed", size=0.8) + ggtitle(paste(tumor, "VAF Density Curve", sep = " ")) + theme(text = element_text(size=11), plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), legend.position = "none") + ylab("Density")
  
  # save combined plot
  grid_plot <- grid.arrange(density_plot, coverage_plot, nrow = 2)
  ggsave(grid_plot, file=paste(tumor, ".pdf", sep=""), width = 11, height = 8)
  #}
}

