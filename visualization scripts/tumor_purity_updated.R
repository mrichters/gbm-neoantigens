# 8.09.19

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
# set wd
setwd("~/Desktop/GBM_plots")

# creating manual color palette for plotting
my_pal <- c(pal_npg("nrc", alpha = 1)(10), pal_simpsons("springfield", alpha = 1)(10))
# for y axis range!!
#+ scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))

# list all sheets in excel file
twck_samples <- excel_sheets(path="~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx")
twck_tumors <- unique(sapply(twck_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
htc_samples <- excel_sheets(path="~/Desktop/GBM/gbm_google_drive/HTC_annotated_variants.xlsx")
htc_tumors <- unique(sapply(htc_samples, function(v) return(unlist(strsplit(v, '-'))[2])))

# all variants / neoantigen spreadsheets
# HTC no rna filter
variants <- "~/Desktop/GBM/gbm_google_drive/HTC_annotated_variants.xlsx"
classI <- "~/Desktop/GBM/gbm_google_drive/Batch1(HTC)_pvacseq/clinical_no_rna_filter/HTC_clinical_nornafilter_classI_pvacseq.xlsx"
classII <- "~/Desktop/GBM/gbm_google_drive/Batch1(HTC)_pvacseq/clinical_no_rna_filter/HTC_clinical_nornafilter_classII_pvacseq.xlsx"
# HTC
variants <- "~/Desktop/GBM/gbm_google_drive/HTC_annotated_variants.xlsx"
classI <- "~/Desktop/GBM/gbm_google_drive/Batch1(HTC)_pvacseq/clinical/HTC_clinical_classI_pvacseq.xlsx"
classII <- "~/Desktop/GBM/gbm_google_drive/Batch1(HTC)_pvacseq/clinical/HTC_clinical_classII_pvacseq.xlsx"
# TWCK no rna filter
variants <- "~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx"
classI <- "~/Desktop/GBM/gbm_google_drive/Batch2(TWCK)_pvacseq/clinical_no_rna_filter/TWCK_clinical_nornafilter_classI_pvacseq.xlsx"
classII <- "~/Desktop/GBM/gbm_google_drive/Batch2(TWCK)_pvacseq/clinical_no_rna_filter/TWCK_clinical_nornafilter_classII_pvacseq.xlsx"
# TWCK
variants <- "~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx"
classI <- "~/Desktop/GBM/gbm_google_drive/Batch2(TWCK)_pvacseq/clinical/TWCK_clinical_classI_pvacseq.xlsx"
classII <- "~/Desktop/GBM/gbm_google_drive/Batch2(TWCK)_pvacseq/clinical/TWCK_clinical_classII_pvacseq.xlsx"

## TUMOR PURITY DENSITY PLOTS ##
# run this once for htc and once for twck to create a combined df
samples_df <- c()
for ( name in twck_samples ) {
  df <- as.data.frame(read_excel("~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx", sheet = name))
  if ( grepl("Re", name) == FALSE ) {
    new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
    tumor <- unlist(str_split(name, '-'))[2]
  } else {
    new_name <- paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
    tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")  
  } 
  df.subset <- df %>% 
    select(CHROM, POS, TUMOR.AD, TUMOR.DP, `Tumor VAF (%)`) 
  df.separate <- df.subset %>%
                    separate(TUMOR.AD, into = c("REF.READS", "VAR.READS"), sep = ',')
  df.final <- cbind(Samples = new_name, Tumors = tumor, df.separate)
  #if ( length(row.names(df.separate)) >= 20 ) {
  samples_df <- rbind(samples_df, df.final)
  #}
}
for ( tumor in as.vector(unique(samples_df$Tumors)) ) {
  samples_df.subset <- samples_df %>%
                        filter(Tumors == tumor)
  # plotting jitter w/ x=VAF, y=TUMOR.DP
  point_name <- paste(tumor, ".coverage.pdf", sep="")
  point_plot <- ggplot(samples_df.subset, aes(`Tumor VAF (%)`, TUMOR.DP, colour = Samples)) + geom_point() + xlim(0,100) + ylab("Coverage") + scale_color_manual(values = my_pal) + theme(legend.position = "bottom")
  #ggsave(point_plot, file=point_name, width = 8, height = 4)
  # plotting density plots w/ median vline x=VAF, y=density
  median_values <- ddply(samples_df.subset, "Samples", summarise, value=median(`Tumor VAF (%)`))
  density_name <- paste(tumor, ".pdf", sep="")
  density_plot <- ggplot(samples_df.subset, aes(`Tumor VAF (%)`, fill = Samples, colour = Samples)) + geom_density(adjust=1/3, alpha=0.1) + xlim(0,100) + scale_fill_manual(values = my_pal) + scale_color_manual(values = my_pal) + geom_vline(data=median_values, aes(xintercept = value, color = Samples), linetype="dashed") + ggtitle(paste(tumor, "VAF Density Curve", sep = " ")) + theme(text = element_text(size=11), plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), legend.position = "none") + ylab("Density")
  #ggsave(density_plot, file=density_name, width = 8, height = 5)
  grid_plot <- grid.arrange(density_plot, point_plot, nrow = 2)
  ggsave(grid_plot, file=density_name, width = 11, height = 8)
}


## VARIANT UPSET PLOTS ##
total_shared_value_df <- data.frame()
for ( tumor in htc_tumors ) {
  tumor_df <- c()
  sample_names <- c()
  for ( name in htc_samples[grep(tumor, htc_samples)] ) {
    df <- as.data.frame(read_excel("~/Desktop/GBM/gbm_google_drive/HTC_annotated_variants.xlsx", sheet = name))
    if ( grepl("Re", name) == FALSE ) {
      new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
      tumor <- unlist(str_split(name, '-'))[2]
    } else {
      new_name <-  paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
      tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
    } 
    sample_names <- c(sample_names, new_name)
    df.subset <- df %>% 
      select(CHROM, POS, REF, ALT) 
    df.subset.cbind <- cbind(df.subset, Sample = new_name, Value = 1)
    df.final <- df.subset.cbind %>%
                  unite(Variant, CHROM, POS, REF, ALT, sep=':')
    tumor_df <- rbind(tumor_df, df.final)
  }
  df.spread <- tumor_df %>% 
                spread(Sample, Value)
  df.spread[is.na(df.spread)] <- 0
  counter = 0
  shared_value_vector <- c()
  for ( line in 1:nrow(df.spread) ) {
    nums_only_df <- df.spread[,-1]
    shared_value <- sum(nums_only_df[line, ])
    print(shared_value)
    shared_value_vector <- c(shared_value_vector, shared_value)
    if ( shared_value >= 2 ) {
      counter = counter + 1
    }
  }
  df.spread$Shared <- shared_value_vector
  clonal <- max(df.spread$Shared)
  subclonal_private <- 1
  subclonal_shared <- (subclonal_private+1):(clonal-1)
  if ( length(subclonal_shared) == 1 ) {
    subclonal_shared <- c(subclonal_shared, 0)
  }
  shared_categories <- c()
  for ( line in 1:nrow(df.spread) ) {
    if ( df.spread$Shared[line] == clonal ) {
      category <- "Clonal"
    } else if ( df.spread$Shared[line] == subclonal_shared[1] ) {
      category <- "Subclonal Shared"
    } else if ( df.spread$Shared[line] == subclonal_shared[2] ) {
      category <- "Subclonal Shared"
    } else if ( df.spread$Shared[line] == subclonal_private ) {
      category <- "Subclonal Private"
    }
    shared_categories <- c(shared_categories, category)
  }
  df.spread$Shared_category <- shared_categories
  shared_df <- cbind(data.frame(tumor), df.spread$Shared, df.spread$Shared_category)
  colnames(shared_df) <- c("Tumor", "Shared Value", "Variant Type")
  add_to_finaldf <- as.data.frame(table(shared_df$`Variant Type`))
  add_to_finaldf$Tumor <- tumor
  total_shared_value_df <- rbind(total_shared_value_df, add_to_finaldf)
  # create upset plots
  #variant_upset <- upset(df.spread, order.by = "degree", mainbar.y.label = "Variant Intersections", sets.x.label = "Variants per Sample", sets = sort(sample_names, decreasing = TRUE), keep.order = TRUE, text.scale = 1.75, point.size = 4, line.size = 1.5)
  #pdf(paste(tumor, ".variant.upset.pdf", sep = ""), width = 10, height = 8)
  #print(variant_upset)
  #dev.off()
}
# change order of the Samples
total_shared_value_df$Tumor <- as.factor(total_shared_value_df$Tumor)
colnames(total_shared_value_df) <- c("Variant Type", "Variant Count", "Tumor")
# for proportions
#function
#df[, -1] <- lapply( df[ , -1], function(x) x/sum(x, na.rm=TRUE) )
total_shared_value_df.spread <- total_shared_value_df %>%
                                  spread(key = Tumor, value = `Variant Count`)
total_shared_value_df.spread[is.na(total_shared_value_df.spread)] <- 0
total_shared_value_df.spread[, -1] <- lapply( total_shared_value_df.spread[ , -1], function(x) x/sum(x, na.rm=TRUE) )
total_shared_value_df.tidy <- total_shared_value_df.spread %>%
                                gather(key = "Tumor", value = "Variant Proportion", -`Variant Type`)
# ordering by proportion
#sample_levels <- names(sort(apply(total_shared_value_df.spread[,2:ncol(total_shared_value_df.spread)], 2, mean),decreasing = TRUE))
# reorder df
#total_shared_value_df.tidy$Samples <- as.factor(total_shared_value_df.tidy$Samples)
#levels(total_shared_value_df.tidy$Samples) <- sample_levels
# raw counts
ggplot(total_shared_value_df, aes(x=Tumor, y=`Variant Count`, fill=`Variant Type`)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=14)) + scale_fill_manual(values=my_pal[c(1,3,4)]) + ylab("Variant Count")
# proportions
ggplot(total_shared_value_df.tidy, aes(x=Tumor, y=`Variant Proportion`, fill=`Variant Type`)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=14)) + scale_fill_manual(values=my_pal[c(1,3,4)]) + ylab("Variant Type Proportion")

# calculate shared variants and create final df for shared variants proportion
counter = 0
shared_value_vector <- c()
for ( line in 1:nrow(df.spread) ) {
  nums_only_df <- df.spread[,-1]
  shared_value <- sum(nums_only_df[line, ])
  print(shared_value)
  shared_value_vector <- c(shared_value_vector, shared_value)
  if ( shared_value >= 2 ) {
    counter = counter + 1
  }
}
if ( counter == 0 ) {
  final_shared_value <- 0
} else {
  final_shared_value <- as.numeric(counter / nrow(df.spread))
}
sample_line <- cbind(data.frame(tumor), data.frame(final_shared_value))
total_shared_value_df <- rbind(total_shared_value_df, sample_line)
colnames(total_shared_value_df) <- c("Tumor", "Proportion of Shared Variants")
total_lines <- c()
for ( line in total_shared_value_df[,1] ) {
  if ( grepl("Br", line) == TRUE ) {
    type <- "Brain Met"
  } else if ( grepl("Re", line) == TRUE ) {
    type <- "Recurrent GBM"
  } else {
    type <- "Primary GBM"
  }
  total_lines <- c(total_lines, type)
}
total_shared_value_df$`Tumor Type` <- total_lines
# order the x axis by proportion
total_df.order <- as.vector(total_shared_value_df[order(total_shared_value_df$`Proportion of Shared Variants`), ][, 1])
# create lollipop plot
ggplot(total_shared_value_df, aes(Tumor, `Proportion of Shared Variants`)) + theme_bw() + geom_segment(aes(x=Tumor, xend=Tumor, y=0, yend=`Proportion of Shared Neoantigens`)) + scale_x_discrete(limits=total_df.order) + geom_point(size=5, aes(color=`Tumor Type`)) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + scale_color_manual(values = my_pal)

## NEOANTIGEN UPSET PLOTS ##
total_shared_value_df <- c()
for ( tumor in twck_tumors ) {
  tumor_df <- c()
  sample_names <- c()
  for ( name in twck_samples[grep(tumor, twck_samples)] ) {
    # combined class I and II neoantigens
    df <- as.data.frame(rbind(read_excel(classI, sheet = name, range = cell_cols("A:T")), read_excel(classII, sheet = name, range = cell_cols("A:T"))))
    # class I only
    #df <- as.data.frame(rbind(read_excel(classI, sheet = name, range = cell_cols("A:T"))))
    # class II only
    #df <- as.data.frame(rbind(read_excel(classII, sheet = name, range = cell_cols("A:T"))))
    if ( grepl("Re", name) == FALSE ) {
      new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
      tumor <- unlist(str_split(name, '-'))[2]
    } else {
      new_name <-  paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
      tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
    } 
    if ( nrow(df) != 0 ) {
      sample_names <- c(sample_names, new_name)
      df.subset <- df %>% 
        # for total neoantigens
        select(Mutation, `Protein Position`, `Gene Name`, `MT Epitope Seq`)
        # for class I and II neoantigens
        #select(Mutation, `Protein Position`, `Gene Name`, `HLA Allele`, `MT Epitope Seq`) 
      df.subset.cbind <- cbind(df.subset, Sample = new_name, Value = 1)
      df.final <- df.subset.cbind %>%
        # for total neoantigens
        unite(Neoantigen, Mutation, `Protein Position`, `Gene Name`, `MT Epitope Seq`, sep='_')
        # for class I and II neoantigens
        #unite(Neoantigen, Mutation, `Protein Position`, `Gene Name`, `HLA Allele`, `MT Epitope Seq`, sep='_')
      tumor_df <- rbind(tumor_df, unique(df.final))
    }
  }
  if ( nrow(count(unique(as.vector(tumor_df$Sample)))) > 0 ) {
  df.spread <- tumor_df %>% 
    spread(Sample, Value)
  df.spread[is.na(df.spread)] <- 0
  # calculate shared neoantigens and create final df
  counter = 0
  shared_value_vector <- c()
  for ( line in 1:nrow(df.spread) ) {
    nums_only_df <- as.data.frame(df.spread[,-1])
    shared_value <- sum(nums_only_df[line, ])
    print(shared_value)
    shared_value_vector <- c(shared_value_vector, shared_value)
    if ( shared_value >= 2 ) {
      counter = counter + 1
    }
  }
  df.spread$Shared <- shared_value_vector
  clonal <- max(df.spread$Shared)
  subclonal_private <- 1
  subclonal_shared <- (subclonal_private+1):(clonal-1)
  if ( length(subclonal_shared) == 1 ) {
    subclonal_shared <- c(subclonal_shared, 0)
  }
  shared_categories <- c()
  for ( line in 1:nrow(df.spread) ) {
    if ( df.spread$Shared[line] == clonal ) {
      category <- "Clonal"
    } else if ( df.spread$Shared[line] == subclonal_shared[1] ) {
      category <- "Subclonal Shared"
    } else if ( df.spread$Shared[line] == subclonal_shared[2] ) {
      category <- "Subclonal Shared"
    } else if ( df.spread$Shared[line] == subclonal_private ) {
      category <- "Subclonal Private"
    }
    shared_categories <- c(shared_categories, category)
  }
  df.spread$Shared_category <- shared_categories
  shared_df <- cbind(data.frame(tumor), df.spread$Shared, df.spread$Shared_category)
  colnames(shared_df) <- c("Tumor", "Shared Value", "Variant Type")
  add_to_finaldf <- as.data.frame(table(shared_df$`Variant Type`))
  add_to_finaldf$Tumor <- tumor
  total_shared_value_df <- rbind(total_shared_value_df, add_to_finaldf)
  }
}
#colnames(total_shared_value_df) <- c("Variant Type", "Variant Count", "Tumor")
total_shared_value_df <- total_shared_value_df[order(total_shared_value_df$Tumor), ]
# create shared stacked barplot
# change order of the Samples
total_shared_value_df$Tumor <- as.factor(total_shared_value_df$Tumor)
colnames(total_shared_value_df) <- c("Neoantigen Type", "Neoantigen Count", "Tumor")
# for proportions
#function
#df[, -1] <- lapply( df[ , -1], function(x) x/sum(x, na.rm=TRUE) )
total_shared_value_df.spread <- total_shared_value_df %>%
  spread(key = Tumor, value = `Neoantigen Count`)
total_shared_value_df.spread[is.na(total_shared_value_df.spread)] <- 0
total_shared_value_df.spread[, -1] <- lapply( total_shared_value_df.spread[ , -1], function(x) x/sum(x, na.rm=TRUE) )
total_shared_value_df.tidy <- total_shared_value_df.spread %>%
  gather(key = "Tumor", value = "Neoantigen Proportion", -`Neoantigen Type`)
# proportions
p1 <- ggplot(total_shared_value_df.tidy, aes(x=Tumor, y=`Neoantigen Proportion`, fill=`Neoantigen Type`)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=14)) + scale_fill_manual(values=my_pal[c(1,3,4)]) + ylab("Neoantigen Type Proportion")
# raw counts
ggplot(total_shared_value_df, aes(x=Tumor, y=`Neoantigen Count`, fill=`Neoantigen Type`)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=14)) + scale_fill_manual(values=my_pal[c(1,3,4)]) + ylab("Neoantigen Count (with RNA filter)")
### COMBINED PLOTS FOR RAW COUNTS / PROPORTIONS
grid.arrange(p2, p1, nrow=2)
########### FACET WRAP TRIED AND FAILED ################
# break into facetwrap
facet_vector <- c()
for ( line in 1:nrow(total_shared_value_df) ) {
  if ( grepl("BrMET", total_shared_value_df$Tumor[line]) == TRUE ) {
    facet <- "Facet1"
  } else if ( grepl("GBM065", total_shared_value_df$Tumor[line]) == TRUE ) {
    facet <- "Facet3"
  } else {
    facet <- "Facet2"
  }
  facet_vector <- c(facet_vector, facet)
}
total_shared_value_df$Facet <- facet_vector
# plot with facet wrap
breaks = 10**(1:10 * 0.5)
ggplot(total_shared_value_df, aes(x=Tumor, y=`Neoantigen Count`, fill=`Neoantigen Type`)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=14), legend.position = "none") + scale_fill_manual(values=my_pal[c(1,3,4)]) + ylab("Neoantigen Count") + coord_flip()

p1 <- ggplot(filter(total_shared_value_df, Facet == "Facet1"), aes(x=Tumor, y=`Neoantigen Count`, fill=`Neoantigen Type`)) + geom_bar(stat = "identity", width = 0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=14), legend.position = "none") + scale_fill_manual(values=my_pal[c(1,3,4)]) + ylab("Neoantigen Count")

p2 <- ggplot(filter(total_shared_value_df, Facet == "Facet2"), aes(x=Tumor, y=`Neoantigen Count`, fill=`Neoantigen Type`)) + geom_bar(stat = "identity", width = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=14), axis.title.y = element_blank(), legend.position = "none") + scale_fill_manual(values=my_pal[c(1,3,4)]) + ylab("Neoantigen Count")

p3 <- ggplot(filter(total_shared_value_df, Facet == "Facet3"), aes(x=Tumor, y=`Neoantigen Count`, fill=`Neoantigen Type`)) + geom_bar(stat = "identity", width = 0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=14), axis.title.y = element_blank(), axis.title.x = element_blank()) + scale_fill_manual(values=my_pal[c(1,3,4)])
  
grid.arrange(p1, p2, p3, nrow=1)
ggarrange(p1, p2, p3, widths = c(1.5,2))

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g <- rbind(g2, g3, size = "first")
g$widths <- unit.pmax(g2$widths, g3$widths)
grid.newpage()
grid.draw(g)
###########################################
# for upset plots/lollipop proportion
  if ( counter == 0 ) {
    final_shared_value <- 0
  } else {
    final_shared_value <- as.numeric(counter / nrow(df.spread))
  }
  sample_line <- cbind(data.frame(tumor), data.frame(final_shared_value))
  total_shared_value_df <- rbind(total_shared_value_df, sample_line)
  # create upset plots
  variant_upset <- upset(df.spread, order.by = "degree", mainbar.y.label = "Neoantigen Intersections", sets.x.label = "Class I Neoantigens", sets = sort(sample_names, decreasing = TRUE), keep.order = TRUE, text.scale = 1.75, point.size = 4, line.size = 1.5)
  pdf(paste(tumor, ".neoantigenI.upset.pdf", sep = ""), width = 10, height = 8)
  print(variant_upset)
  dev.off()
  }
}
colnames(total_shared_value_df) <- c("Tumor", "Proportion of Shared Neoantigens")
total_lines <- c()
for ( line in total_shared_value_df[,1] ) {
  if ( grepl("Br", line) == TRUE ) {
    type <- "Brain Met"
  } else if ( grepl("Re", line) == TRUE ) {
    type <- "Recurrent GBM"
  } else {
    type <- "Primary GBM"
  }
  total_lines <- c(total_lines, type)
}
total_shared_value_df$`Tumor Type` <- total_lines
# order the x axis by proportion
total_df.order <- as.vector(total_shared_value_df[order(total_shared_value_df$`Proportion of Shared Neoantigens`), ][, 1])
# create lollipop plot
ggplot(total_shared_value_df, aes(Tumor, `Proportion of Shared Neoantigens`)) + theme_bw() + geom_segment(aes(x=Tumor, xend=Tumor, y=0, yend=`Proportion of Shared Neoantigens`)) + scale_x_discrete(limits=total_df.order) + geom_point(size=5, aes(color=`Tumor Type`)) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + scale_color_manual(values = my_pal)



## CREATE X=TOTAL NEOANTIGENS, Y=VARIANTS ##
# initialize df: sample, variants, total neoantigens
scatter_df <- c()
pvacseq_scatter_df <- c()
all_scatter_df <- c()
for ( tumor in twck_tumors ) {
  tumor_df <- c()
  sample_names <- c()
  for ( name in twck_samples[grep(tumor, twck_samples)] ) {
    variants_df <- as.data.frame(read_excel(variants, sheet = name))
    neoantigens_df <- as.data.frame(rbind(read_excel(classI, sheet = name, range = cell_cols("A:T")), read_excel(classII, sheet = name, range = cell_cols("A:T"))))
    neoantigensI_df <- as.data.frame(rbind(read_excel(classI, sheet = name, range = cell_cols("A:T"))))
    neoantigensII_df <- as.data.frame(rbind(read_excel(classII, sheet = name, range = cell_cols("A:T"))))
    if ( grepl("Re", name) == FALSE ) {
      new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
      tumor <- unlist(str_split(name, '-'))[2]
    } else { 
      new_name <- paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
      tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
    } 
    sample_names <- c(sample_names, new_name)
    # for creating master figure w/ mutational signatures: need PROPORTION of each variant type that i'm getting from the variants_df
    # for df, every consequence is a column and value is relative proportion (count / total)
    variants_by_conseq = 0
    vbc_vector <- c()
    for ( consequence in variants_df$Consequence ) {
      if ( grepl("missense", consequence) || grepl("inframe", consequence) || grepl("frameshift", consequence) || grepl("protein_altering", consequence) ) {
        variants_by_conseq = variants_by_conseq + 1
        vbc_vector <- c(vbc_vector, consequence)
      } 
    }
    # for scatter_df
    final_df_line <- cbind(data.frame(new_name), data.frame(tumor), data.frame(nrow(variants_df)), data.frame(nrow(neoantigens_df)))
    colnames(final_df_line) <- c("Samples", "Tumors", "Variants", "Total Neoantigens")
    scatter_df <- rbind(scatter_df, final_df_line)
    # for pvacseq_scatter_df
    pvacseq_final_df_line <- cbind(data.frame(new_name), data.frame(tumor), data.frame(variants_by_conseq), data.frame(nrow(neoantigens_df)))
    colnames(pvacseq_final_df_line) <- c("Samples", "Tumors", "pVACseq-compatible Variants", "Total Neoantigens")
    pvacseq_scatter_df <- rbind(pvacseq_scatter_df, pvacseq_final_df_line)
    # all_scatter_df
    all_scatter_df_line <- cbind(data.frame(new_name), data.frame(tumor), data.frame(variants_by_conseq), data.frame(nrow(neoantigensI_df)), data.frame(nrow(neoantigensII_df)))
    colnames(all_scatter_df_line) <- c("Samples", "Tumors", "pVACseq-compatible Variants", "Class I Neoantigens", "Class II Neoantigens")
    all_scatter_df <- rbind(all_scatter_df, all_scatter_df_line)
  }
}
colnames(scatter_df) <- c("Samples", "Tumors", "Variants", "Total Neoantigens")
colnames(pvacseq_scatter_df) <- c("Samples", "Tumors", "pVACseq-compatible Variants", "Total Neoantigens")
colnames(all_scatter_df) <- c("Samples", "Tumors", "pVACseq-compatible Variants", "Class I Neoantigens", "Class II Neoantigens")
# order df by Samples
scatter_order <- scatter_df[order(as.character(scatter_df$Tumors)),][,2]
levels(scatter_df$Tumors) <- as.vector(unique(scatter_order))
# create total variants scatter plot
ggplot(scatter_df, aes(x=`Total Neoantigens`, y=Variants, color = Tumors)) + geom_jitter(size=5) + theme(text = element_text(size=16)) + scale_color_manual(values = my_pal) + ylab("Total Variants")

# order df by Samples
scatter_order <- pvacseq_scatter_df[order(as.character(pvacseq_scatter_df$Tumors)),][,2]
levels(pvacseq_scatter_df$Tumors) <- as.vector(unique(scatter_order))
# create pvacseq compatible variants scatter plot
ggplot(pvacseq_scatter_df, aes(x=`Total Neoantigens`, y=`pVACseq-compatible Variants`)) + geom_jitter(size=5, aes(color=Tumors)) + theme(text = element_text(size=16)) + scale_color_manual(values = my_pal) + xlab("Neoantigens") #+ geom_smooth(method=lm)
# create combined total/pvacseq scatter df
combined_scatter_df <- cbind(all_scatter_df[,c(1:3)], `Total Variants` = scatter_df[,3], all_scatter_df[,c(4:5)], `Total Neoantigens` = scatter_df[,4])
write_excel_csv(combined_scatter_df, path = "~/Desktop/GBM_plots/master_df.csv")
# same plot but no colors and best fit line
ggplot(pvacseq_scatter_df, aes(x=`Total Neoantigens`, y=`pVACseq-compatible Variants`)) + geom_jitter(size=5) + theme(text = element_text(size=16)) + xlab("Neoantigens") + geom_smooth(method=lm)
ggplot(scatter_df, aes(x=`Total Neoantigens`, y=Variants)) + geom_jitter(size=3) + theme(text = element_text(size=16)) + ylab("Total Variants") + geom_smooth(method=lm)
## FIND R SQUARED VALUE OF VARIANTS / NEOANTIGENS CORRELATION
scatter_df.lm = lm(Variants ~ `Total Neoantigens`, data=scatter_df)
summary(scatter_df.lm)$adj.r.squared
summary(scatter_df.lm) #<-- tells us p-value and whether correlation is significant

## GGPAIRS VAF SCATTER PLOTS ##
for ( tumor in twck_tumors ) {
  tumor_df <- c()
  for ( name in twck_samples[grep(tumor, twck_samples)] ) {
    df <- as.data.frame(read_excel("~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx", sheet = name))
    if ( grepl("Re", name) == FALSE ) {
      new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
      tumor <- unlist(str_split(name, '-'))[2]
    } else {
      new_name <-  paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
      tumor <- str_sub(unlist(str_split(name, '-'))[2], 1, 9)  
    } 
    df.subset <- df %>% 
      select(CHROM, POS, REF, ALT, `Tumor VAF (%)`)
    df.subset$`Tumor VAF (%)` <- df.subset$`Tumor VAF (%)` / 100
    names(df.subset)[5] <- "VAF"
    df.subset.cbind <- cbind(df.subset, Sample = new_name)
    df.final <- df.subset.cbind %>%
      unite(Variant, CHROM, POS, REF, ALT, sep=':')
    tumor_df <- rbind(tumor_df, df.final)
  }
  df.spread <- tumor_df %>% 
    spread(Sample, VAF)
  df.spread[is.na(df.spread)] <- 0
  df.ggpairs <- df.spread[,-c(1)]
  #pdf(paste(tumor, ".vaf.scatter.pdf", sep = ""), width = 10, height = 8)
  pm <- ggpairs(df.ggpairs, diag=list(continuous="barDiag"), axisLabels='show')
  pm2 <- pm
  for(i in 2:pm$nrow) {
    for(j in 1:(i-1)) {
      pm2[i,j] <- pm[i,j] +
        scale_x_continuous(limits = c(0, 1)) +
        scale_y_continuous(limits = c(0, 1))
    }
  }
  pdf(paste(tumor, "_vaf_scatter.pdf", sep = ""), width = 10, height = 8)
  print(pm2)
  dev.off()
}


## load combined_scatter_df back as an excel file! ~/Desktop/GBM_plots/gbm_variant_neoantigen_counts.xlsx
combined_df <- read_excel("~/Desktop/GBM_plots/most_recent_plots/gbm_variant_neoantigen_counts.xlsx")
combined_df.gather <- combined_df %>%
                        gather(`pVACseq-compatible Variants`, `Class I Neoantigens`, `Class II Neoantigens`, key = Category, value = Count)
ggplot(combined_df.gather, aes(x = Samples, y = Count)) + scale_color_manual(values = my_pal) + geom_bar(stat = "identity", aes(fill = Category))  + theme(legend.position="top", axis.text.x = element_text(angle = 45, hjust = 1))

