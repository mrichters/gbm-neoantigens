library(UpSetR)
library(readxl)
library(tidyverse)
library(GGally)
setwd("~/Desktop/GBM/")
# for variants - use 'gbm_google_drive/GBM_variants_per_sample_v3.xlsx'
# for neoantigens class I - use 'gbm_google_drive/GBM_pvactools_clinicalclassI_11.20.18.xlsx'
# for neoantigens class II - use 'gbm_google_drive/GBM_pvactools_clinicalclassII_11.28.18.xlsx'

### OLD for variants ###

----------------
# 6.30.19
# variants
HTC_samples <- excel_sheets(path='gbm_google_drive/HTC_annotated_variants.xlsx')
HTC_tumors <- unique(sapply(HTC_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
TWCK_samples <- excel_sheets(path='gbm_google_drive/TWCK_annotated_variants.xlsx')
TWCK_tumors <- unique(sapply(TWCK_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
# create combined df for variants - run 
samples_df <- c()
for ( i in TWCK_samples ) {
  df <- as.data.frame(read_excel('gbm_google_drive/TWCK_annotated_variants.xlsx', range = cell_cols("A:B"), sheet = i))
  if ( is.null(dim(df)) ) {
    invisible()
  }
  df[[i]] <- 1
  df.gather <- gather(df, Sample, Value, 3)
  df.final <- unite(df.gather, Variant, CHROM, POS, sep=':')
  samples_df <- rbind(samples_df, df.final)
} 

#NEOANTIGENS#
HTC_samples <- excel_sheets(path='../gbm_google_drive/HTC_annotated_variants.xlsx')
HTC_tumors <- unique(sapply(HTC_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
TWCK_samples <- excel_sheets(path='../clinical_no_rna_filter/TWCK_clinical_nornafilter_classII_pvacseq_results.xlsx')
TWCK_tumors <- unique(sapply(TWCK_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
#i <- 'TWCK-BrMET018-BrMET018_Tumor_2'
for ( i in TWCK_samples ) {
  df <- as.data.frame(read_excel('../clinical_no_rna_filter/TWCK_clinical_nornafilter_classII_pvacseq_results.xlsx', sheet = i))
  if ( nrow(df) != 0 ) {
  df <- df[, c(10,11,12,15,19)]
  df[[i]] <- 1
  df <- df %>% gather(i, key = Sample, value = Value)
  df.final <- unite(df, Neoantigen, Mutation, `Protein Position`, `Gene Name`, `HLA Allele`, `MT Epitope Seq`, sep='_')
  assign(i, df.final)
  } 
}
# couldn't get this to work so doing manually below
#samples <- c(HTC_sample_names[1:4])
#samples.new <- sapply(samples, function(x) return(gsub('"', '`', x)))

#one for each tumor
sample <- distinct(rbind(`TWCK-GBM055-GBM055_Tumor_1`, `TWCK-GBM055-GBM055_Tumor_2`))#, `TWCK-GBM031-GBM031.Re_Tumor_4`))
# create separate columns for each sample (and removes redundancy of variant column)
sample <-spread(sample, Sample, Value)
# change NA values to 0
sample[is.na(sample)] <- 0
# shortening column names
sam_1 <- "GBM055-Tumor-1"
sam_2 <- "GBM055-Tumor-2"
#sam_3 <- "GBM031.Re-Tumor-4"
#sam_4 <- "BrMET018-Tumor-1"
colnames(sample) <- c("Variant", sam_1, sam_2)#, sam_3)
# write data frame to file
#write.csv(new_df, 'BrMET008_variant_intersection.csv')
# creating upset plot!
upset(sample, order.by = "degree", mainbar.y.label = "Neoantigen Intersections", sets.x.label = "Class II Neoantigens", sets = c(sam_2, sam_1), keep.order = TRUE, text.scale = 1.75, point.size = 4, line.size = 1.5)
------------------

## ? another way to handle above ? ##
# for loop to create combined data frame for all samples per tumor
#sample_names <- vector()

sample <- HTC_sample_names[1:4]
total_df <- data.frame()
for (i in sample) {
  sample_name <- paste(unlist(strsplit(i, '-'))[2], unlist(strsplit(i, '-'))[5], sep = "-")
  df <- as.data.frame(read_excel('../gbm_google_drive/GBM_HTC_annotated_variants.xlsx', sheet=sheet_name))
  df.new <- data.frame(df[, c(1,2,3,4,17)])
  df.new[[i]] <- 1
  colnames(df.new) <- c("CHROM", "POS", "REF", "ALT", "Gene", sample_name)
  df.gather <- gather(df.new, Sample, Value, 6)
  df.final <- unite(df.gather, Variant, CHROM, POS, sep=':')
  sample_names <- c(sample_names, sample_name)
  total_df <- rbind(total_df, df.final)
}
# create separate columns for each sample (and removes redundancy of variant column)
df.spread <- spread(total_df, Sample, Value)
# change NA values to 0
df.spread[is.na(df.spread)] <- 0
# write data frame to file
#write.csv(df.spread, '_variant_intersection.csv')

# create upset plots
upset(df.spread, order.by = "degree", mainbar.y.label = "Variant Intersections", sets.x.label = "Variants per Sample", sets = sample_names, keep.order = TRUE, text.scale = 1.75, point.size = 4, line.size = 1.5)

#create pairwise scatter plots
ggpairs(BrMET008.spr[,-c(1)], diag=list(continuous="bar"), upper=list(axisLabels='show'))




### NEW for variants ###
variant_df <- as.data.frame(read.csv('~/Desktop/GBM_data_plots_dfs/upset_plots_input_tables/Variant_Upset_Input/GBM052_variant_intersection.csv'))

variant_df$Sums <- rowSums(variant_df[, c(3:5)])
shared_var <- sum(variant_df$Sums == '3')
shared_partial_var <- sum(variant_df$Sums == '2') #+ sum(variant_df$Sums == '3')
unique_var <- sum(variant_df$Sums == '1')

table_df <- data.frame("H_TC-GBM052-052", shared_var, shared_partial_var, unique_var, nrow(variant_df))
colnames(table_df) <- c("Tumor", "Shared Variants", "Partially Shared Variants", "Unique Variants", "Total Variants")

final_table_df <- rbind(final_table_df, table_df)

# remove a row or column
final_table_df <- final_table_df[-c(3),]




### for neoantigens ###
colnames(df) <- c("Gene_Name", "Protein_Position", "Mutation", "HLA_Allele", "Neoantigen", sample_name)
df <- gather(df, Sample, Value, -Gene_Name, -Protein_Position, -Mutation, -HLA_Allele, -Neoantigen)

# defining variables
sample_nums <- c(4, 3, 2, 4, 3, 3, 3)
names(sample_nums) <- c("H_TC-BrMET008-008-Tumor-", "H_TC-BrMET009-009-Tumor-", "H_TC-GBM030-030-Tumor-", "H_TC-GBM032-032-Tumor-", "H_TC-GBM047-047-Tumor-", "H_TC-GBM051-051-Tumor-", "H_TC-GBM052-052-Tumor-")

TWCK_samples <- excel_sheets(path='../gbm_google_drive/GBM_TWCK_annotated_variants.xlsx')
TWCK_sample_nums <- c(3,3,3,3,4,3,3,2,3,3,3,3,3)
names(TWCK_sample_nums) <- unique(unlist(sapply(TWCK_samples, function(v) return(substr(v,1,nchar(v)-1)))))
#sheet_name <- "TWCK-GBM058-GBM058_Tumor_4"
# for loop to create combined data frame for all samples per tumor
sample_names <- vector()
total_df <- data.frame()
for ( tumor in TWCK_samples ) {
  #total_df <- data.frame()
  #if ( tumor != "TWCK-GBM058-GBM058_Tumor_" ) { 
  #for ( i in 1:TWCK_sample_nums[[tumor]] ) {
    sheet_name <- tumor#paste(tumor, i, sep = "")
    #sample_name <- gsub('-Tumor-', '', tumor)
    df <- as.data.frame(read_excel('../gbm_google_drive/GBM_TWCK_annotated_variants.xlsx', sheet=sheet_name))
    #df <- as.data.frame(read_excel('../gbm_google_drive/GBM_pvactools_clinicalclassII_11.28.18.xlsx', sheet=sheet_name))
    ### df <- rbind(df_I, df_II)
    df <- data.frame(df[, c(1,2,3,4,15)])
    #df$Sample <- sheet_name
    colnames(df) <- c("CHROM", "POS", "REF", "ALT", "VAF")#, "Sample")
    df$VAF <- df$VAF / 100
    #total_df <- rbind(total_df, df)
    #sample_names <- c(sample_names, sample_name)
    #assign(sample_name, total_df)
  #}
#new_df <- spread(total_df, Sample, VAF)
  new_df <- df
  new_df.sorted <- new_df %>% arrange(readr::parse_number(CHROM))
  new_df.sorted[is.na(new_df.sorted)] <- 0
  assign(sheet_name, new_df.sorted)
}
#}
scatterplot_names <- unique(sample_names)

# make ggpairs pairwise scatter plots from dfs

scatterplot_names <- unique(sample_names)
dataframes <- mget(scatterplot_names, envir=.GlobalEnv) 
result_list <- llply(dataframes, function(x) {
  pm <- ggpairs(new_df.sorted[,-c(1:4)], diag=list(continuous="bar"), upper=list(axisLabels='show'))
  pm2 <- pm
  for(i in 2:pm$nrow) {
    for(j in 1:(i-1)) {
      pm2[i,j] <- pm[i,j] +
        scale_x_continuous(limits = c(0, 1)) +
        scale_y_continuous(limits = c(0, 1))
    }
  }
  return(pm2)
})


# save data
write.csv(`H_TC-BrMET008-008`, 'BrMET008_pairscatter_variants.csv', row.names = FALSE)

#find variants with all these
brmet008.1 <- filter(`H_TC-BrMET008-008`, `H_TC-BrMET008-008-Tumor-1` > .8, `H_TC-BrMET008-008-Tumor-2` > .8)
brmet008.2 <- filter(`H_TC-BrMET008-008`, `H_TC-BrMET008-008-Tumor-1` > .8, `H_TC-BrMET008-008-Tumor-3` > .8)
brmet008.3 <- filter(`H_TC-BrMET008-008`, `H_TC-BrMET008-008-Tumor-1` > .8, `H_TC-BrMET008-008-Tumor-4` > .8)
brmet008.4 <- filter(`H_TC-BrMET008-008`, `H_TC-BrMET008-008-Tumor-2` > .8, `H_TC-BrMET008-008-Tumor-3` > .8)
brmet008.5 <- filter(`H_TC-BrMET008-008`, `H_TC-BrMET008-008-Tumor-2` > .8, `H_TC-BrMET008-008-Tumor-4` > .8)
brmet008.6 <- filter(`H_TC-BrMET008-008`, `H_TC-BrMET008-008-Tumor-3` > .8, `H_TC-BrMET008-008-Tumor-4` > .8)
brmet008 <- distinct(rbind(brmet008.1,brmet008.2,brmet008.3,brmet008.4,brmet008.5,brmet008.6))

write.csv(df, 'BrMET008_high_vafpairs.csv', row.names = FALSE)
# create separate columns for each sample (and removes redundancy of variant column)
#new_df <- spread(total_df, Sample, VAF)
# change NA values to 0
#new_df[is.na(new_df)] <- 0

dataframes <- mget(TWCK_samples, envir=.GlobalEnv)
histogram_list <- llply(dataframes, function(x) {
plot <- ggplot(x, aes(x=VAF)) + geom_histogram(binwidth = 0.02) + xlim(0,1) + ggtitle(x)
return(plot)
})


ggplot(`H_TC-GBM052-052`, aes(x=`H_TC-GBM052-052`$`H_TC-GBM052-052-Tumor-2`, y=`H_TC-GBM052-052`$`H_TC-GBM052-052-Tumor-3`)) + geom_point(size=2.5, color='navy') + xlab("GBM052-Tumor-2 VAF") + ylab("GBM052-Tumor-3 VAF") + theme(text = element_text(size=16))

# write dfs to excel file
write_xlsx(list(`H_TC-BrMET008-008`, `H_TC-BrMET009-009`, `H_TC-GBM030-030`, `H_TC-GBM032-032`, `H_TC-GBM047-047`, `H_TC-GBM051-051`, `H_TC-GBM052-052`), path = "~/Desktop/GBM_data_plots_dfs/pairwise_scatter_VAF_data.xlsx", col_names = TRUE, format_headers = TRUE)


new_df$Sums <- rowSums(new_df[, c(6:9)])

final_df <- data.frame()
for (row in 1:nrow(new_df)) {
  sum <- new_df[row, 'Sums']
  if ( sum == 4) {
   final_df <- rbind(final_df, new_df[row, ])
  }
}

####### extra code to make sumary table #######

#calculate the number of na that share all or 3,2,unique
new_df$Sums <- rowSums(new_df[, c(6:9)])
shared_na <- sum(new_df$Sums == '4')
shared_partial_na <- sum(new_df$Sums == '2') + sum(new_df$Sums == '3')
unique_na <- sum(new_df$Sums == '1')

#summary_df_I <- data.frame("H_TC-GBM052-052", shared_na, shared_partial_na, unique_na, nrow(new_df))
#colnames(summary_df_I) <- c("Tumor", "Shared Neoantigens Class I", "Partially Shared Neoantigens Class I", "Unique Neoantigens Class I", "Total Neoantigens Class I")

summary_df_II <- data.frame("H_TC-GBM052-052", shared_na, shared_partial_na, unique_na, nrow(new_df))
colnames(summary_df_II) <- c("Tumor", "Shared Neoantigens Class II", "Partially Shared Neoantigens Class II", "Unique Neoantigens Class II", "Total Neoantigens Class II")

table_df <- merge(summary_df_I, summary_df_II, by="Tumor")

final_table_df <- rbind(final_table_df, table_df)

final_table_df$'Number of Samples' <- c(4, 3, 2, 4, 3, 3, 3)
final_table_df <- final_table_df[,c(1,10,2,3,4,5,6,7,8,9)]
# also look at percentages

############ end of extra code #############

# creating upset plot!
upset(new_df, order.by = "degree", mainbar.y.label = "Neoantigen Intersections", sets.x.label = "Neoantigens per Sample", sets = c(sample_names[3], sample_names[2], sample_names[1]), keep.order = TRUE, text.scale = 1.75, point.size = 4, line.size = 1.5)
# save df

write.csv(final_table_df, 'variant_upset_plots_summary.csv')
