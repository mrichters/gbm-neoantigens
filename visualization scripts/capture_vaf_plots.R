setwd("~/Desktop/GBM/")

new_vafs <- as.data.frame(read_excel("variants.annotated.tsv.variants_annotated_vaf.xlsx"))
old_vafs <- as.data.frame(read_excel("gbm_google_drive/Batch2_annotated_variants.xlsx"), sheet = "TWCK-BrMET010-BrMET010_Tumor_3")

vaf_df <- data.frame("Capture" = new_vafs$`Tumor VAF (%)`, "Original" = old_vafs$`Tumor VAF (%)`)

vafs_df <- as.data.frame(read_excel("vafs.xlsx"))

vafs_df.gather <- gather(vafs_df, VAF, Value)

ggplot(vafs_df.gather, aes(x=Value)) + geom_histogram(data=subset(vafs_df.gather, VAF== 'Original'), fill="red", alpha = 0.4) + geom_histogram(data=subset(vafs_df.gather, VAF== 'Capture'), fill="blue", alpha = 0.4) + xlab("VAF") + geom_vline(xintercept = c(22.19982,22.74742), linetype = "dashed")
