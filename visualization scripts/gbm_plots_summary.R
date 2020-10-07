# dimensions - ? currently vary depending on figure 

# theme - theme_bw() + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=13))

# color palettes - ggsci [ https://nanx.me/ggsci/articles/ggsci.html ] , R colors [ http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf ] , scico [ https://github.com/thomasp85/scico ]

# heatmap package - ComplexHeatmap, geom_tile

### Figure 1 ###
# variant counts #
# color palette - ggsci JAMA
# jitter and violin
ggplot(final_tumor_counts_df, aes(x=`Tumor Type`, y=Unique_Counts)) + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=13),axis.title.x = element_blank()) + geom_violin(fill = "grey90") + geom_jitter(position=position_jitter(0.1),aes(color=`Cancer Type`), size=3) + ylab("Variant Count") + scale_color_manual(values = pal_jama("default", alpha=0.9)(7)[1:7]) + coord_trans(y="log2") + scale_y_continuous(breaks=c(100,500,1000,2000,4000,6000))

# waterfall #
# color palette - picked manual colors for each variant type and clinical data info (met type, TERT promoter status)
GenVisR::waterfall(gbm_df, fileType = "Custom", variant_class_order = c("missense_variant", "frameshift_variant", "inframe_deletion", "start_lost", "stop_gained", "stop_lost", "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant"), mainPalette = c("stop_gained"='#4f00A8', "frameshift_variant"='#A80100', "inframe_deletion"='#ff9b34', "missense_variant"='#009933'), mainRecurCutoff = 0.05, maxGenes = 30, mainXlabel = TRUE, rmvSilent = TRUE, mutBurden = gbm_mut, mainDropMut = TRUE, sampOrder = sample_ordering_waterfall, clinLegCol=1, clinVarCol=c('Brain Metastasis'='#c2ed67', 'Primary GBM'='#E63A27', 'Recurrent GBM'='#e69127'), clinVarOrder=c('Brain Metastasis', 'Primary GBM', 'Recurrent GBM'))

# CNV by chromosome, tumor type #
# cp - R "red" and "blue"
# line and area
ggplot(data=gbm, aes(x=Start, y=G_score, color=Type)) + geom_line() + geom_area(aes(fill=Type)) + geom_hline(yintercept = c(0.24, -0.23), linetype="dashed") + facet_wrap( ~ Chromosome, scales="free_x", nrow = 1) + theme_bw() + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") + ylab("G-score") + scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), limits=c(-1.5,2.7)) + scale_color_manual(values = c("red", "blue")) + scale_fill_manual(values = c("red", "blue")) + ggtitle("GBM cohort")

## Figure 2,3 ##
# clonal proportion split boxplot
# color palette - R "firebrick4",  "royalblue3", "forestgreen" - how to match python colors?
# boxplot, dotplot
ggplot(clonal_df_spr, aes(x=`Tumor Type`, y=Proportion, fill=`Proportion Category`)) + geom_boxplot(position=position_dodge(0.75)) + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize=0.8) + theme_bw() + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=14), legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank()) + scale_y_continuous(expand=c(0,0)) + scale_fill_manual(values = c("firebrick4",  "royalblue3", "forestgreen"))

# upset plots #
# theme, colors NA

# CNV heatmaps by tumor type
# color palette - scico vikO
# ggtile
ggplot(df.g, aes(x=Samples, y=CNV, fill=`CNV Threshold`)) + geom_tile() + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.title.x.bottom=element_blank(), axis.ticks.y=element_blank(), legend.position = "right", text = element_text(size=13), plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio=0.75) + scale_fill_scico(palette = "vikO", limit = c(-2,2)) + ggtitle("BrMET NSCLC")

# CNV clonality
# will be done in python like variant/neoantigen clonality #

## Figure 4 ##
# Danaher immune heatmap #
# color palette - default of package
# ComplexHeatmap package
Heatmap(final_spread_matrix, name = "log2(TPM + 1)", row_order = row_order, width = unit(28, "cm"), height = unit(8, "cm"))

## Figure 6 ##
# GBM065 pair wise scatter plots 
# color palette - R "grey30" and ggsci Locus Zoom
# ggpairs package

# df.ggpairs$Gene <- factor(df.ggpairs$Gene, levels = c("Other", "TP53 R280G", "TP53 T125T", "IDH1 R132H", "MSH6 T1219I"))
# df.ggpairs <- df.ggpairs[order(df.ggpairs$Gene),]
# my_palette <- c("grey30", pal_locuszoom()(5)[c(1:4)])
pm <- ggpairs(df.ggpairs, diag=list(continuous="barDiag"), axisLabels='show', upper='blank')#, columns=1:5, mapping = ggplot2::aes(color=Gene), legend=c(2,1))
for(i in 2:pm$nrow) {
    for(j in 1:(i-1)) {
        pm[i,j] <- pm[i,j] +
            scale_x_continuous(limits = c(0, 1)) +
            scale_y_continuous(limits = c(0, 1)) #+
        #scale_color_manual(values = my_palette)
    }
}
# for(i in 1:5) {
#     for(j in 1:5) {
#         pm[i,j] <- pm[i,j] #+
#             #scale_fill_manual(values = my_palette)
#     }
# }
pm

# GBM065 sciclone / clonevol - in progress w/ ts data
# each package has default plots built in