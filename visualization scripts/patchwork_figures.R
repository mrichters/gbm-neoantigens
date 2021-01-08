library(tidyverse)
library(patchwork)

setwd("~/Desktop/GBM Figures/Rdata/Figure3")

for (f in list.files()) {
    load(f)
}

# Figure 1 done in Adobe bc wf formats didn't work w/ layout

# Figure 2
load("3D.Rdata")
# cnv_clonality
# cnv_heatmap
# BrMET019_upset
# GBM047.Re_upset
# clonality
# clonality_boxplot

(clonality + clonality_boxplot + plot_layout(guides = "collect")) / (Fig2c + Fig2d + plot_spacer() + plot_layout(widths = c(1,2))) / cnv_heatmap / (cnv_clonality + plot_spacer() + plot_layout(ncol = 4)) + plot_layout(ncol = 1, heights = c(1.5,1,1.8,1)) + plot_annotation(tag_levels = 'A')

fig2_ab <- clonality + clonality_boxplot + plot_layout(guides = "collect")
fig2_cd <- Fig2c + Fig2d + plot_layout(widths = c(1,2))
cnv_heatmap
cnv_clonality 

# Figure 3
na_clonality + na_clonal_prop + plot_layout(widths = c(2,1.5))

Fig3c + Fig3d + plot_layout(widths = c(1,2))


