library(tidyverse)
library(gridExtra)
library(readxl)
library(wesanderson)
library(viridis)
library(gtools)

get_gbm_order <- function(df) {
    recurrent_samples <- as.character(unique(df[grep("Re", df$Samples), ]$Samples))
    not_recurrent_samples <- as.character(unique(filter(df, !Samples %in% recurrent_samples)$Samples))
    sample_order <- c(not_recurrent_samples, recurrent_samples)
    tumor_order <- unique(sapply(sample_order, function(x) return(substr(x, 1, nchar(x)-4))))
    tumor_order_short <- sapply(tumor_order, function(x) return(unlist(strsplit(x, '-'))[2]))
    order_list <- tumor_order
    return(order_list)
}

# for input into this function, need to first fix file header
create_cnv_df <- function(file) {
    scores <- read.table(file, header=T)
    scores_tidy <- mutate(scores, G_score = ifelse(test = (Type == "Del"), yes = G_score * -1, no = G_score))
    for (row in 1:nrow(scores_tidy)) {
        if (scores_tidy[row, "q_value"] >= -log(0.10)) {
            scores_tidy[row, "Type"] == paste(scores_tidy[row, "Type"], "sig", sep = " ")
        }
    }
    return(scores_tidy)
}

setwd("~/Desktop/GBM/cnvkit/07.2020/")

files=c("gbm_0.1/scores.gistic", "brmet_lung/scores.gistic", "brmet_breast/scores.gistic")

sig_value <- 0.10
#G_value <- 0.24
gbm <- create_cnv_df(files[1])
gbm_amp <-filter(gbm, Type == "Amp" & q_value >= -log(sig_value))
gbm_amp_G <- mean(filter(gbm_amp, q_value == min(gbm_amp$q_value))$G_score)
gbm_del <- filter(gbm, Type == "Del" & q_value >= -log(sig_value))
gbm_del_G <- mean(filter(gbm_del, q_value == min(gbm_del$q_value))$G_score)

brmet <- create_cnv_df(files[2])
brmet_amp <-filter(brmet, Type == "Amp" & q_value >= -log(sig_value))
brmet_amp_G <- mean(filter(brmet_amp, q_value == min(brmet_amp$q_value))$G_score)
brmet_del <- filter(brmet, Type == "Del" & q_value >= -log(sig_value))
brmet_del_G <- mean(filter(brmet_del, q_value == min(brmet_del$q_value))$G_score)


p1 <- ggplot(data=gbm, aes(x=Start, y=G_score, color=Type)) + geom_line() + geom_area(aes(fill=Type)) + geom_hline(yintercept = c(0.24, -0.23), linetype="dashed") + facet_wrap( ~ Chromosome, scales="free_x", nrow = 1) + theme_bw() + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") + ylab("G-score") + scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), limits=c(-1.5,2.7)) + scale_color_manual(values = c("red", "blue")) + scale_fill_manual(values = c("red", "blue")) + ggtitle("GBM cohort")

p2 <- ggplot(data=brmet, aes(x=Start, y=G_score, color=Type)) + geom_line() + geom_area(aes(fill=Type)) + geom_hline(yintercept = c(0.78, -0.8), linetype="dashed") + facet_wrap( ~ Chromosome, scales="free_x", nrow = 1) + theme_bw() + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") + ylab("G-score") + scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3,4,5,6)) + scale_color_manual(values = c("red", "blue")) + scale_fill_manual(values = c("red", "blue")) + ggtitle("BrMET NSCLC cohort")

p_left <- grid.arrange(p1, p2, nrow = 2)

#### heatmap with threshold per sample (qvalue corresponds to peak itself) ####

cn_values <- as.data.frame(read_excel("brmet_breast/all_lesions.conf_75.xlsx"))
cn_values <- filter(cn_values, grepl("CN values",`Unique Name`))
cn_values <- cn_values[,c(2,10:ncol(cn_values))]
# gather values by peak descriptor
cn_values.gather <- gather(cn_values, key="Samples", value="CN Value", -Descriptor)
# convert columns to factors (Sample and Descriptor)
cn_values.gather$Samples <- factor(cn_values.gather$Samples, levels=unique(sort(cn_values.gather$Samples)))
#sorted_descriptor <- read.table("cn_values_brmet_sorted.tsv", col.names = "Descriptor")
cn_values.gather$Descriptor <- factor(cn_values.gather$Descriptor, levels=unique(mixedsort(cn_values.gather$Descriptor, decreasing = T)))
cn_values_brmet <- cn_values.gather 

pal = wes_palette("Zissou1", 100, type = "continuous")

p3 <- ggplot(cn_values_gbm, aes(x=Samples, y=Descriptor, fill=`CN Value`)) + geom_tile() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5), text = element_text(size=13)) + scale_fill_viridis(discrete=F) + ylab("Peaks") + ggtitle("GBM cohort") + labs(fill="CN Change")

p4 <- ggplot(cn_values_brmet, aes(x=Samples, y=Descriptor, fill=`CN Value`)) + geom_tile() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "right", plot.title = element_text(hjust = 0.5), text = element_text(size=13)) + scale_fill_viridis(discrete=F) + ggtitle("BrMET breast cancer cohort") + labs(fill="CN Change")

p_right <- grid.arrange(p3, p4, nrow = 1)
