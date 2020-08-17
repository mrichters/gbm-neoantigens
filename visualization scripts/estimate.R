# import packages
library(estimate)
library(gridExtra)
library(ggsci)
library(readxl)
library(tidyverse)
library(BBmisc)

# functions
get_gbm_order <- function(df) {
    recurrent_samples <- as.character(unique(df[grep("Re", df$Samples), ]$Samples))
    not_recurrent_samples <- as.character(unique(filter(df, !Samples %in% recurrent_samples)$Samples))
    sample_order <- c(not_recurrent_samples, recurrent_samples)
    tumor_order <- unique(sapply(sample_order, function(x) return(substr(x, 1, nchar(x)-2))))
    order_list <- list(sample_order, tumor_order)
    return(order_list)
}

# set wd
setwd("~/Desktop/GBM/immune_infiltration/estimate/")
# custom palette
my_pal <- c(pal_npg("nrc", alpha = 1)(10), pal_simpsons("springfield", alpha = 1)(10))
# looking at sample data
#help(package="estimate")
#OvarianCancerExpr <- system.file("extdata", "sample_input.txt", package="estimate")
#sample <- read.table(OvarianCancerExpr)

# read in input file
input_file <- "~/Desktop/GBM/immune_infiltration/matrix.abundance_gbm.txt"
input_df <- read.table(input_file, header=T)
# get sample names
sample_names <- colnames(input_df)
# step 1 - output_file: "estiamate.gct"
filterCommonGenes(input.f=input_file, output.f="estiamate.gct", id="GeneSymbol")
# step 2 - output_file: "estimate_score.gct"
estimateScore(input.ds="estiamate.gct", output.ds="estimate_score.gct", platform="illumina")
# this function doesn't work with illumina data
#plotPurity(scores = "estimate_score.gct", platform = "illumina", output.dir = "estimated_purity_plots/")

# create df from final output file
scores_df <- read.table(file = "estimate_score.gct", header=T, skip=2)
scores_df <- scores_df[, -1]
# change samples from "BrMET008.1" "BrMET008-1"
colnames(scores_df) <- sapply(colnames(scores_df), function(x) return(stri_replace_last_fixed(x, '.', '-')))
# gather samples
scores_df.tidy <- gather(scores_df, Samples, Score, -Description)
# spread 3 scores
scores_df.spread <- spread(scores_df.tidy, Description, Score)
scores_df.spread <- subset(scores_df.spread, select = -c(ESTIMATEScore))
#colnames(scores_df.spread) <- c("Samples", "Impurity Score", "Immunity Score", "Stromal Score")

df.normalize <- normalize(scores_df.spread, margin=2, range = c(0,1), method="range")
#df.normalize$Samples <- Samples
#df.normalize$`Purity Score` <- 1 - df.normalize$`Impurity Score` 

df.normalize.gather <- df.normalize %>%
                    gather(key = Score_Type, value = Score, -Samples)

gbm_order <- get_gbm_order(df.normalize)
df.normalize$Samples <- factor(df.normalize$Samples, levels = gbm_order[[1]])
df.normalize.order <- df.normalize[order(df.normalize$Samples), ]
df.normalize.order$Tumor_Type <- sapply(df.normalize.order$Samples, function(x) return(str_extract(x, "[a-zA-Z]+" )))
# now order the samples by value
df.normalize.scoreorder <- df.normalize[order(df.normalize$ImmuneScore), ]
df.normalize.scoreorder$Samples <- factor(df.noramlize.scoreorder$Samples, levels=df.noramlize.scoreorder$Samples)

# paired barplot
ggplot(df.normalize.scoreorder, aes(x = Samples, y = ImmuneScore, fill=Tumor_Type)) + geom_bar(stat = "identity", colour="black") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=12), axis.title.x=element_blank(), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), legend.position = "bottom") + scale_fill_manual(values=my_pal[3:4]) + ylab("Normalized Immune Score") 

#p2 <- ggplot(filter(scores_df.tidy2, Score != "Impurity Score" & Score != "Purity Score"), aes(x = Samples, y = Value, fill = Score)) + scale_fill_manual(values = my_pal[3:4]) + geom_bar(stat = "identity", colour="black") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Normalized Value")

#grid_plot <- grid.arrange(p1, p2, nrow = 2)
ggsave(grid_plot, filename = "estimate.pdf", width = 11, height = 8)
