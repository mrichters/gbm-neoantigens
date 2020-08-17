library(tidyverse)
library(readxl)

data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
          sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
}

setwd("~/Desktop/GBM")

df <- as.data.frame(read_excel("ELISPOT_Data.xlsx", sheet = "BrMET008"))
col_order <- colnames(df)

df.gather <- gather(df, key = "Genes", value = "Spots")

df.gather <- filter(df.gather, !is.na(Spots))

df2 <- data_summary(df.gather, varname = "Spots", groupnames = "Genes")
df2$Genes <- factor(df2$Genes, levels = col_order)
df2 <- df2[order(df2$Genes),]
df2$color <- c(rep("grey30", length(col_order)-1), "red")

cols <- c("grey30" = "grey30", "red" = "red")

ggplot(df2, aes(x = Genes, y = Spots, fill=color)) + geom_bar(stat = 'identity', color = 'black') + geom_errorbar(aes(ymin=Spots, ymax=Spots+sd), width=.2, position=position_dodge(.9)) + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=14), legend.position = "none" , axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + ylab('Spots per 20,000 TIL') + ggtitle("BrMET008 TIL Class I") + scale_fill_manual(values = cols) + scale_y_continuous(expand = c(0,0), limits = c(0, 45))

