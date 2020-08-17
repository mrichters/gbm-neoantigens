# 1.16.2020
# reduce the number of viral probes for targeted sequencing
library(readxl)
library(tidyverse)

probes <- read_excel("~/Desktop/targeted_resequencing/IDT_probe_designs/viral_probes/02839966_Fronick_Disco_JT_probes_v2.xlsx", sheet = "Probes")

probes.filter <- filter(probes, Risk == "LOW")

virus_names <- unique(probes.filter$`Group Name`)
df.final <- data.frame()
for ( virus in virus_names ) {
  df.subset <- filter(probes.filter, `Group Name` == virus) 
  #print(virus); print(nrow(df.subset))
  df.subset.order <- df.subset[order(df.subset$Start), ]
  if ( grepl("SV40", virus) ) {
    df.new <- df.subset.order[seq(1, nrow(df.subset.order), 5), ]
  } else {
    df.new <- df.subset.order[seq(1, nrow(df.subset.order), 10), ]
  }
  print(virus); print(nrow(df.new))
  #df.final <- rbind(df.final, df.new)
}

write_tsv(df.final, "final_viral_probes_list.tsv")
