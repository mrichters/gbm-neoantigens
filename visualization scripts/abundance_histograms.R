
#open transcript abundance file
tx_data <- read.table('GBM052-3_abundance.tsv', sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE)
#filter tpm values to remove zeros
filtered_tx <- tx_data %>% filter(tpm > 0.0)
#non-zero transcript abundance distribution
ggplot(filtered_tx, aes(x=tpm)) + geom_histogram(binwidth = 2000) + ggtitle("Transcript Abundance (non-zero values)") + xlab("TPM")
#create log2 transformed graph
filtered_tx[, 5] <- log2(filtered_tx[5])
#log2 abundance distribution
ggplot(filtered_tx, aes(x=tpm)) + geom_histogram(binwidth = 1) + ggtitle("Transcript Abundance") + xlab("log2(TPM)")

#open gene abundance file
gx_data <- read.table('GBM052-3_gene_abundance.tsv', sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE)
#filter tpm values to remove zeros
filtered_gx <- gx_data %>% filter(abundance > 0.0)
#non-zero gene abundance distribution
ggplot(filtered_gx, aes(x=abundance)) + geom_histogram(binwidth = 2000) + ggtitle("Gene Abundance (non-zero value)") + xlab("TPM")
#create log2 transformed graph
filtered_gx[, 3] <- log2(filtered_gx[3])
#log2 abundance distribution
ggplot(filtered_gx, aes(x=abundance)) + geom_histogram(binwidth = 1) + ggtitle("Gene Abundance") + xlab("log2(TPM)")                                             
