#files <- list.files(path="coverage_files", pattern="*coverage.tsv", full.names=TRUE, recursive=FALSE)
#for (chr in files){

histogram <- read.table("coverage_files/chr21.coverage.tsv", sep="\t", quote="", stringsAsFactors=FALSE, header=FALSE)

colnames(histogram) <- c('Chr', 'Start', 'End', 'Coverage')

#highcov <- histogram %>% filter(Coverage > 100000)
#highcov <- histogram %>% filter(Start > 52995300 & Start < 52997000)
sorted_hg <- histogram[ order(histogram[,4]), ] 
highcov <- histogram %>% filter(Start > 8200000 & Start < 82500000)
sum_hg_14 <- sum(highcov$Coverage)
sum_hg_14
#highcov <- histogram %>% top_n(1, Coverage)


#print(highcov)
#}
plot(histogram$Start, histogram$Coverage)
library(lattice, pos=10)
pdf(file="chr14_highcov.pdf",width=11,height=8)
xyplot(Coverage ~ Start, type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), scales=list(x=list(relation='same'), y=list(relation='same')), data=highcov, main="Chr14 High Coverage - BrMET008")
dev.off()

