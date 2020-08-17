library(clonevol)
library(fishplot)

setwd("~/Desktop/GBM/GBM065/sciclone")


driver_list <- list("chr17:7675994", "chr17:7673782", "chr2:208248388", "chr2:47806213")
names(driver_list) <- c("TP53-1", "TP53-2", "IDH1", "MSH6")

my_data <- read.table("mindepth_50_no_cn/GBM065_clusters.tsv", header=T)
merge(my_data, c("GBM065.re1", "GBM065.re3", "GBM065.re4", "GBM065.re5"), "GBM065.Re")
# shorten vaf column names as they will be
vaf.col.names <- grep('.vaf', colnames(my_data), value=T) 
sample.names <- gsub('.vaf', '', vaf.col.names)
my_data[, sample.names] <- my_data[, vaf.col.names] 
vaf.col.names <- sample.names
# prepare sample grouping
sample.groups <- c('P', 'R', 'R', 'R', 'R')
names(sample.groups) <- vaf.col.names
# setup the order of clusters to display in various plots (later)
my_data_tidy <- my_data[order(my_data$cluster),]

# merge samples into 1 sample
#my_data_tidy <- merge.samples(my_data_tidy, c("GBM065.re1", "GBM065.re3")), "GBM065.Re", "Re",
      ref.cols = NULL, var.cols = NULL)

# omit rows that didn't pass depth filter (cluster == NA)
my_data_tidy <- na.omit(my_data_tidy)
# create variant column to identify drivers
my_data <- unite(my_data, variant, "chr", "st", sep=':')

#my_data_tidy$gene <- '-'
#my_data_tidy$is.driver <- FALSE

# clone colors
clone.colors <- c('#999793', '#8d4891', '#f8e356', '#fe9536', '#d7352e', "deepskyblue")

pdf('box.pdf', width = 8, height = 8, useDingbats = FALSE, title='') 
pp <- plot.variant.clusters(my_data_tidy,
     cluster.col.name = 'cluster',
     show.cluster.size = FALSE,
     cluster.size.text.color = 'blue',
     vaf.col.names = vaf.col.names,
     vaf.limits = 70,
     sample.title.size = 10,
     violin = FALSE,
     box = FALSE,
     jitter = TRUE,
     jitter.shape = 1,
     jitter.color = clone.colors,
     jitter.size = 3,
     jitter.alpha = 1,
     jitter.center.method = 'median',
     jitter.center.size = 1,
     jitter.center.color = 'darkgray',
     jitter.center.display.value = 'none',
     #highlight = 'is.driver',
     #highlight.shape = 21,
     #highlight.color = 'blue',
     #highlight.fill.color = 'green',
     #highlight.note.col.name = 'gene',
     #highlight.note.size = 2,
     order.by.total.vaf = FALSE)
dev.off()

plot.pairwise(my_data_tidy, col.names = vaf.col.names, out.prefix = 'variants.pairwise.plot',
              colors = clone.colors)

pdf('flow.pdf', width=8, height=8, useDingbats=FALSE, title='')
plot.cluster.flow(my_data_tidy, vaf.col.names = vaf.col.names, sample.names = c('Primary', 'R1', 'R2', 'R3', 'R5'), colors = clone.colors)
dev.off()

y = infer.clonal.models(variants = my_data_tidy, cluster.col.name = 'cluster', vaf.col.names = vaf.col.names, sample.groups = sample.groups, cancer.initiation.model='monoclonal', subclonal.test = 'bootstrap', subclonal.test.model = 'non-parametric', num.boots = 1000,
    founding.cluster = 1,
    cluster.center = 'mean',
    ignore.clusters = NULL,
    clone.colors = clone.colors,
    min.cluster.vaf = 0.01,
    # min probability that CCF(clone) is non-negative
    sum.p = 0.05,
    # alpha level in confidence interval estimate for CCF(clone) 
    alpha = 0.05)

y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')

plot.clonal.models(y,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver', fancy.variant.boxplot.highlight.shape = 21, fancy.variant.boxplot.highlight.fill.color = 'red', fancy.variant.boxplot.highlight.color = 'black', fancy.variant.boxplot.highlight.note.col.name = 'gene', fancy.variant.boxplot.highlight.note.color = 'blue', fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25, cell.border.color = 'black', clone.grouping = 'horizontal', 
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE, show.score = FALSE,
                   cell.frac.ci = TRUE, disable.cell.frac = FALSE, 
                   # output figure parameters
                   out.dir = 'output', out.format = 'pdf', overwrite.output = TRUE,
                   width = 8,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,5,3))




pdf('trees.pdf', width = 5, height = 5, useDingbats = FALSE) 
plot.all.trees.clone.as.branch(y, branch.width = 0.5, node.size = 1, node.label.size = 0.25)
dev.off()


f <- generateFishplotInputs(y, rescale = TRUE, samples = NULL)
fishes = createFishPlotObjects(f)
pdf("fishplot.pdf", width=8, height=5)
for (i in 1:length(fishes)){
    fish = layoutClones(fishes[[i]])
    fish = setCol(fish,f$clonevol.clone.colors)
    fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
             vlines=seq(1, length(sample.names)), vlab=sample.names, pad.left=0.5)
    
}
dev.off()
