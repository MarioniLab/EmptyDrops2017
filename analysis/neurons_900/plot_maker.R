# This generates pretty plots of the t-SNEs.

library(scran)
sce <- readRDS("sce.rds")
coords <- reducedDim(sce, "TSNE")
dir.create("pics", showWarning=FALSE)

# Defining arrow coordinates.
interneuron <- colMeans(coords[sce$Cluster=="9",]) 
FUN <- function(coloration, ...) {
    plot(coords[,1], coords[,2], col=coloration, pch=16, xlab="t-SNE1", ylab="t-SNE2", cex.axis=1.2, cex.lab=1.4, ...)
    SHIFT <- c(0, -16)
    WIDTH <- c(0, 4)
    arrows(interneuron[1] - SHIFT[1] - WIDTH[1], interneuron[2] - SHIFT[2], 
           interneuron[1] - SHIFT[1], interneuron[2] - SHIFT[2] - WIDTH[2], angle=20, length=0.1, lwd=2)
}

COLBAR <- function(FUN, high.text="High", low.text="Low") {
    x <- max(coords[,1]) + diff(range(coords[,1])) * 0.1
    y <- 1:100/5
    y <- y - mean(y)
    rect(x, y, x+5, y+diff(y)[1]*1.5, col=FUN(length(y)), border=NA)
    text(x+2.5, max(y), high.text, pos=3)
    text(x+2.5, min(y), low.text, pos=1)
}

# Making a plot of detection status.
coloration <- rep("grey", nrow(coords))
coloration[sce$Detection=="EmptyDrops"] <- "salmon"
coloration[sce$Detection=="CellRanger"] <- "dodgerblue"

pdf("pics/by_detection.pdf")
FUN(coloration)
legend("bottomright", legend=c("Both", "EmptyDrops"), col=c("grey", "salmon"), pch=16)
dev.off()

# Making a plot of interneuron marker expression.
library(viridis)

by.gene <- Matrix::colSums(logcounts(sce)[c("Gad1", "Gad2", "Slc6a1"),])
by.gene <- pmin(by.gene, quantile(by.gene, 0.99)) # just to bring out the contrast a bit more.
by.segment <- cut(by.gene, 100)
coloration <- viridis(100)[by.segment]

pdf("pics/by_interneuron.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Interneuron gene expression", cex.main=1.4)
COLBAR(viridis)
dev.off()

# Making a plot of ribosomal gene expression.
all.ribo <- grep("^Rp(l|s)[0-9]+", rownames(sce))
by.gene <- Matrix::colSums(logcounts(sce)[all.ribo,])
by.segment <- cut(by.gene, 100)
coloration <- inferno(100)[by.segment]

pdf("pics/by_ribo.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Ribosomal protein expression", cex.main=1.4)
COLBAR(inferno)
dev.off()

# Making a plot of mitochondrial gene expression.
by.gene <- sce$pct_counts_Mito
by.gene <- pmin(by.gene, 20)
by.segment <- cut(by.gene, 100)
coloration <- plasma(100)[by.segment]

pdf("pics/by_mito.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Mitochondrial proportion", cex.main=1.4)
COLBAR(magma, high.text=expression("" >= "20%"), low.text="0%")
dev.off()


