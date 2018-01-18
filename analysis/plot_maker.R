# This generates pretty plots of the t-SNEs.

library(scran)
sce <- readRDS("sce.rds")
coords <- reducedDim(sce, "TSNE")
dir.create("pics", showWarning=FALSE)

# Defining arrow coordinates.
platelets <- colMeans(coords[sce$Cluster=="10",]) 
lowqual <- colMeans(coords[sce$Cluster=="5",]) 
FUN <- function(coloration, ...) {
    plot(coords[,1], coords[,2], col=coloration, pch=16, xlab="t-SNE1", ylab="t-SNE2", cex.axis=1.2, cex.lab=1.4, ...)
    arrows(platelets[1] - 5, platelets[2], platelets[1] - 1, angle=20, length=0.1, lwd=2)
    arrows(lowqual[1] - 7, lowqual[2], lowqual[1] - 3, angle=20, length=0.1, lwd=2)
}

COLBAR <- function(FUN) {
    x <- max(coords[,1]) + diff(range(coords[,1])) * 0.1
    y <- 1:100/5
    y <- y - mean(y)
    rect(x, y, x+5, y+diff(y)[1]*1.5, col=FUN(length(y)), border=NA)
    text(x+2.5, max(y), "High", pos=3)
    text(x+2.5, min(y), "Low", pos=1)
}

# Making a plot of detection status.
coloration <- rep("grey", nrow(coords))
coloration[sce$Detection=="EmptyDrops"] <- "salmon"
coloration[sce$Detection=="CellRanger"] <- "dodgerblue"

pdf("pics/by_detection.pdf")
FUN(coloration)
legend("bottomright", legend=c("Both", "EmptyDrops", "CellRanger"), col=c("grey", "salmon", "dodgerblue"), pch=16)
dev.off()

# Making a plot of PF4 and PPBP expression.
library(viridis)

by.gene <- Matrix::colSums(logcounts(sce)[c("GP9", "PPBP", "PF4"),])
by.segment <- cut(by.gene, 100)
coloration <- viridis(100)[by.segment]

pdf("pics/by_platelet.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Platelet gene expression", cex.main=1.4)
COLBAR(viridis)
dev.off()

# Making a plot of ribosomal gene expression.
all.ribo <- grep("^RP(L|S)[0-9]+", rownames(sce))
by.gene <- Matrix::colSums(logcounts(sce)[all.ribo,])
by.segment <- cut(by.gene, 100)
coloration <- inferno(100)[by.segment]

pdf("pics/by_ribo.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Ribosomal protein expression", cex.main=1.4)
COLBAR(inferno)
dev.off()

# Making a plot of mitochondrial gene expression.
all.mito <- grep("^MT-", rownames(sce))
by.gene <- Matrix::colSums(counts(sce)[all.mito,])/Matrix::colSums(counts(sce))
by.segment <- cut(by.gene, 100)
coloration <- magma(100)[by.segment]

pdf("pics/by_mito.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Mitochondrial proportion", cex.main=1.4)
COLBAR(magma)
dev.off()


