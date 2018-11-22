# This generates pretty plots of the t-SNEs.

library(scran)
sce <- readRDS("sce.rds")
coords <- reducedDim(sce, "TSNE")
dir.create("pics", showWarning=FALSE)

mo.clust <- metadata(sce)$Monocytes
monocytes <- colMeans(coords[sce$Cluster==mo.clust,]) 

# Defining arrow coordinates.
FUN <- function(coloration, ...) {
  plot(coords[,1], coords[,2], col=coloration, pch=16, xlab="t-SNE1", ylab="t-SNE2", cex.axis=1.2, cex.lab=1.4, ...)
  SHIFT <- c(11.5, -4)
  WIDTH <- c(0, 4)
  arrows(monocytes[1] - SHIFT[1] - WIDTH[1], monocytes[2] - SHIFT[2],
         monocytes[1] - SHIFT[1]+5, monocytes[2] - SHIFT[2] - WIDTH[2],
         angle=20, length=0.1, lwd=2, col="black")
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
legend("bottomright", legend=c("Both", "EmptyDrops", "CellRanger"), col=c("grey", "salmon", "dodgerblue"), pch=16)
dev.off()

# Making a plot of KCNA5, STX11, CFP, and S100A12 expression.
library(viridis)

by.gene <- Matrix::colSums(logcounts(sce)[c("KCNA5", "STX11", "CFP", "S100A12"),])
by.segment <- cut(by.gene, 100)
coloration <- viridis(100)[by.segment]

pdf("pics/by_monocyte.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Monocyte gene expression", cex.main=1.4)
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
by.gene <- sce$pct_counts_Mito
by.gene <- pmin(by.gene, 10)
by.segment <- cut(by.gene, 100)
coloration <- plasma(100)[by.segment]

pdf("pics/by_mito.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Mitochondrial proportion", cex.main=1.4)
COLBAR(plasma, high.text=expression("" >= "10%"), low.text="0%")
dev.off()