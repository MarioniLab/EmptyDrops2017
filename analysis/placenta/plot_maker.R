# This generates pretty plots of the t-SNEs.

library(scran)
sce <- readRDS("sce.rds")
coords <- reducedDim(sce, "TSNE")
dir.create("pics", showWarning=FALSE)

# Defining arrow coordinates.
FUN <- function(coloration, loc, SHIFT, WIDTH, orientation=-5, ...) {
  plot(coords[,1], coords[,2], col=coloration, pch=16, xlab="t-SNE1", ylab="t-SNE2", cex.axis=1.2, cex.lab=1.4, ...)
  arrows(loc[1] + SHIFT[1] + WIDTH[1], loc[2] + SHIFT[2],
         loc[1] + SHIFT[1] + orientation, loc[2] + SHIFT[2] + WIDTH[2],
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

#######################################################
# Monocyte arrows.

mo.clust <- metadata(sce)$Monocytes
monocytes <- colMeans(coords[sce$Cluster==mo.clust,]) 
SHIFT <- c(11.5, 8)
WIDTH <- c(0, -4)

# Making a plot of detection status.
coloration <- rep("grey", nrow(coords))
coloration[sce$Detection=="EmptyDrops"] <- "salmon"
coloration[sce$Detection=="CellRanger"] <- "dodgerblue"

pdf("pics/by_detection.pdf")
FUN(coloration, loc=monocytes, SHIFT=SHIFT, WIDTH=WIDTH)
legend("bottomright", legend=c("Both", "EmptyDrops", "CellRanger"), col=c("grey", "salmon", "dodgerblue"), pch=16)
dev.off()

# Making a plot of Monocyte gene expression.
library(viridis)

by.gene <- sce$monocyte_genes_total 
by.segment <- cut(by.gene, 100)
coloration <- viridis(100)[by.segment]

pdf("pics/by_monocyte.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Monocyte gene expression", cex.main=1.4, loc=monocytes, SHIFT=SHIFT, WIDTH=WIDTH)
COLBAR(viridis)
dev.off()

#######################################################
# T cell arrows.

t.clust <- metadata(sce)$TCells
tcells <- colMeans(coords[sce$Cluster==t.clust,]) 
SHIFT <- c(-10, 4)
WIDTH <- c(0, -4)

by.gene <- sce$t_genes_total 
by.segment <- cut(by.gene, 100)
coloration <- viridis(100)[by.segment]

pdf("pics/by_tcell.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="T cell gene expression", cex.main=1.4, loc=tcells, SHIFT=SHIFT, WIDTH=WIDTH, orientation=5)
COLBAR(viridis)
dev.off()

#######################################################
# Stripped nuclei arrows.

stripped <- metadata(sce)$StrippedNucleus
striploc <- colMeans(coords[sce$Cluster==stripped,]) 
SHIFT <- c(-10, 4)
WIDTH <- c(0, 0)

by.gene <- sce$mito_genes_total 
by.segment <- cut(by.gene, 100)
coloration <- plasma(100)[by.segment]

pdf("pics/by_mito.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Mitochondrial gene expression", cex.main=1.4, loc=striploc, SHIFT=SHIFT, WIDTH=WIDTH, orientation=5)
COLBAR(plasma)
dev.off()

# Making a plot of ribosomal gene expression.
all.ribo <- grep("^RP(L|S)[0-9]+", rownames(sce))
by.gene <- Matrix::colSums(logcounts(sce)[all.ribo,])
by.segment <- cut(by.gene, 100)
coloration <- inferno(100)[by.segment]

pdf("pics/by_ribo.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
FUN(coloration, main="Ribosomal protein expression", cex.main=1.4, loc=striploc, SHIFT=SHIFT, WIDTH=WIDTH)
COLBAR(inferno)
dev.off()


