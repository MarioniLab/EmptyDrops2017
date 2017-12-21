# This is a negative control check to ensure that no barcodes are 
# rejected from what we've assumed to be the ambient pool.

ALLFILES <- c("pbmc4k/raw_gene_bc_matrices/GRCh38",
              "neurons_900/raw_gene_bc_matrices/mm10/",
              "293t/matrices_mex/hg19/",
              "ercc/matrices_mex/ercc92/")

library(DropletUtils)
library(Matrix)
library(viridis)
fout <- "pics-negcheck"
dir.create(fout, showWarning=FALSE)
set.seed(1000)

for (fname in ALLFILES) { 
    sce <- read10xResults(file.path("..", "data", fname))
    totals <- colSums(counts(sce))
    ambient <- counts(sce)[,totals<=100 & totals > 0]
    out <- testEmptyDrops(ambient, test.ambient=TRUE)

    stub <- sub("/.*", "", fname)
    png(file.path(fout, paste0("dev_", stub, ".png")), res=150, width=7, height=7, units='in', pointsize=12)
    plot(out$Total, out$Deviance, log="xy", xlab="Total UMI count", ylab="Deviance", pch=16, cex=0.5, cex.axis=1.2, cex.lab=1.4)
    o <- order(out$Total)
    lines(out$Total[o], out$Expected[o], col="dodgerblue", lwd=2)
    is.sig <- out$FDR <= 0.05
    points(out$Total[is.sig], out$Deviance[is.sig], col="red", pch=16)
    dev.off()

    pdf(file.path(fout, paste0("hist_", stub, ".pdf")))
    plot(0,0,type="n", xlim=c(0,1), ylim=c(0, 3), xlab="p-value", ylab="Density", cex.axis=1.2, cex.lab=1.4)
    minvals <- c(0, 10, 20, 40, 60, 80)
    colors <- viridis(length(minvals))
    for (m in seq_along(minvals)) { 
        X <- hist(out$PValue[out$Total >= minvals[m]], plot=FALSE, breaks=20)
        lines(X$mids, X$density, col=colors[m], lwd=3) 
    }
    abline(h=1, col="red", lty=2)
    legend("topright", col=colors, lwd=3, legend=paste("At least", minvals), cex=1.1, bg="white")
    dev.off()
}
