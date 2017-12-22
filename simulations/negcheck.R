# This is a negative control check to ensure that no barcodes are 
# rejected from what we've assumed to be the ambient pool.

library(DropletUtils)
library(Matrix)
library(viridis)
source("functions.R")
fout <- "pics-negcheck"
dir.create(fout, showWarning=FALSE)

# Defining the files to use.
ALLFILES <- c("pbmc4k/raw_gene_bc_matrices/GRCh38",
              "neurons_900/raw_gene_bc_matrices/mm10/",
              "293t/matrices_mex/hg19/",
              "jurkat/matrices_mex/hg19/",
              "t_4k/raw_gene_bc_matrices/GRCh38/",
              "neuron_9k/raw_gene_bc_matrices/mm10/")

# Looping through the files.
set.seed(1000)
plot.legend <- TRUE
for (fname in ALLFILES) { 
    sce <- read10xResults(file.path("..", "data", fname))
    totals <- colSums(counts(sce))
    ambient <- counts(sce)[,totals<=100 & totals > 0]
    out <- testEmptyDrops(ambient, test.ambient=TRUE)

    # Drawing the standard deviance vs total plot.
    stub <- sub("/.*", "", fname)
    png(file.path(fout, paste0("dev_", stub, ".png")), res=150, width=7, height=7, units='in', pointsize=12)
    plot(out$Total, out$Deviance, log="xy", xlab="Total UMI count", ylab="Deviance", pch=16, cex=0.5, cex.axis=1.2, cex.lab=1.4)
    o <- order(out$Total)
    lines(out$Total[o], out$Expected[o], col="dodgerblue", lwd=2)
    is.sig <- out$FDR <= 0.05
    points(out$Total[is.sig], out$Deviance[is.sig], col="red", pch=16)
    dev.off()

    # Drawing a very complicated histogram.
    pdf(file.path(fout, paste0("hist_", stub, ".pdf")))
#    par(mar=c(5.1, 4.1, 2.1, 4.5))
    N <- 20
    breaks <- seq(0, 1, length.out=N+1)
    choice <- pmax(1, pmin(N, findInterval(out$PValue, breaks)))
    collected <- aggregate(out$Total, by=list(choice), FUN=sum)

    all.percents <- collected[,2]/sum(out$Total)*100
    max.percent <- max(all.percents)
    plot(0, 0, type="n", xlim=c(0,1), ylim=c(0, max.percent), bty="n", main=stub, cex.main=1.4,
         xlab="p-value", ylab="Contribution to ambient profile (%)", cex.axis=1.2, cex.lab=1.4)
    rect(breaks[-N-1], 0, breaks[-1], all.percents, bty="n", col="grey90")
    axis(2, labels=FALSE)

#    minvals <- c(0, 20, 40, 60, 80, 100)
#    N2 <- length(minvals)-1
#    densities <- vector("list", N2)
#    for (m in seq_len(N2)) { 
#        keep <- out$Total > minvals[m] & out$Total <= minvals[m+1]
#        X <- hist(out$PValue[keep], plot=FALSE, breaks=breaks)
#        densities[[m]] <- X$density
#    }
#
#    max.dens <- max(sapply(densities, max))
#    rescale <- max.percent/max.dens * 0.95
#    colors <- viridis(N2)
#    shift <- 0
#    for (m in seq_len(N2)) { 
#        plotHistogramOutline(breaks+shift, densities[[m]]*rescale,  col=colors[m], lwd=2)
#        shift <- shift + 0.002
#    }
#
#    prettified <- pretty(c(0,max.dens))
#    prettified <- prettified[prettified < max.dens]
#    axis(4, at=prettified*rescale, labels=prettified, cex.axis=1.2)
#    mtext("Density", side=4, line=2.5, cex=1.4)
#
#    if (plot.legend){ 
#        legend(1, max.percent, col=colors, lwd=3, 
#               legend=sprintf("(%i, %i]", minvals[-N2-1], minvals[-1]), 
#               cex=1.1, bg="white", xjust=1, yjust=1)
#        plot.legend <- FALSE
#    }
    dev.off()
}
