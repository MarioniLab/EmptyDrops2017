# Functions to generate the simulated data.

library(DropletUtils)

SIMFUN <- function(raw.mat, group1=500, group2=500, down.rate=0.1, lower=100) {
    stats <- barcodeRanks(raw.mat)

    # Assuming all cells below the inflection point are empty droplets.
    totals <- stats$total
    empty.limit <- stats$inflection
    ambient.prof <- rowSums(raw.mat[,totals <= empty.limit])
    empty.totals <- totals[totals <= empty.limit]

    # Scrambling to generate resampled empty droplets.
    gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
    gene.ids <- sample(gene.ids, sum(empty.totals), replace=TRUE)
    cell.ids <- rep(seq_along(empty.totals), empty.totals)
    resampled <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat))

    # Sampling real cells in group 1.
    # All things above the knee point are assumed to be real.
    is.real <- which(empty.limit <= totals)
    mat1 <- raw.mat[,sample(is.real, group1, replace=TRUE)]

    # Sampling cells in group 2, with downsampling.
    mat2 <- downsampleMatrix(raw.mat[,sample(is.real, group2, replace=TRUE)], prop=down.rate)
    
    # Completed.
    return(list(counts=cbind(resampled, mat1, mat2),
                identity=rep(0:2, c(ncol(resampled), ncol(mat1), ncol(mat2)))))
}

assessMethod <- function(keep, identity) {
    g1 <- sum(keep & identity==1, na.rm=TRUE)/sum(identity==1)
    g2 <- sum(keep & identity==2, na.rm=TRUE)/sum(identity==2)
    fdr <- sum(identity==0 & keep, na.rm=TRUE)/sum(keep, na.rm=TRUE)
    return(c(G1=g1, G2=g2, FDR=fdr))
}

plotBarcodes <- function(ranks, totals, fitted=NULL, add=FALSE, ...) {
    # Dropping non-unique plots to save space.
    keep <- !duplicated(totals)
    Rx <- ranks[keep]
    Tx <- totals[keep]

    if (!add) {
        plot(Rx, Tx, log="xy", xlab="Rank", ylab="Total", ...)
    } else {
        points(Rx, Tx, ...)
    }
   
    if (!is.null(fitted)) {  
        Fx <- fitted[keep]
        o <- order(Rx)
        lines(Rx[o], Fx[o], col="red")
    }
    return(invisible(NULL))
}

plotHistogramOutline <- function(breaks, heights, ...) {
    x <- rep(breaks, each=2)
    y <- c(0, rep(heights, each=2), 0)
    lines(x, y, ...)
}

colors <- c("emptyDrops"="salmon", "CellRanger"="dodgerblue", "Knee point"="orange")
