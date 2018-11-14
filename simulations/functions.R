# Functions to generate the simulated data.

library(DropletUtils)

SIMFUN <- function(raw.mat, group1=500, group2=500, reorder.rate=0.1, down.rate=0.1, lower=100) {
    stats <- barcodeRanks(raw.mat)

    # Assuming all cells below the inflection point are empty droplets.
    totals <- stats$total
    empty.limit <- stats$inflection
    ambient.prof <- rowSums(raw.mat[,totals <= empty.limit])
    empty.totals <- totals[totals <= empty.limit]

    # Scrambling to generate resampled empty droplets.
    gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
    gene.ids <- sample(gene.ids)
    cell.ids <- rep(seq_along(empty.totals), empty.totals)
    resampled <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat))

    # Sampling real cells in group 1.
    # All things above the inflection point are assumed to be real.
    # We shuffle 10% of the genes to break any remaining "ambience".
    is.real <- which(empty.limit <= totals)
    mat1 <- raw.mat[,sample(is.real, group1, replace=TRUE)]
    mat1 <- reorderRows(mat1, round(nrow(mat1) * reorder.rate))

    # Sampling cells in group 2, with downsampling.
    mat2 <- raw.mat[,sample(is.real, group2, replace=TRUE)]
    mat2 <- reorderRows(mat2, round(nrow(mat2) * reorder.rate))
    mat2 <- downsampleMatrix(mat2, prop=down.rate)
    
    # Completed.
    return(list(counts=cbind(resampled, mat1, mat2),
                identity=rep(0:2, c(ncol(resampled), ncol(mat1), ncol(mat2)))))
}

reorderRows <- function(mat, N) {
    i <- seq_len(nrow(mat))
    choice <- sample(nrow(mat), N)
    i[choice] <- sample(choice)
    mat[i,,drop=FALSE]
}

assessMethod <- function(keep, identity) {
    g1 <- sum(keep & identity==1, na.rm=TRUE)/sum(identity==1)
    g2 <- sum(keep & identity==2, na.rm=TRUE)/sum(identity==2)
    fdr <- sum(identity==0 & keep, na.rm=TRUE)/sum(keep, na.rm=TRUE)
    return(c(G1=g1, G2=g2, FDR=fdr))
}

createRocPts <- function(stats, identity, n=50) {
    o <- order(stats)
    stats <- stats[o]
    identity <- factor(identity[o])

    chunked <- cut(stats, n)
    by.chunk <- split(identity, chunked)
    N <- do.call(rbind, lapply(by.chunk, table))
    
    apply(N, 2, cumsum)
}

plotBarcodes <- function(ranks, totals, fitted=NULL, subset=NULL, ...) {
    xlim <- range(ranks[ranks>0])
    ylim <- range(totals[totals>0])

    # Dropping non-unique points to save space.
    # Subsetting performed after range() to create comparable plots, if desired.
    keep <- !duplicated(totals)
    if (!is.null(subset)) {
        alt.keep <- keep
        keep[] <- FALSE
        keep[subset] <- alt.keep[subset] 
    }
    Rx <- ranks[keep]
    Tx <- totals[keep]

    # Actually making the plot, also plotting the fitted values if requested.
    plot(Rx, Tx, log="xy", xlab="Rank", ylab="Total count", xlim=xlim, ylim=ylim, ...)
    if (!is.null(fitted)) {  
        Fx <- fitted[keep]
        o <- order(Rx)
        lines(Rx[o], Fx[o], col="red", lwd=2)
    }
    return(invisible(NULL))
}

plotHistogramOutline <- function(breaks, heights, ...) {
    x <- rep(breaks, each=2)
    y <- c(0, rep(heights, each=2), 0)
    lines(x, y, ...)
}

saveRaster <- function(...) {
    tmp <- tempfile(fileext=".png")
    png(tmp, width=7, height=7, units="in", res=300, pointsize=12)
    par(mar=c(0,0,0,0))
    options(bitmapType="cairo")
    plot(..., xlab="", ylab="", axes=FALSE)
    dev.off()
    tmp
}

loadRaster <- function(fpath, log=FALSE) {
    limits <- par("usr")
    img <- png::readPNG(fpath)
    if (log) { limits <- 10^limits }
    rasterImage(img, limits[1], limits[3], limits[2], limits[4])
    invisible(limits)
}

colors <- c("EmptyDrops"="salmon", "CellRanger"="dodgerblue", "Knee point"="orange")
