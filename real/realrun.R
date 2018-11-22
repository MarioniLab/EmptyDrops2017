library(DropletUtils)
library(Matrix)
library(UpSetR)
source("../simulations/functions.R")
ppath <- "pics-real"
dir.create(ppath)

# Specifying all the input files.
ALLFILES <- c("pbmc4k/raw_gene_bc_matrices/GRCh38",
              "neurons_900/raw_gene_bc_matrices/mm10/",
              "293t/matrices_mex/hg19/",
              "jurkat/matrices_mex/hg19/",
              "t_4k/raw_gene_bc_matrices/GRCh38/",
              "neuron_9k/raw_gene_bc_matrices/mm10/",
              sprintf("placenta%i/placenta%i_raw_gene_bc_matrices/GRCh38/", seq_len(6), seq_len(6)))
expected <- c(4000,
              900,
              2800,
              3200,
              4000,
              9000,
              rep(5000, 6))

for (i in seq_along(ALLFILES)) { 
    fname <- ALLFILES[i]
    fpath <- file.path("..", "data", fname)
    if (!file.exists(fpath)) {
        message("missing data files for '", fname, "'")
        next
    }

    sce <- read10xCounts(fpath)
    stub <- sub("/.*", "", fname)

    final <- counts(sce)
    stats <- barcodeRanks(final)

    # Making a barcode plot for diagnostics.
    pdf(file.path(ppath, paste0("rank_", stub, ".pdf")))
    plotBarcodes(stats$rank, stats$total, fitted=stats$fitted, cex.axis=1.2, cex.lab=1.4, main=stub, col="grey70", pch=16)
    abline(h=stats$knee, col="dodgerblue", lty=2, lwd=2)
    abline(h=stats$inflection, col="forestgreen", lty=2, lwd=2)
    legend("bottomleft", lty=c(1,2,2), col=c("red", "dodgerblue", "forestgreen"),
           legend=c("Fitted spline", "Knee", "Inflection"), lwd=2, cex=1.2)
    dev.off()

    # Testing emptyDrops.
    e.out <- emptyDrops(final)
    e.keep <- e.out$FDR <= 0.001
    e.keep[is.na(e.keep)] <- FALSE
    
    # Keeping everything above the knee point.
    k.keep <- stats$knee <= stats$total

    # Using the CellRanger approach.
    c.keep <- defaultDrops(final, expected=expected[i])

    ############################
    # Quantifying the intersections.
    at.least.one <- k.keep | e.keep | c.keep
    collected <- data.frame(EmptyDrops=as.integer(e.keep), 
                            `Knee point`=as.integer(k.keep), 
                            CellRanger=as.integer(c.keep),
                            check.names=FALSE)[at.least.one,]

    pdf(file.path(ppath, paste0("intersect_", stub, ".pdf")), onefile=FALSE)
    upset(collected, sets.bar.color = "#56B4E9", point.size=5, order.by = "freq",
          text.scale=c(2, 1.5, 1, 1, 1.5, 1.5))
    dev.off()

    ############################
    # Having a look at the distribution of retained cells in more detail.
    limits <- log10(range(stats$total[at.least.one]))
    breaks <- seq(limits[1], limits[2], length.out=21)
    modes <- list("EmptyDrops"=e.keep, "Knee point"=k.keep, "CellRanger"=c.keep)

    collected.x <- collected.y <- list()
    for (mode in names(modes)) {
        keep <- modes[[mode]]
        d <- hist(log10(stats$total[keep]), breaks=breaks, plot=FALSE)
        collected.x[[mode]] <-d$mid
        collected.y[[mode]] <-d$count
    }
    xrange <- range(sapply(collected.x, range))
    yrange <- range(sapply(collected.y, range))
    
    pdf(file.path(ppath, paste0("kept_", stub, ".pdf")))
    plot(0,0,type="n", xlim=xrange, ylim=yrange, xlab=expression(Log[10]~"Total count"), 
         ylab="Number of cells", cex.axis=1.2, cex.lab=1.4, main=stub, cex.main=1.4)

    shift <- 0
    for (mode in names(modes)) { 
        plotHistogramOutline(breaks+shift, collected.y[[mode]], col=colors[mode], lwd=2)
        shift <- shift + 0.01
    }

    legend("topleft", col=colors[names(modes)], legend=names(modes), lwd=2, cex=1.2)
    dev.off()

    ############################
    # Examining the distribution of deviances.
    X <- e.out$Total/1000
    Y <- -e.out$LogProb/1000
    xlim <- range(X[!is.na(e.out$LogProb)]) 
    tmp <- saveRaster(X, Y, pch=16, cex=0.5, col=ifelse(e.keep, "red", "grey80"), log="xy", xlim=xlim)

    pdf(file.path(ppath, paste0("dev_", stub, ".pdf")))
    par(mar=c(5.1, 5.1, 4.1, 1.1))
    
    plot(X, Y, log="xy", 
         xlab=expression("Total count (x "*10^3*")"), 
         ylab=expression("-Log likelihood (x "*10^3*")"), 
         xlim=xlim, type="n", cex.axis=1.2, cex.lab=1.4,
         main=stub, cex.main=1.4)
    loadRaster(tmp, log=TRUE)
    box()

    legend("bottomright", col=c("red", "grey80"), pch=16, cex=1.2,
           legend=c("Putative cell", "Empty droplet"))
    dev.off()

    ############################
    # Creating a MA plot.
    ambient.prof <- rowSums(final[,colSums(final) <= 100])
    cell.prof <- rowSums(final[,e.keep])
    summed <- cbind(ambient.prof, cell.prof)
    vals <- edgeR::cpm(summed, log=TRUE, prior.count=5)
    M <- vals[,1] - vals[,2]
    A <- (vals[,1] + vals[,2])/2

    tmp <- saveRaster(A, M, pch=16, cex=0.5)

    pdf(file.path(ppath, paste0("ma_", stub, ".pdf")))
	par(mar=c(5.1, 4.1, 2.1, 5.1))
    plot(A, M, xlab="A", ylab="M (ambient - cells)", cex.main=1.4, main=stub,
        cex.axis=1.2, cex.lab=1.4, pch=16, cex=0.2, type="n")
    limits <- loadRaster(tmp)

	# Choosing the top 10 DE genes to show, by fitting a Poisson GLM.
    library(edgeR)
    y <- DGEList(summed)
    retained <- which(aveLogCPM(y) > 1)
    y <- y[retained,]
    y <- calcNormFactors(y)
    design <- cbind(1, c(0, 1))
    fit <- glmFit(y, design, dispersion=0) 
    res <- glmLRT(fit)
    o <- order(res$table$PValue)
	extreme <- retained[head(o, 10)]

    ypos <- rank(M[extreme]) / length(extreme) * diff(range(M)) + min(M)
    xpos <- limits[2] + 0.01 * (limits[2] - limits[1])
    text(xpos, ypos, col="red", pos=4, xpd=TRUE, offset=0.1,
        rowData(sce)$Symbol[match(names(M)[extreme], rownames(sce))])
    segments(xpos, ypos, A[extreme], M[extreme], col="red", xpd=TRUE)
    box()
    dev.off()
}
