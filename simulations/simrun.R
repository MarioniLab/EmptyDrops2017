library(DropletUtils)
library(Matrix)
source("functions.R")

opath <- "results-sim"
dir.create(opath)
ppath <- "pics-sim"
dir.create(ppath)

ALLFILES <- c("pbmc4k/raw_gene_bc_matrices/GRCh38",
              "neurons_900/raw_gene_bc_matrices/mm10/",
              "293t/matrices_mex/hg19/",
              "jurkat/matrices_mex/hg19/",
              "t_4k/raw_gene_bc_matrices/GRCh38/",
              "neuron_9k/raw_gene_bc_matrices/mm10/")

for (fname in ALLFILES) { 
    sce <- read10xResults(file.path("..", "data", fname))
    stub <- sub("/.*", "", fname)
    ofile <- file.path(opath, paste0(stub, ".tsv"))
    unlinker <- TRUE 
    
    for (g1 in c(500, 2000)){ 
        for (g2 in c(500, 2000)) {
            first <- TRUE 

            for (it in 1:10) { 
                out <- SIMFUN(counts(sce), group1=g1, group2=g2)
                final <- out$counts
                totals <- colSums(final)
               
                if (first) {  
                    o <- order(totals, decreasing=TRUE)
                    r <- seq_along(o)
                    ot <- totals[o]
                    oi <- out$identity[o]
    
                    pdf(file.path(ppath, paste0(stub, "_", g1, "_", g2, ".pdf")), width=10, height=8)
                    par(mfrow=c(2,2), cex.axis=1.2, cex.lab=1.4, cex.main=1.4, mar=c(4.1, 4.1, 2.1, 1.1))
                    plotBarcodes(r, ot, pch=16, main="All barcodes")
                    plotBarcodes(r, ot, pch=16, subset=(oi==0), col="grey80", main="Empty droplets")
                    plotBarcodes(r, ot, pch=16, subset=(oi==1), col="blue", main="Large cells")
                    plotBarcodes(r, ot, pch=16, subset=(oi==2), col="red", main="Small cells")
                    dev.off()
                    first <- FALSE
                }
                
                # Testing emptyDrops.
                e.out <- emptyDrops(final)
                is.sig <- e.out$FDR <= 0.01
                emp.res <- assessMethod(is.sig, out$identity)
                
                # Keeping everything above the knee point.
                K <- barcodeRanks(final)$knee
                knee.res <- assessMethod(totals >= K, out$identity)
                
                # Using the CellRanger approach.
                expected <- sum(out$identity!=0)
                expected.totals <- sort(totals, decreasing=TRUE)[seq_len(expected)]
                threshold <- quantile(expected.totals, 0.99)/10
                cell.res <- assessMethod(totals >= threshold, out$identity)
    
                # Saving the current set of results to file.
                write.table(cbind(G1Size=g1, G2Size=g2, 
                                  Method=c("emptyDrops", "Knee point", "CellRanger"),
                                  rbind(emp.res, knee.res, cell.res)),
                            file=ofile, append=!unlinker, col.names=unlinker, 
                            row.names=FALSE, quote=FALSE, sep="\t")
                unlinker <- FALSE
            }
        }
    }
}
