library(DropletUtils)
library(Matrix)
source("functions.R")
set.seed(73953024)

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
    sce <- read10xCounts(file.path("..", "data", fname))
    stub <- sub("/.*", "", fname)
    ofile <- file.path(opath, paste0(stub, "-res.tsv"))
    ffile <- file.path(opath, paste0(stub, "-fdr.tsv"))
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
    
                    pdf(file.path(ppath, paste0(stub, "_", g1, "_", g2, "-sim.pdf")), width=10, height=8)
                    par(mfrow=c(2,2), cex.axis=1.2, cex.lab=1.4, cex.main=1.4, mar=c(4.1, 4.1, 2.1, 1.1))
                    plotBarcodes(r, ot, pch=16, main="All barcodes")
                    plotBarcodes(r, ot, pch=16, subset=(oi==0), col="grey80", main="Empty droplets")
                    plotBarcodes(r, ot, pch=16, subset=(oi==1), col="blue", main="Large cells")
                    plotBarcodes(r, ot, pch=16, subset=(oi==2), col="red", main="Small cells")
                    dev.off()
                }

                # Testing emptyDrops.
                e.out <- emptyDrops(final)
                is.sig <- e.out$FDR <= 0.001
                emp.res <- assessMethod(is.sig, out$identity)

                # Keeping everything above the knee point.
                K <- barcodeRanks(final)$knee
                knee.res <- assessMethod(totals >= K, out$identity)

                # Not using the knee point safeguard.
                E2 <- p.adjust(e.out$PValue, method="BH")
                emp2.res <- assessMethod(E2 <= 0.001, out$identity)

                # Using the CellRanger approach.
                expected <- sum(out$identity!=0)
                c.keep <- defaultDrops(final, expected=expected)
                cell.res <- assessMethod(c.keep, out$identity)

                # Saving the default assessment results to file.
                write.table(cbind(G1Size=g1, G2Size=g2, 
                                  Method=c("emptyDrops", "Knee point", "No knee", "CellRanger"),
                                  rbind(emp.res, knee.res, emp2.res, cell.res)),
                            file=ofile, append=!unlinker, col.names=unlinker, 
                            row.names=FALSE, quote=FALSE, sep="\t")

                if (first) {
                    # Creating a ROC plot.
                    emp.roc <- createRocPts(log(e.out$FDR+1e-3), out$identity)
                    lib.roc <- createRocPts(-log(totals+1), out$identity)

                    emp.fdr <- emp.roc[,"0"]/rowSums(emp.roc)
                    emp.tp1 <- emp.roc[,"1"]/g1
                    emp.tp2 <- emp.roc[,"2"]/g2
                    lib.fdr <- lib.roc[,"0"]/rowSums(lib.roc)
                    lib.tp1 <- lib.roc[,"1"]/g1
                    lib.tp2 <- lib.roc[,"2"]/g2

                    emp.col <- colors["emptyDrops"]
                    lib.col <- "grey50"

                    pdf(file.path(ppath, paste0(stub, "_", g1, "_", g2, "-roc.pdf")), width=8, height=6)
                    par(mar=c(5.1, 4.1, 4.1, 12.1))
                    plot(c(0, emp.fdr), c(0, emp.tp1), xlim=c(0, 0.01), ylim=c(0, 1), type="l", col=emp.col, lwd=3,
                        xlab="False discovery rate", ylab="True positive rate", cex.axis=1.2, cex.lab=1.4, cex.main=1.4)
                    lines(c(0, lib.fdr), c(0, lib.tp1), col=lib.col,lwd=3)
                    lines(c(0, emp.fdr), c(0, emp.tp2), col=emp.col, lwd=3, lty=2)
                    lines(c(0, lib.fdr), c(0, lib.tp2), col=lib.col, lwd=3, lty=2)

                    par(xpd=TRUE)
                    legend(0.011, 1, lwd=3, col=rep(c(emp.col, lib.col), each=2), lty=rep(1:2, 2),
                        legend=c("emptyDrops (large)", "emptyDrops (small)", "Total count (large)", "Total count (small)"))
                    dev.off()
                }

                # Reporting on FDR control.
                all.thresholds <- c(0.001, 0.002, 0.005, 0.01) 
                collected <- numeric(length(all.thresholds))
                for (i in seq_along(all.thresholds)) {
                    collected[[i]] <- assessMethod(e.out$FDR <= all.thresholds[i], out$identity)["FDR"]
                }

                collected <- rbind(collected)
                colnames(collected) <- all.thresholds
                write.table(cbind(G1Size=g1, G2Size=g2, collected), file=ffile, append=!unlinker,
                        col.names=unlinker, row.names=FALSE, quote=FALSE, sep="\t")

                unlinker <- FALSE
                first <- FALSE
            }
        }
    }
}
