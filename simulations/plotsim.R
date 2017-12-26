# This makes plots of the simulation results for all metrics.
    
source("functions.R")
resdir <- "results-sim"
newdir <- file.path(resdir, "pics")
dir.create(newdir, showWarning=FALSE)

for (f in list.files(resdir, pattern="sim.*\\.tsv$", full=TRUE)) { 
    tab <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    stub <- sub("sim_(.+).tsv", "\\1", basename(f))
    scenario <- paste0(tab$G1Size, "/", tab$G2Size, "_", tab$Method)

    for (s in c("G1", "G2", "FDR")) {
        all.values <- tab[[s]]
        if (s=="G1") {
            ylab <- "Recall (large cells)"
            ylim <- range(all.values)
        } else if (s=="G2") {
            ylab <- "Recall (small cells)"
            ylim <- range(all.values)
        } else {
            ylab <- s
            ylim <- c(0, 0.01)
        }
        by.scenario <- split(all.values, scenario)

        # Creating a bar plot, with points.
        scenarios <-  sub("_.*", "", names(by.scenario))
        methods <- sub(".*_", "", names(by.scenario))
        Xcoords <- seq_along(by.scenario) + cumsum(as.integer(methods==methods[1])*2)
        pdf(file.path(newdir, paste0(stub, "_", s, ".pdf")))
        plot(0,0,type="n", xlim=range(Xcoords), ylim=ylim, ylab=ylab, xlab="", xaxt="n", cex.axis=1.2, cex.lab=1.4)

        # Creating the other axis.
        Xgroup <- split(Xcoords, scenarios)
        left <- unlist(lapply(Xgroup, min)) - 0.4
        right <- unlist(lapply(Xgroup, max)) + 0.4
        middle <- unlist(lapply(Xgroup, mean))
        Yline <- ylim[1] - diff(ylim) * 0.07
        for (g in names(Xgroup)) { 
            segments(left[[g]], Yline, right[[g]], lwd=2, xpd=TRUE)
            text(middle[[g]], Yline, labels=g, xpd=TRUE, pos=1, cex=1.2)
            rect(left[[g]]-0.5, ylim[1]*0.5, right[[g]]+0.5, ylim[2]*2, border=NA, col="grey95")
        }
        box()

        for (i in seq_along(by.scenario)) {
            current <- by.scenario[[i]]
            curcol <- colors[methods[i]]
            curX <- Xcoords[i]
            points(rep(curX, length(current)), current, col=curcol)

            # Specifying the mean.
            mean.val <- mean(current)
            se.val <- sd(current)/sqrt(length(current))
            segments(curX-0.4, mean.val, curX+0.4, col=curcol, lwd=3)

            # Adding error bars.
            lower <- mean.val - se.val
            upper <- mean.val + se.val
            segments(curX, lower, curX, upper, col=curcol, lwd=2)
            segments(curX-0.25, lower, curX+0.25, col=curcol, lwd=2)
            segments(curX-0.25, upper, curX+0.25, col=curcol, lwd=2)
        }

        dev.off()
    }
}

pdf(file.path(newdir, "legend.pdf"))
plot.new()
re.colors <- colors[order(names(colors))]
legend("topleft", col=re.colors, legend=names(re.colors), cex=1.2, lwd=2, pch=1)
dev.off()

