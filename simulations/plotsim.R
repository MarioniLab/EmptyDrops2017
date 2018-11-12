# This makes plots of the simulation results for all metrics.
    
source("functions.R")
resdir <- "results-sim"
newdir <- file.path(resdir, "pics")
dir.create(newdir, showWarning=FALSE)

for (f in list.files(resdir, pattern=".*-res\\.tsv$", full=TRUE)) { 
    tab <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    tab <- tab[tab$Method %in% c("emptyDrops", "CellRanger"),]
    stub <- sub("(.+)-res.tsv", "\\1", basename(f))
    scenario <- paste0(tab$G1Size, "/", tab$G2Size)

    combined.recall <- c(tab$G1, tab$G2)
    combined.scenario <- c(paste0(scenario, "_G1_", tab$Method), paste0(scenario, "_G2_", tab$Method))
    by.scenario <- split(combined.recall, combined.scenario)

    # Creating a bar plot, with points.
    scenarios <-  sub("_.*", "", names(by.scenario))
    methods <- sub(".*_", "", names(by.scenario))
    poptype <- sub(".*_(.*)_.*", "\\1", names(by.scenario))
    Xcoords <- seq_along(by.scenario) + cumsum(as.integer(methods==methods[1]))
    ylim <- range(combined.recall)

    pdf(file.path(newdir, paste0(stub, "-recall.pdf")), width=8, height=6)
    par(mar=c(5.1, 4.1, 2.1, 11.1))
    CONSTANT <- 0.4
    plot(0,0,type="n", xlim=range(Xcoords) + c(-CONSTANT, CONSTANT)*2, ylim=ylim, 
        ylab="Recall", xlab="", xaxt="n", cex.axis=1.2, cex.lab=1.4, xaxs="i")

    # Creating the other axis.
    Xgroup <- split(Xcoords, scenarios)
    left <- unlist(lapply(Xgroup, min)) - CONSTANT
    right <- unlist(lapply(Xgroup, max)) + CONSTANT
    middle <- unlist(lapply(Xgroup, mean))
    Yline <- ylim[1] - diff(ylim) * 0.07
    for (g in names(Xgroup)) { 
        segments(left[[g]], Yline, right[[g]], lwd=2, xpd=TRUE)
        text(middle[[g]], Yline, labels=g, xpd=TRUE, pos=1, cex=1.2)
        rect(left[[g]]-CONSTANT, ylim[1]-1, right[[g]]+CONSTANT, ylim[2]*2, border=NA, col="grey95")
    }
    box()

    for (i in seq_along(by.scenario)) {
        current <- by.scenario[[i]]
        curcol <- colors[methods[i]]
        pch <- ifelse(poptype[i]=="G1", 16, 4)
        curX <- Xcoords[i]
        points(rep(curX, length(current)), current, col=curcol, pch=pch)

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

    par(xpd=TRUE)
    legend(max(Xcoords)+CONSTANT*3, max(combined.recall),
        col=rep(colors[c("CellRanger", "emptyDrops")], 2),
        pch=rep(c(16, 4), each=2),
        legend=c("CellRanger (large)", "emptyDrops (large)", "CellRanger (small)", "emptyDrops (small)"),
        lwd=2)

    dev.off()
}


