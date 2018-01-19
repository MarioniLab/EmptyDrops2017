---
title: Analyzing the 4K PBMC dataset
author: Aaron Lun and others
date: 28 December 2017
output: 
  BiocStyle::html_document:
    fig_caption: no
---



# Introduction

This document performs a brief analysis of the 4K PBMC dataset, with particular focus on the cells that are uniquely detected by EmptyDrops or CellRanger.
First, we load in the raw count matrix:


```r
library(DropletUtils)
fname <- "../../data/pbmc4k/raw_gene_bc_matrices/GRCh38"
sce <- read10xCounts(fname, col.names=TRUE)
sce
```

```
## class: SingleCellExperiment 
## dim: 33694 737280 
## metadata(0):
## assays(1): counts
## rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
##   ENSG00000268674
## rowData names(2): ID Symbol
## colnames(737280): AAACCTGAGAAACCAT-1 AAACCTGAGAAACCGC-1 ...
##   TTTGTCATCTTTAGTC-1 TTTGTCATCTTTCCTC-1
## colData names(2): Sample Barcode
## reducedDimNames(0):
## spikeNames(0):
```

... and define the cells with `emptyDrops`:


```r
set.seed(100)
e.out <- emptyDrops(counts(sce))
e.keep <- e.out$FDR <= 0.01
summary(e.keep)
```

```
##    Mode   FALSE    TRUE    NA's 
## logical     836    4453  731991
```

... and CellRanger (with an expected 4000 cells, as the name of the dataset suggests):


```r
c.keep <- defaultDrops(counts(sce), expected=4000)
summary(c.keep)
```

```
##    Mode   FALSE    TRUE 
## logical  732945    4335
```

We retain all cells that are detected by either EmptyDrops or CellRanger.
This makes it easier to compare the two methods later, rather than having to perform two separate analyses.


```r
keep <- c.keep | (e.keep & !is.na(e.keep))
detection <- rep("Both", length(keep))
detection[c.keep & !e.keep] <- "CellRanger"
detection[!c.keep & e.keep] <- "EmptyDrops"
table(detection)
```

```
## detection
##       Both CellRanger EmptyDrops 
##     736540        311        429
```

Storing this in the metadata for future use.


```r
sce$Detection <- detection
sce$PValue <- e.out$PValue
sce <- sce[,keep]
```

# Adding gene-level annotation

We add some gene-level annotation.


```r
library(EnsDb.Hsapiens.v86)
symb <- mapIds(EnsDb.Hsapiens.v86, keys=rownames(sce), keytype="GENEID", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))
```

```
## DataFrame with 6 rows and 4 columns
##                ID        Symbol         ENSEMBL        SYMBOL
##       <character>   <character>     <character>   <character>
## 1 ENSG00000243485  RP11-34P13.3 ENSG00000243485     MIR1302-2
## 2 ENSG00000237613       FAM138A ENSG00000237613       FAM138A
## 3 ENSG00000186092         OR4F5 ENSG00000186092         OR4F5
## 4 ENSG00000238009  RP11-34P13.7 ENSG00000238009  RP11-34P13.7
## 5 ENSG00000239945  RP11-34P13.8 ENSG00000239945  RP11-34P13.8
## 6 ENSG00000239906 RP11-34P13.14 ENSG00000239906 RP11-34P13.14
```

We relabel the rows with the gene symbols for easier reading.


```r
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)
head(rownames(sce))
```

```
## [1] "MIR1302-2"     "FAM138A"       "OR4F5"         "RP11-34P13.7" 
## [5] "RP11-34P13.8"  "RP11-34P13.14"
```

We also determine the chromosomal location for each gene.


```r
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ENSEMBL, 
    column="SEQNAME", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="MT")
```

```
##    Mode   FALSE    TRUE    NA's 
## logical   33537      13     144
```

# Quality control on the cells

Cell detection can be considered an implicit quality control step, so no extra steps are needed.
Nonetheless, we examine some commonly used metrics.


```r
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=="MT")))
par(mfrow=c(1,3))
hist(log10(sce$total_counts), breaks=20, col="grey80")
hist(log10(sce$total_features_by_counts), breaks=20, col="grey80")
hist(sce$pct_counts_Mito, breaks=20, col="grey80")
```

<img src="analysis_files/figure-html/qchist-1.png" width="100%"  class="widefigure" />

Interestingly, a large number of the features with low total counts also have high mitochondrial proportions.


```r
par(mfrow=c(1,2))
plot(sce$total_features_by_counts, sce$pct_counts_Mito)
plot(sce$total_counts, sce$pct_counts_Mito)
```

<img src="analysis_files/figure-html/qcscatter-1.png" width="100%"  class="widefigure" />

This may indicate that the cells uniquely detected by EmptyDrops are, in fact, damaged.
We'll have a look at this in more detail during the clustering step.

# Examining gene expression

We have a look at the average expression of each gene.


```r
ave <- calcAverage(sce)
rowData(sce)$AveCount <- ave
hist(log10(ave), col="grey80")
```

<img src="analysis_files/figure-html/abhist-1.png" width="100%" />

We also examine the top-most expressed genes.
This contains ribosomal protein genes and other usual suspects.


```r
plotHighestExprs(sce)
```

<img src="analysis_files/figure-html/highexpr-1.png" width="100%"  class="widefigure" />

# Normalizing for cell-specific biases

All cells with outlier values for the library size are defined as one cluster.
This is necessary to avoid problems when normalizing very small libraries with the large libraries.
Everything else is used in clustering with `quickCluster`.


```r
library(scran)
low.lib <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
clusters <- numeric(length(low.lib))
clusters[!low.lib] <- quickCluster(sce[,!low.lib], method="igraph",
    subset.row=ave>=0.1, irlba.args=list(maxit=1000)) # for convergence.
table(clusters)
```

```
## clusters
##    0    1    2    3    4 
##  363  203  608 2488 1102
```

We then use the deconvolution method to compute size factors for each cell.


```r
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
summary(sizeFactors(sce))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  0.007114  0.682433  0.909410  1.000000  1.171121 13.484090
```

We can plot these against the library sizes to see how much of a difference it makes.


```r
plot(sce$total_counts, sizeFactors(sce), log="xy")
```

<img src="analysis_files/figure-html/sfplot-1.png" width="100%" />

Finally, we compute normalized log-expresion values.


```r
sce <- normalize(sce)
```

# Modelling the mean-variance trend

We assume that the technical noise is Poisson and create a fitted trend on that basis.


```r
means <- c(0:10, seq(11, max(ave), length.out=20))
tol <- 1e-8
collected.means <- collected.vars <- numeric(length(means))
for (i in seq_along(means)) {
    m <- means[i]
    lower <- qpois(tol, lambda=m)
    upper <- qpois(tol, lambda=m, lower=FALSE)
    ranged <- lower:upper
    p <- dpois(ranged, lambda=m)
    lvals <- log2(ranged + 1)
    lmean <- sum(lvals * p) / sum(p)
    collected.means[i] <- lmean
    collected.vars[i] <- sum((lvals - lmean)^2 * p) / sum(p)
}
new.trend <- splinefun(collected.means, collected.vars)
```

We actually estimate the variances and plot the trend against the original variances as well.


```r
fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))
plot(fit$mean, fit$var, pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)
curve(new.trend(x), col="red", add=TRUE)
```

<img src="analysis_files/figure-html/trendplot-1.png" width="100%" />

We decompose the variance and have a look at the genes with the highest residual.


```r
dec <- decomposeVar(fit=fit)
top.dec <- dec[order(dec$bio, decreasing=TRUE),] 
head(top.dec)
```

```
##             mean    total      bio      tech p.value FDR
## LYZ     1.907466 5.057544 4.084819 0.9727253       0   0
## S100A9  1.877589 4.571260 3.604340 0.9669195       0   0
## S100A8  1.682848 4.504609 3.574976 0.9296337       0   0
## HLA-DRA 1.998496 3.603509 2.613700 0.9898088       0   0
## CD74    2.726807 3.403452 2.307571 1.0958811       0   0
## CST3    1.405345 2.815407 1.945260 0.8701464       0   0
```

We can plot the genes with the largest biological components, to verify that they are indeed highly variable.


```r
plotExpression(sce, feature=rownames(top.dec)[1:10])
```

<img src="analysis_files/figure-html/hvgplot-1.png" width="100%"  class="widefigure" />

# Dimensionality reduction

We use the `denoisePCA` function to perform PCA, using the assumed Poisson technical trend.


```r
sce <- denoisePCA(sce, technical=new.trend, approx=TRUE)
ncol(reducedDim(sce, "PCA"))
```

```
## [1] 100
```

```r
plot(attr(reducedDim(sce), "percentVar"))
```

<img src="analysis_files/figure-html/unnamed-chunk-15-1.png" width="100%" />

We can plot the first few components.


```r
plotPCA(sce, ncomponents=3, colour_by="Detection")
```

<img src="analysis_files/figure-html/pcaplot-1.png" width="100%" />

Same with using _t_-SNE for visualization.


```r
sce <- runTSNE(sce, use_dimred="PCA", perplexity=20, rand_seed=100)
plotTSNE(sce, colour_by="Detection")
```

<img src="analysis_files/figure-html/unnamed-chunk-16-1.png" width="100%" />

# Clustering with graph-based methods

We use the shared nearest neighbour method for clustering.


```r
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
table(clusters$membership)
```

```
## 
##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
##   60 1015  400 1314  238  219  513   56   67   54  452  174   92   23   51 
##   16 
##   36
```

Plotting them out to verify separateness.


```r
sce$Cluster <- factor(clusters$membership)
plotTSNE(sce, colour_by="Cluster")
```

<img src="analysis_files/figure-html/unnamed-chunk-18-1.png" width="100%" />

Also examining their modularity scores.
We look at the ratio of the observed and expected edge weights, as the raw modularity varies by orders of magnitudes across clusters.


```r
cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.values=TRUE)
log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)
library(pheatmap)
pheatmap(log.ratio, cluster_rows=FALSE, cluster_cols=FALSE, 
    color=colorRampPalette(c("white", "blue"))(100))
```

<img src="analysis_files/figure-html/unnamed-chunk-19-1.png" width="100%" />

# Marker gene detection

Detecting marker genes for each cluster.


```r
marker.out <- findMarkers(sce, clusters=sce$Cluster)
```

Having a look at how the clusters interact with the detection status, so we can focus on EmptyDrops-unique clusters.


```r
table(sce$Cluster, sce$Detection)
```

```
##     
##      Both CellRanger EmptyDrops
##   1    57          0          3
##   2   903         12        100
##   3   356         33         11
##   4  1106        203          5
##   5     6          0        232
##   6   197         22          0
##   7   499          6          8
##   8    55          1          0
##   9    67          0          0
##   10    3          0         51
##   11  418         31          3
##   12  161          0         13
##   13   91          0          1
##   14   20          1          2
##   15   49          2          0
##   16   36          0          0
```

Focusing on cluster 10, which seems to be made of platelets:


```r
current <- marker.out[["10"]]
chosen <- rownames(current)[current$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts", 
    zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
    colour_columns_by="Cluster", columns=order(sce$Cluster))
```

<img src="analysis_files/figure-html/heatmap8-1.png" width="100%"  class="widefigure" />

... and 5, which seems to contain damaged cells with high mitochondrial content:


```r
current <- marker.out[["5"]]
chosen <- rownames(current)[current$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts", 
    zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
    colour_columns_by="Cluster", columns=order(sce$Cluster))
```

<img src="analysis_files/figure-html/heatmap6-1.png" width="100%"  class="widefigure" />

# Wrapping up

We save the various bits and pieces for further plotting.


```r
saveRDS(sce, file="sce.rds")
```

Printing the session information.


```r
sessionInfo()
```

```
## R Under development (unstable) (2017-10-27 r73632)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.5 LTS
## 
## Matrix products: default
## BLAS: /home/cri.camres.org/lun01/Software/R/trunk/lib/libRblas.so
## LAPACK: /home/cri.camres.org/lun01/Software/R/trunk/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] pheatmap_1.0.8             scran_1.7.13              
##  [3] scater_1.7.4               ggplot2_2.2.1             
##  [5] EnsDb.Hsapiens.v86_2.99.0  ensembldb_2.3.7           
##  [7] AnnotationFilter_1.3.0     GenomicFeatures_1.31.1    
##  [9] AnnotationDbi_1.41.4       DropletUtils_0.99.14      
## [11] SingleCellExperiment_1.1.2 SummarizedExperiment_1.9.7
## [13] DelayedArray_0.5.16        matrixStats_0.52.2        
## [15] Biobase_2.39.1             GenomicRanges_1.31.6      
## [17] GenomeInfoDb_1.15.1        IRanges_2.13.10           
## [19] S4Vectors_0.17.22          BiocGenerics_0.25.1       
## [21] BiocParallel_1.13.1        knitr_1.18                
## [23] BiocStyle_2.7.5           
## 
## loaded via a namespace (and not attached):
##  [1] ProtGenerics_1.11.0      bitops_1.0-6            
##  [3] bit64_0.9-7              RColorBrewer_1.1-2      
##  [5] progress_1.1.2           httr_1.3.1              
##  [7] rprojroot_1.3-2          dynamicTreeCut_1.63-1   
##  [9] tools_3.5.0              backports_1.1.2         
## [11] irlba_2.3.2              DT_0.2                  
## [13] R6_2.2.2                 vipor_0.4.5             
## [15] DBI_0.7                  lazyeval_0.2.1          
## [17] colorspace_1.3-2         gridExtra_2.3           
## [19] prettyunits_1.0.2        RMySQL_0.10.13          
## [21] bit_1.1-12               curl_3.1                
## [23] compiler_3.5.0           rtracklayer_1.39.5      
## [25] labeling_0.3             bookdown_0.5            
## [27] scales_0.5.0             stringr_1.2.0           
## [29] digest_0.6.14            Rsamtools_1.31.1        
## [31] rmarkdown_1.8            XVector_0.19.5          
## [33] pkgconfig_2.0.1          htmltools_0.3.6         
## [35] limma_3.35.5             htmlwidgets_0.9         
## [37] rlang_0.1.6              RSQLite_2.0             
## [39] FNN_1.1                  shiny_1.0.5             
## [41] DelayedMatrixStats_1.1.4 bindr_0.1               
## [43] dplyr_0.7.4              RCurl_1.95-4.10         
## [45] magrittr_1.5             GenomeInfoDbData_1.1.0  
## [47] Matrix_1.2-12            Rcpp_0.12.14            
## [49] ggbeeswarm_0.6.0         munsell_0.4.3           
## [51] Rhdf5lib_1.1.5           viridis_0.4.1           
## [53] stringi_1.1.6            yaml_2.1.16             
## [55] edgeR_3.21.6             zlibbioc_1.25.0         
## [57] Rtsne_0.13               rhdf5_2.23.4            
## [59] plyr_1.8.4               grid_3.5.0              
## [61] blob_1.1.0               shinydashboard_0.6.1    
## [63] lattice_0.20-35          cowplot_0.9.2           
## [65] Biostrings_2.47.2        locfit_1.5-9.1          
## [67] pillar_1.1.0             igraph_1.1.2            
## [69] rjson_0.2.15             reshape2_1.4.3          
## [71] biomaRt_2.35.8           XML_3.98-1.9            
## [73] glue_1.2.0               evaluate_0.10.1         
## [75] data.table_1.10.4-3      httpuv_1.3.5            
## [77] gtable_0.2.0             assertthat_0.2.0        
## [79] mime_0.5                 xtable_1.8-2            
## [81] viridisLite_0.2.0        tibble_1.4.1            
## [83] GenomicAlignments_1.15.4 beeswarm_0.2.3          
## [85] memoise_1.1.0            tximport_1.7.4          
## [87] bindrcpp_0.2             statmod_1.4.30
```
