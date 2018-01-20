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

Cell detection can be considered an implicit quality control step, so technically, no extra steps are needed.
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

We perform some pre-clustering to break up obvious clusters.


```r
library(scran)
clusters <- quickCluster(sce, method="igraph", subset.row=ave>=0.1, 
    irlba.args=list(maxit=1000)) # for convergence.
table(clusters)
```

```
## clusters
##    1    2    3 
##  415 2347 2002
```

We then use the deconvolution method to compute size factors for each cell.


```r
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
summary(sizeFactors(sce))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -0.004299  0.660967  0.913900  0.999158  1.162356 10.936650
```

We can plot these against the library sizes to see how much of a difference it makes.


```r
plot(sce$total_counts, sizeFactors(sce), log="xy")
```

<img src="analysis_files/figure-html/sfplot-1.png" width="100%" />

Note that some size factors are very small and negative. 
This represents cells that have so few expressed features that it is not possible to obtain a sensible size factor.


```r
neg.sf <- sizeFactors(sce)<0
summary(neg.sf)
```

```
##    Mode   FALSE    TRUE 
## logical    4760       4
```

Instead, we replace the size factor with the (scaled) library size.


```r
library(Matrix)
lib.sizes <- colSums(counts(sce))
scaled.lib.sizes <- lib.sizes/mean(lib.sizes)
sizeFactors(sce)[neg.sf] <- scaled.lib.sizes[neg.sf]
```

Finally, we compute normalized log-expresion values.


```r
sce <- normalize(sce)
```

# Modelling the mean-variance trend

We assume that the technical noise is Poisson and create a fitted trend on that basis.


```r
new.trend <- makeTechTrend(x=sce)
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
## LYZ     1.913139 5.068595 4.056640 1.0119545       0   0
## S100A9  1.886379 4.602469 3.595864 1.0066042       0   0
## S100A8  1.693172 4.561085 3.592223 0.9688621       0   0
## HLA-DRA 2.010263 3.686198 2.655423 1.0307744       0   0
## CD74    2.739638 3.496692 2.359398 1.1372939       0   0
## CST3    1.407413 2.824938 1.922024 0.9029148       0   0
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
## [1] 20
```

```r
plot(attr(reducedDim(sce), "percentVar"))
```

<img src="analysis_files/figure-html/unnamed-chunk-17-1.png" width="100%" />

We can plot the first few components.


```r
plotPCA(sce, ncomponents=3, colour_by="Detection")
```

<img src="analysis_files/figure-html/pcaplot-1.png" width="100%" />

Same with using _t_-SNE for visualization.


```r
sce <- runTSNE(sce, use_dimred="PCA", perplexity=30, rand_seed=100)
plotTSNE(sce, colour_by="Detection")
```

<img src="analysis_files/figure-html/unnamed-chunk-18-1.png" width="100%" />

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
##  233  968  538  377  619  125 1151  174  102  211   93   35   32   55   15 
##   16 
##   36
```

Plotting them out to verify separateness.


```r
sce$Cluster <- factor(clusters$membership)
plotTSNE(sce, colour_by="Cluster")
```

<img src="analysis_files/figure-html/unnamed-chunk-20-1.png" width="100%" />

Also examining their modularity scores.
We look at the ratio of the observed and expected edge weights, as the raw modularity varies by orders of magnitudes across clusters.


```r
cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.values=TRUE)
log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)
library(pheatmap)
pheatmap(log.ratio, cluster_rows=FALSE, cluster_cols=FALSE, 
    color=colorRampPalette(c("white", "blue"))(100))
```

<img src="analysis_files/figure-html/unnamed-chunk-21-1.png" width="100%" />

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
##   1    13          0        220
##   2   939         14         15
##   3   519         10          9
##   4   356         18          3
##   5   554         55         10
##   6   124          1          0
##   7   964        181          6
##   8   161          0         13
##   9    18          0         84
##   10  179         31          1
##   11   93          0          0
##   12   35          0          0
##   13   29          1          2
##   14    4          0         51
##   15    0          0         15
##   16   36          0          0
```

Focusing on cluster 14, which seems to be made of platelets:


```r
current <- marker.out[["14"]]
chosen <- rownames(current)[current$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts", 
    zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
    colour_columns_by="Cluster", columns=order(sce$Cluster))
```

<img src="analysis_files/figure-html/heatmap8-1.png" width="100%"  class="widefigure" />

... and 5, which seems to contain damaged cells with high mitochondrial content:


```r
current <- marker.out[["1"]]
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
##  [1] pheatmap_1.0.8             Matrix_1.2-12             
##  [3] scran_1.7.13               scater_1.7.4              
##  [5] ggplot2_2.2.1              EnsDb.Hsapiens.v86_2.99.0 
##  [7] ensembldb_2.3.7            AnnotationFilter_1.3.0    
##  [9] GenomicFeatures_1.31.1     AnnotationDbi_1.41.4      
## [11] DropletUtils_0.99.14       SingleCellExperiment_1.1.2
## [13] SummarizedExperiment_1.9.7 DelayedArray_0.5.16       
## [15] matrixStats_0.52.2         Biobase_2.39.1            
## [17] GenomicRanges_1.31.6       GenomeInfoDb_1.15.1       
## [19] IRanges_2.13.10            S4Vectors_0.17.22         
## [21] BiocGenerics_0.25.1        BiocParallel_1.13.1       
## [23] knitr_1.18                 BiocStyle_2.7.5           
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
## [47] Rcpp_0.12.14             ggbeeswarm_0.6.0        
## [49] munsell_0.4.3            Rhdf5lib_1.1.5          
## [51] viridis_0.4.1            stringi_1.1.6           
## [53] yaml_2.1.16              edgeR_3.21.6            
## [55] zlibbioc_1.25.0          Rtsne_0.13              
## [57] rhdf5_2.23.4             plyr_1.8.4              
## [59] grid_3.5.0               blob_1.1.0              
## [61] shinydashboard_0.6.1     lattice_0.20-35         
## [63] cowplot_0.9.2            Biostrings_2.47.2       
## [65] locfit_1.5-9.1           pillar_1.1.0            
## [67] igraph_1.1.2             rjson_0.2.15            
## [69] reshape2_1.4.3           biomaRt_2.35.8          
## [71] XML_3.98-1.9             glue_1.2.0              
## [73] evaluate_0.10.1          data.table_1.10.4-3     
## [75] httpuv_1.3.5             gtable_0.2.0            
## [77] assertthat_0.2.0         mime_0.5                
## [79] xtable_1.8-2             viridisLite_0.2.0       
## [81] tibble_1.4.1             GenomicAlignments_1.15.4
## [83] beeswarm_0.2.3           memoise_1.1.0           
## [85] tximport_1.7.4           bindrcpp_0.2            
## [87] statmod_1.4.30
```
