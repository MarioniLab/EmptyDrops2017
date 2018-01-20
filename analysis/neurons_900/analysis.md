---
title: Analyzing the 4K PBMC dataset
author: Aaron Lun and others
date: 28 December 2017
output: 
  BiocStyle::html_document:
    fig_caption: no
---



# Introduction

This document performs a brief analysis of the 900 neuron dataset, with particular focus on the cells that are uniquely detected by EmptyDrops or CellRanger.
First, we load in the raw count matrix:


```r
library(DropletUtils)
fname <- "../../data/neurons_900/raw_gene_bc_matrices/mm10"
sce <- read10xCounts(fname, col.names=TRUE)
sce
```

```
## class: SingleCellExperiment 
## dim: 27998 737280 
## metadata(0):
## assays(1): counts
## rownames(27998): ENSMUSG00000051951 ENSMUSG00000089699 ...
##   ENSMUSG00000096730 ENSMUSG00000095742
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
## logical    1580    1974  733726
```

... and CellRanger (with an expected 4000 cells, as the name of the dataset suggests):


```r
c.keep <- defaultDrops(counts(sce), expected=900)
summary(c.keep)
```

```
##    Mode   FALSE    TRUE 
## logical  736488     792
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
##       Both EmptyDrops 
##     736098       1182
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
library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))
```

```
## DataFrame with 6 rows and 4 columns
##                   ID      Symbol            ENSEMBL      SYMBOL
##          <character> <character>        <character> <character>
## 1 ENSMUSG00000051951        Xkr4 ENSMUSG00000051951        Xkr4
## 2 ENSMUSG00000089699      Gm1992 ENSMUSG00000089699          NA
## 3 ENSMUSG00000102343     Gm37381 ENSMUSG00000102343          NA
## 4 ENSMUSG00000025900         Rp1 ENSMUSG00000025900         Rp1
## 5 ENSMUSG00000109048         Rp1 ENSMUSG00000109048          NA
## 6 ENSMUSG00000025902       Sox17 ENSMUSG00000025902       Sox17
```

We relabel the rows with the gene symbols for easier reading.


```r
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)
head(rownames(sce))
```

```
## [1] "Xkr4"               "ENSMUSG00000089699" "ENSMUSG00000102343"
## [4] "Rp1"                "ENSMUSG00000109048" "Sox17"
```

We also determine the chromosomal location for each gene.


```r
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
    column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")
```

```
##    Mode   FALSE    TRUE    NA's 
## logical   21767      13    6218
```

# Quality control on the cells

Cell detection can be considered an implicit quality control step, so technically, no further work is needed.
Nonetheless, we examine some commonly used metrics.


```r
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=="chrM")))
par(mfrow=c(1,3))
hist(sce$log10_total_counts, breaks=20, col="grey80")
hist(sce$log10_total_features_by_counts, breaks=20, col="grey80")
hist(sce$pct_counts_Mito, breaks=20, col="grey80")
```

<img src="analysis_files/figure-html/qchist-1.png" width="100%"  class="widefigure" />

A large number of the features with low total counts have high mitochondrial proportions.
This is unsurprising and probably corresponds to damaged cells - we'll filter them out.


```r
par(mfrow=c(1,2))
plot(sce$total_features_by_counts, sce$pct_counts_Mito)
plot(sce$total_counts, sce$pct_counts_Mito)
```

<img src="analysis_files/figure-html/mitoscatter-1.png" width="100%"  class="widefigure" />

```r
discard <- isOutlier(sce$pct_counts_Mito)
sce <- sce[,!discard]
summary(discard)
```

```
##    Mode   FALSE    TRUE 
## logical    1862     112
```

We remove a few cells that have high haemoglobin expression; these are probably red blood cells that slipped through.
We'll remove them because (i) we're looking at brain cells and (ii) RBCs have so few features that they cause problems later.


```r
library(Matrix)
haem.genes <- grep("^Hbb", rownames(sce))
haem.prop <- colSums(counts(sce)[haem.genes,])/colSums(counts(sce))
plot(sce$total_counts, haem.prop)
```

<img src="analysis_files/figure-html/haemscatter-1.png" width="100%" />

```r
discard <- haem.prop > 0.1
sce <- sce[,!discard]
summary(discard)
```

```
##    Mode   FALSE    TRUE 
## logical    1790      72
```

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

We cluster cells together to improve the accuracy of the pooling later.


```r
library(scran)
clusters <- quickCluster(sce, method="igraph", subset.row=ave >= 0.1)
table(clusters)
```

```
## clusters
##   1   2   3 
## 313 594 883
```

We then use the deconvolution method to compute size factors for each cell.


```r
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
summary(sizeFactors(sce))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  0.004727  0.111507  0.607179  1.000000  1.316854 13.350225
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
##            mean     total      bio     tech       p.value           FDR
## Malat1 5.055957 10.098168 7.103930 2.994237  0.000000e+00  0.000000e+00
## Fabp7  2.063064  6.206551 4.153429 2.053122  0.000000e+00  0.000000e+00
## Dbi    2.030756  5.004354 2.974121 2.030233 2.189246e-221 1.474813e-219
## Vim    1.258004  3.297017 1.889667 1.407350 1.252528e-193 7.119403e-192
## Hmgb2  1.200557  3.204983 1.851686 1.353298 2.370964e-199 1.411370e-197
## Tubb3  2.535708  4.155165 1.835223 2.319942  2.106622e-83  5.201419e-82
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
## [1] 5
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
sce <- runTSNE(sce, use_dimred="PCA", perplexity=10, rand_seed=100)
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
##   1   2   3   4   5   6   7   8   9 
## 287  90 188 314 183 182 336  84 126
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
##     Both EmptyDrops
##   1    0        287
##   2    0         90
##   3    4        184
##   4  146        168
##   5  129         54
##   6  154         28
##   7  329          7
##   8    0         84
##   9    9        117
```

Focusing on cluster 1:


```r
current <- marker.out[["1"]]
chosen <- rownames(current)[current$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts", 
    zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
    colour_columns_by="Cluster", columns=order(sce$Cluster))
```

<img src="analysis_files/figure-html/heatmap8-1.png" width="100%"  class="widefigure" />

... and 9, which seems to be interneurons:


```r
current <- marker.out[["9"]]
chosen <- rownames(current)[current$Top <= 20]
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
##  [1] pheatmap_1.0.8                        
##  [2] scran_1.7.13                          
##  [3] Matrix_1.2-12                         
##  [4] TxDb.Mmusculus.UCSC.mm10.ensGene_3.4.0
##  [5] GenomicFeatures_1.31.1                
##  [6] scater_1.7.4                          
##  [7] ggplot2_2.2.1                         
##  [8] org.Mm.eg.db_3.5.0                    
##  [9] AnnotationDbi_1.41.4                  
## [10] DropletUtils_0.99.14                  
## [11] SingleCellExperiment_1.1.2            
## [12] SummarizedExperiment_1.9.7            
## [13] DelayedArray_0.5.16                   
## [14] matrixStats_0.52.2                    
## [15] Biobase_2.39.1                        
## [16] GenomicRanges_1.31.6                  
## [17] GenomeInfoDb_1.15.1                   
## [18] IRanges_2.13.10                       
## [19] S4Vectors_0.17.22                     
## [20] BiocGenerics_0.25.1                   
## [21] BiocParallel_1.13.1                   
## [22] knitr_1.18                            
## [23] BiocStyle_2.7.5                       
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7             
##  [3] RColorBrewer_1.1-2       progress_1.1.2          
##  [5] httr_1.3.1               rprojroot_1.3-2         
##  [7] dynamicTreeCut_1.63-1    tools_3.5.0             
##  [9] backports_1.1.2          irlba_2.3.2             
## [11] DT_0.2                   R6_2.2.2                
## [13] vipor_0.4.5              DBI_0.7                 
## [15] lazyeval_0.2.1           colorspace_1.3-2        
## [17] gridExtra_2.3            prettyunits_1.0.2       
## [19] RMySQL_0.10.13           bit_1.1-12              
## [21] compiler_3.5.0           labeling_0.3            
## [23] rtracklayer_1.39.5       bookdown_0.5            
## [25] scales_0.5.0             stringr_1.2.0           
## [27] digest_0.6.14            Rsamtools_1.31.1        
## [29] rmarkdown_1.8            XVector_0.19.5          
## [31] pkgconfig_2.0.1          htmltools_0.3.6         
## [33] limma_3.35.5             htmlwidgets_0.9         
## [35] rlang_0.1.6              RSQLite_2.0             
## [37] FNN_1.1                  shiny_1.0.5             
## [39] DelayedMatrixStats_1.1.4 bindr_0.1               
## [41] dplyr_0.7.4              RCurl_1.95-4.10         
## [43] magrittr_1.5             GenomeInfoDbData_1.1.0  
## [45] Rcpp_0.12.14             ggbeeswarm_0.6.0        
## [47] munsell_0.4.3            Rhdf5lib_1.1.5          
## [49] viridis_0.4.1            stringi_1.1.6           
## [51] yaml_2.1.16              edgeR_3.21.6            
## [53] zlibbioc_1.25.0          Rtsne_0.13              
## [55] rhdf5_2.23.4             plyr_1.8.4              
## [57] grid_3.5.0               blob_1.1.0              
## [59] shinydashboard_0.6.1     lattice_0.20-35         
## [61] cowplot_0.9.2            Biostrings_2.47.2       
## [63] locfit_1.5-9.1           pillar_1.1.0            
## [65] igraph_1.1.2             rjson_0.2.15            
## [67] reshape2_1.4.3           biomaRt_2.35.8          
## [69] XML_3.98-1.9             glue_1.2.0              
## [71] evaluate_0.10.1          data.table_1.10.4-3     
## [73] httpuv_1.3.5             gtable_0.2.0            
## [75] assertthat_0.2.0         mime_0.5                
## [77] xtable_1.8-2             viridisLite_0.2.0       
## [79] tibble_1.4.1             GenomicAlignments_1.15.4
## [81] beeswarm_0.2.3           memoise_1.1.0           
## [83] tximport_1.7.4           bindrcpp_0.2            
## [85] statmod_1.4.30
```
