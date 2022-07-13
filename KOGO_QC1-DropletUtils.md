KOGO_QC1(DropletUtils)
================
2022-07-04

### **Library**

``` r
library(devtools)
library(Seurat)
library(scRNAseq)
library(scater)
library(harmony)
library(DropletUtils)
library(scran)
```

### **Set Directory Path**
``` r
rdatadir = './kogo2022/QC/'
```

### **Load Data**

``` r
rawsce_1 <- read10xCounts(paste0(rdatadir, "20094_0001_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_2 <- read10xCounts(paste0(rdatadir, "20094_0002_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_3 <- read10xCounts(paste0(rdatadir, "20094_0003_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_4 <- read10xCounts(paste0(rdatadir, "20094_0004_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_5 <- read10xCounts(paste0(rdatadir, "20094_0005_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_6 <- read10xCounts(paste0(rdatadir, "20094_0006_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_7 <- read10xCounts(paste0(rdatadir, "20094_0007_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_8 <- read10xCounts(paste0(rdatadir, "20094_0008_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_9 <- read10xCounts(paste0(rdatadir, "20094_0009_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_12 <- read10xCounts(paste0(rdatadir, "20094_0012_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
```

### **Check Data**

``` r
rawsce_1
```

    ## class: SingleCellExperiment 
    ## dim: 36604 6794880 
    ## metadata(1): Samples
    ## assays(1): counts
    ## rownames(36604): ENSG00000243485 ENSG00000237613 ... Htag2 Htag3
    ## rowData names(3): ID Symbol Type
    ## colnames: NULL
    ## colData names(2): Sample Barcode
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):

### **Run DropletUtils**

``` r
rawsce_list <- c('rawsce_1', 'rawsce_2', 'rawsce_3', 'rawsce_4', 'rawsce_5',
                 'rawsce_6', 'rawsce_7', 'rawsce_8', 'rawsce_9', 'rawsce_12')

for (i in 1:length(rawsce_list)) {
  rawsce <- get(rawsce_list[i])
  
  br.out <- barcodeRanks(counts(rawsce))

  png(paste0('DropletUtils_', rawsce_list[i], '.png'))
  
  plot(br.out$rank, br.out$total, log = "xy", xlab = "Rank", ylab = "Total", main = rawsce_list[i])
  
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col = "red")
  
  set.seed(2022)
  e.out <- emptyDrops(counts(rawsce))  ## Cells that have UMI counts lower than 100 are empty cells.
  table(Sig=e.out$FDR <= 0.05, Limited=e.out$Limited)
  is.cell <- e.out$FDR <= 0.05
  
  print(sum(is.cell, na.rm=TRUE))
  print(table(br.out$rank == sum(is.cell, na.rm=TRUE)))
  
  abline(h=min(br.out$fitted[o], na.rm=TRUE), col="red", lty=2)
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), legend=c("knee", "inflection", "FDR_0.05"))
  dev.off()
  
  colnames(rawsce) = colData(rawsce)$Barcode
  rawsce <- rawsce[,which(e.out$FDR <= 0.05)]
  
  assign(paste0('DropletUtils_', rawsce_list[i]), rawsce)
}
```

    ## [1] 12732
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 7273
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 6490
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 4485
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 4091
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 2671
    ## 
    ##   FALSE    TRUE 
    ## 6794877       3

    ## [1] 2001
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 5136
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 2839
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 7164
    ## 
    ##   FALSE 
    ## 6794880

### **Check DropletUtils data**

``` r
DropletUtils_rawsce_1
```

    ## class: SingleCellExperiment 
    ## dim: 36604 12732 
    ## metadata(1): Samples
    ## assays(1): counts
    ## rownames(36604): ENSG00000243485 ENSG00000237613 ... Htag2 Htag3
    ## rowData names(3): ID Symbol Type
    ## colnames(12732): AAACCCAAGGCTGAAC-1 AAACCCAAGGGTGAAA-1 ...
    ##   TTTGTTGTCTTCGACC-1 TTTGTTGTCTTGGTCC-1
    ## colData names(2): Sample Barcode
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):

### **Reference**
Lun, A. T., Riesenfeld, S., Andrews, T., Gomes, T., & Marioni, J. C. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome biology, 20(1), 1-9.

Pekayvaz, K., Leunig, A., Kaiser, R., Joppich, M., Brambs, S., Janjic, A., ... & Nicolai, L. (2022). Protective immune trajectories in early viral containment of non-pneumonic SARS-CoV-2 infection. Nature communications, 13(1), 1-21.
