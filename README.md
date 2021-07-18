LinDeconSeq
===================================================

Cell Type Deconvolution using gene expression data for bulk samples.
![LinDeconSeq\_pipeline](data/pipeline.jpg)

# Installing the package
---------------------
You can install the package using devtools::install_github:

``` r
devtools::install_github("lihuamei/LinDeconSeq")
```

# Getting started with LinDeconSeq
----------------------------
In this tutorial we will use GSE19830 (mixture of Liver, Brain and Lung) as an example.

``` r
library(LinDeconSeq)
pures <- shen_orr$data[, rowSums(shen_orr$annotation$pure) != 0 ] %>% t
markerRes <- findMarkers(pures[, colnames(shen_orr$phenotype)], shen_orr$phenotype, min.group = 100, max.group = 300, norm.method = 'QN', data.type = 'MA')


```

# Deconvolution
----------------------------
To deconvolve the dataset, signature marker genes must be known in advance.

```r

fractions <- deconSeq(shen_orr$data %>% t, markerRes$sigMatrix$sig.mat, verbose = TRUE)

```
![LinDeconSeq\_fractions](data/fractions.png)

# Citation
----------------------------
Please cite the publication: ***Li H, Sharma A, Ming W, et al. A deconvolution method and its application in analyzing the cellular fractions in acute myeloid leukemia samples[J]. BMC genomics, 2020, 21(1): 1-15.***<br>
