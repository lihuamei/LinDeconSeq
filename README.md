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
markerRes <- findMarkers(pures[, colnames(shen_orr$phenotype)], shen_orr$phenotype, QN = TRUE)
[INFO] 9 samples and 6669 genes in the reference profile
[INFO] 2873 candidate cell type-specific genes detected with q.cut 0.01
[INFO] Select seed genes and random permutating...
[INFO] Assigning cell type-specific genes to each cell subset
[INFO] Optimizing cell type-specific genes to derive signature matrix...
[INFO] Group size -> 50, condition number -> 1.518792
[INFO] 979 low confidence marker genes have been filtered out.

```

# Deconvolution
----------------------------
To deconvolve the dataset, signature marker genes must be known in advance.

```r

fractions <- deconSeq(shen_orr$data %>% t, markerRes$sigMatrix$sig.mat, verbose = TRUE)
[INFO] There are 42 bulk samples need to be deconvoluted
[INFO] Show deconvolution results....
                 Liver        Brain         Lung
GSM495209 9.995712e-01 0.0004288179 0.000000e+00
GSM495210 9.995061e-01 0.0004938515 0.000000e+00
GSM495211 9.956582e-01 0.0036376592 7.041480e-04
GSM495212 0.000000e+00 1.0000000000 0.000000e+00
GSM495213 9.291706e-04 0.9986788620 3.919674e-04
GSM495214 5.796167e-04 0.9993892521 3.113119e-05
GSM495215 4.002953e-06 0.0000000000 9.999960e-01
GSM495216 0.000000e+00 0.0000000000 1.000000e+00
GSM495217 0.000000e+00 0.0000000000 1.000000e+00
GSM495218 8.380125e-02 0.3066126236 6.095861e-01
GSM495219 9.013318e-02 0.3073880816 6.024787e-01
GSM495220 8.759564e-02 0.3053518980 6.070525e-01
GSM495221 6.544839e-01 0.0914730557 2.540430e-01
GSM495222 6.490951e-01 0.0914795481 2.594253e-01
GSM495223 6.523207e-01 0.0887322675 2.589470e-01
...

```
![LinDeconSeq\_fractions](data/fractions.png)

# More Information
--------------------
Please see [Tutorial](https://lihuamei.github.io//LinDeconSeq/inst/LinDeconSeq.html).

# Citation
----------------------------
Please cite the publication: ***Li H, Sharma A, Ming W, et al. A deconvolution method and its application in analyzing the cellular fractions in acute myeloid leukemia samples[J]. BMC genomics, 2020, 21(1): 1-15.***<br>
