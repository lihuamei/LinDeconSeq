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
GSM495224 2.952270e-01 0.6509099846 5.386302e-02
GSM495225 2.943407e-01 0.6491998940 5.645944e-02
GSM495226 3.059917e-01 0.6407750804 5.323327e-02
GSM495227 6.280094e-01 0.3051263064 6.686430e-02
GSM495228 6.344095e-01 0.3010209770 6.456955e-02
GSM495229 6.319751e-01 0.3010156579 6.700925e-02
GSM495230 4.464379e-01 0.4494756881 1.040864e-01
GSM495231 4.688540e-01 0.4321675433 9.897847e-02
GSM495232 4.735214e-01 0.4321313647 9.434719e-02
GSM495233 5.304824e-01 0.2391622392 2.303553e-01
GSM495234 5.329330e-01 0.2389048131 2.281621e-01
GSM495235 5.331426e-01 0.2363373910 2.305200e-01
GSM495236 4.927669e-01 0.3235672910 1.836658e-01
GSM495237 4.988868e-01 0.3172577046 1.838555e-01
GSM495238 4.907264e-01 0.3250466413 1.842269e-01
GSM495239 5.356287e-01 0.3215945800 1.427767e-01
GSM495240 5.317241e-01 0.3218192892 1.464566e-01
GSM495241 5.305241e-01 0.3285337741 1.409422e-01
GSM495242 4.945030e-01 0.4050589522 1.004381e-01
GSM495243 5.092395e-01 0.3925026497 9.825782e-02
GSM495244 5.008963e-01 0.3986931994 1.004105e-01
GSM495245 5.743211e-01 0.3657359346 5.994296e-02
GSM495246 5.720378e-01 0.3648056213 6.315654e-02
GSM495247 5.770924e-01 0.3629614611 5.994618e-02
GSM495248 6.057886e-01 0.3694192456 2.479219e-02
GSM495249 6.110662e-01 0.3662804061 2.265338e-02
GSM495250 6.090999e-01 0.3666393370 2.426080e-02


```
![LinDeconSeq\_fractions](data/fractions.png)

# Citation
----------------------------
Please cite the publication: ***Li H, Sharma A, Ming W, et al. A deconvolution method and its application in analyzing the cellular fractions in acute myeloid leukemia samples[J]. BMC genomics, 2020, 21(1): 1-15.***<br>
