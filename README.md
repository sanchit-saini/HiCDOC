HiCDOC normalizes intrachromosomal Hi-C matrices and detects A/B compartments
with multiple replicates using unsupervised learning.

## Prerequisites

- R >= 3.5

## Installation

- Open `R`
- Install `devtools` if not already installed

    ```R
    install.packages("devtools")
    ```
        
- Install the package:

    ```R
    library(devtools)
    install_github("mzytnicki/HiCDOC")
     ```

- Load the package

    ```R
    library(HiCDOC)
    ```


HiCDOC is segmented in a collection of functions, which can be used to construct a
custom pipeline.

1. Normalize technical biases with using a cyclic loess<sup>[[publication][multihiccompare-publication]][[implementation][cyclic-loess-implementation]]</sup>
2. Normalize biological biases with Knight-Ruiz method<sup>[[publication][knight-ruiz-publication]][[implementation][knight-ruiz-implementation]]</sup>
3. Normalize distance effect with an MD loess
4. Detect compartments with constrained k-means<sup>[[publication][constrained-k-means-publication]][[implementation][constrained-k-means-implementation]]</sup>
5. Plot compartment changes and measures 

## Input format

Each script accepts a tab-separated multi-replicate sparse matrix with a header
line and optional comment lines.

    # Optional comment lines
    # ...
    chromosome    position 1    position 2    replicate 1.1    replicate 1.2    replicate 2.1    ...
    3             1500000       7500000       145              184              72               ...
    ...

The interaction proportions between `position 1` and `position 2` are reported
in each replicate column, named `replicate <condition.replicate>`. There is no
limit to the number of replicates and conditions.

## Start with HiCDOC

A small dataset is shipped with the package:

```R
object <- HiCDOCExample()
```

The usual pipeline is then:

```R
object <- parseInteractionMatrix3Columns(object)
plotInteractionMatrix(object, log = TRUE)
object <- normalizeCyclicLoess(object)
plotInteractionMatrix(object, log = TRUE)
object <- normalizeKnightRuiz(object)
plotInteractionMatrix(object, log = TRUE)
plotMD(object)
object <- normalizeDistanceCombined(object)
plotMD(object)
plotInteractionMatrix(object, log = FALSE)
object <- detectConstrainedKMeans(object)
object <- findPValues(object)
plotConcordances(object)
plotCompartmentChanges(object)
DIR(object, pvalue = 1)
concordances(object)
compartments(object)
```


## References

Philip A. Knight, Daniel Ruiz, A fast algorithm for matrix balancing, _IMA
Journal of Numerical Analysis_, Volume 33, Issue 3, July 2013, Pages 1029–1047,
https://doi.org/10.1093/imanum/drs019

Rajendra Kumar, Haitham Sobhy, Per Stenberg, Ludvig Lizana, Genome contact map
explorer: a platform for the comparison, interactive visualization and analysis
of genome contact maps, _Nucleic Acids Research_, Volume 45, Issue 17, 29
September 2017, Page e152, https://doi.org/10.1093/nar/gkx644

John C Stansfield, Kellen G Cresswell, Mikhail G Dozmorov, multiHiCcompare:
joint normalization and comparative analysis of complex Hi-C experiments,
_Bioinformatics_, 2019, https://doi.org/10.1093/bioinformatics/btz048

Kiri Wagstaff, Claire Cardie, Seth Rogers, Stefan Schrödl, Constrained K-means
Clustering with Background Knowledge, _Proceedings of 18th International
Conference on Machine Learning_, 2001, Pages 577-584,
https://pdfs.semanticscholar.org/0bac/ca0993a3f51649a6bb8dbb093fc8d8481ad4.pdf

[multihiccompare-publication]: https://doi.org/10.1093/bioinformatics/btz048
[multihiccompare-installation]: https://bioconductor.org/packages/release/bioc/html/multiHiCcompare.html
[gcmapexplorer-publication]: https://doi.org/10.1093/nar/gkx644
[gcmapexplorer-installation]: https://gcmapexplorer.readthedocs.io/en/latest/install.html
[orca-installation]: https://github.com/plotly/orca#installation
[cyclic-loess-implementation]: https://bioconductor.org/packages/release/bioc/vignettes/multiHiCcompare/inst/doc/multiHiCcompare.html#cyclic-loess-normalization
[knight-ruiz-publication]: https://doi.org/10.1093/imanum/drs019
[knight-ruiz-implementation]: https://gcmapexplorer.readthedocs.io/en/latest/commands/normKR.html
[rnr-implementation]: https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.RadiusNeighborsRegressor.html
[interaction-mean-implementation]: https://gcmapexplorer.readthedocs.io/en/latest/commands/normMCFS.html
[constrained-k-means-publication]: https://pdfs.semanticscholar.org/0bac/ca0993a3f51649a6bb8dbb093fc8d8481ad4.pdf
[constrained-k-means-implementation]: https://github.com/Behrouz-Babaki/COP-Kmeans
[silhouette-implementation]: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_samples.html
