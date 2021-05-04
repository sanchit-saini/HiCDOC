# HiCDOC: Compartments prediction and differential analysis with multiple replicates

HiCDOC normalizes intrachromosomal Hi-C matrices, uses unsupervised learning to
predict A/B compartments from multiple replicates, and detects significant
compartment changes between experiment conditions.

It provides a collection of functions assembled into a pipeline:

1. [Filter](#filtering-data):
    1. Remove chromosomes which are too small to be useful.
    2. Filter sparse replicates to remove uninformative replicates with few
       interactions.
    3. Filter positions (*bins*) which have too few interactions.
2. [Normalize](#normalizing-biases):
    1. Normalize technical biases using
       [cyclic loess normalization][multihiccompare-publication], so that
       matrices are comparable.
    2. Normalize biological biases using
       [Knight-Ruiz matrix balancing][knight-ruiz-publication], so that
       all the bins are comparable.
    3. Normalize the distance effect, which results from higher interaction
       proportions between closer regions, with a MD loess.
3. [Predict](#predicting-compartments-and-differences):
    1. Predict compartments using
       [constrained K-means][constrained-k-means-publication].
    2. Detect significant differences between experiment conditions.
4. [Visualize](#visualizing-data-and-results):
    1. Plot the interaction matrices of each replicate.
    2. Plot the overall distance effect on the proportion of interactions.
    3. Plot the compartments in each chromosome, along with their concordance
       (confidence measure) in each replicate, and significant changes between
       experiment conditions.
    4. Plot the overall distribution of concordance differences.
    5. Plot the result of the PCA on the compartments' centroids.
    6. Plot the boxplots of self interaction ratios (differences between self
       interactions and the medians of other interactions) of each compartment,
       which is used for the A/B classification.


# Table of contents

* [Installation](#installation)
* [Quick Start](#quick-start)
* [Usage](#quick-start)
    * [Importing Hi-C data](#importing-hi-c-data)
        * [Tabular files](#tabular-files)
        * [Cooler files](#cooler-files)
        * [Juicer files](#juicer-files)
        * [HiC-Pro files](#hic-pro-files)
    * [Running the HiCDOC pipeline](#running-the-hicdoc-pipeline)
        * [Filtering data](#filtering-data)
        * [Normalizing biases](#normalizing-biases)
        * [Predicting compartments and differences](#predicting-compartments-and-differences)
    * [Visualizing data and results](#visualizing-data-and-results)
* [References](#references)


# Installation

To install, execute the following commands in your console:

```bash
Rscript -e 'install.packages("devtools")'
Rscript -e 'devtools::install_github("mzytnicki/HiCDOC")'
```

After installation, the package can be loaded in R >= 4.0:

```r
library("HiCDOC")
```


# Quick Start

To try out HiCDOC, load the simulated toy data set:

```r
data(exampleHiCDOCDataSet)
hic.experiment <- exampleHiCDOCDataSet
```

Then run the default pipeline on the created object:

```r
hic.experiment <- HiCDOC(hic.experiment)
```

And plot some results:

```r
plotCompartmentChanges(hic.experiment, chromosome = 'Y')
```


# Usage

## Importing Hi-C data

HiCDOC can import Hi-C data sets in various different formats:
- Tabular `.tsv` files.
- Cooler `.cool` or `.mcool` files.
- Juicer `.hic` files.
- HiC-Pro `.matrix` and `.bed` files.

### Tabular files

A tabular file is a tab-separated multi-replicate sparse matrix with a header:

    chromosome    position 1    position 2    C1.R1    C1.R2    C2.R1    ...
    3             1500000       7500000       145      184      72       ...
    ...

The interaction proportions between `position 1` and `position 2` of
`chromosome` are reported in each `condition.replicate` column. There is no
limit to the number of conditions and replicates.

To load Hi-C data in this format:

```r
hic.experiment <- HiCDOCDataSetFromTabular('path/to/data.tsv')
```

### Cooler files

To load `.cool` or `.mcool` files generated by [Cooler][cooler-documentation]:

```r
# Path to each file
paths = c(
  'path/to/condition-1.replicate-1.cool',
  'path/to/condition-1.replicate-2.cool',
  'path/to/condition-2.replicate-1.cool',
  'path/to/condition-2.replicate-2.cool',
  'path/to/condition-3.replicate-1.cool'
)

# Replicate and condition of each file. Can be names instead of numbers.
replicates <- c(1, 2, 1, 2, 1)
conditions <- c(1, 1, 2, 2, 3)

# Resolution to select in .mcool files
binSize = 500000

# Instantiation of data set
hic.experiment <- HiCDOCDataSetFromCool(
  paths,
  replicates = replicates,
  conditions = conditions,
  binSize = binSize # Specified for .mcool files.
)
```

### Juicer files

To load `.hic` files generated by [Juicer][juicer-documentation]:

```r
# Path to each file
paths = c(
  'path/to/condition-1.replicate-1.hic',
  'path/to/condition-1.replicate-2.hic',
  'path/to/condition-2.replicate-1.hic',
  'path/to/condition-2.replicate-2.hic',
  'path/to/condition-3.replicate-1.hic'
)

# Replicate and condition of each file. Can be names instead of numbers.
replicates <- c(1, 2, 1, 2, 1)
conditions <- c(1, 1, 2, 2, 3)

# Resolution to select
binSize <- 500000

# Instantiation of data set
hic.experiment <- HiCDOCDataSetFromHiC(
  paths,
  replicates = replicates,
  conditions = conditions,
  binSize = binSize
)
```

### HiC-Pro files

To load `.matrix` and `.bed` files generated by [HiC-Pro][hicpro-documentation]:

```r
# Path to each matrix file
matrixPaths = c(
  'path/to/condition-1.replicate-1.matrix',
  'path/to/condition-1.replicate-2.matrix',
  'path/to/condition-2.replicate-1.matrix',
  'path/to/condition-2.replicate-2.matrix',
  'path/to/condition-3.replicate-1.matrix'
)

# Path to each bed file
bedPaths = c(
  'path/to/condition-1.replicate-1.bed',
  'path/to/condition-1.replicate-2.bed',
  'path/to/condition-2.replicate-1.bed',
  'path/to/condition-2.replicate-2.bed',
  'path/to/condition-3.replicate-1.bed'
)

# Replicate and condition of each file. Can be names instead of numbers.
replicates <- c(1, 2, 1, 2, 1)
conditions <- c(1, 1, 2, 2, 3)

# Instantiation of data set
hic.experiment <- HiCDOCDataSetFromHiCPro(
  matrixPaths = matrixPaths,
  bedPaths = bedPaths,
  replicates = replicates,
  conditions = conditions
)
```

## Running the HiCDOC pipeline

Once your data is loaded, you can run all the filtering, normalization, and
prediction steps with:

```r
hic.experiment <- HiCDOC(hic.experiment)
```

This one-liner runs all the steps detailed below.

### Filtering data

Remove small chromosomes of length smaller than 100 positions:

```r
hic.experiment <- filterSmallChromosomes(hic.experiment, threshold = 100)
```

Remove sparse replicates filled with less than 5% non-zero interactions:

```r
hic.experiment <- filterSparseReplicates(hic.experiment, threshold = 0.05)
```

Remove weak positions with less than 1 interaction in average:

```r
hic.experiment <- filterWeakPositions(hic.experiment, threshold = 1)
```

### Normalizing biases

Normalize technical biases such as sequencing depth:

```r
hic.experiment <- normalizeTechnicalBiases(hic.experiment)
```

Normalize biological biases (such as GC content, number of restriction sites,
etc.):

```r
hic.experiment <- normalizeBiologicalBiases(hic.experiment)
```

Normalize the distance effect resulting from higher interaction proportions
between closer regions:

```r
hic.experiment <- normalizeDistanceEffect(hic.experiment, loessSampleSize = 20000)
```

### Predicting compartments and differences

Predict A and B compartments and detect significant differences:

```r
hic.experiment <- detectCompartments(
  hic.experiment,
  kMeansDelta = 0.0001,
  kMeansIterations = 50,
  kMeansRestarts = 20
)
```


## Visualizing data and results

Plot the interaction matrix of each replicate:

```r
plotInteractions(hic.experiment, chromosome = '3')
```

Plot the overall distance effect on the proportion of interactions:

```r
plotDistanceEffect(hic.experiment)
```

List and plot compartments with their concordance (confidence measure) in each
replicate, and significant changes between experiment conditions:

```r
compartments(hic.experiment)

concordances(hic.experiment)

differences(hic.experiment)

plotCompartmentChanges(hic.experiment, chromosome = '3')
```

Plot the overall distribution of concordance differences:

```r
plotConcordanceDifferences(hic.experiment)
```

Plot the result of the PCA on the compartments' centroids:

```r
plotCentroids(hic.experiment, chromosome = '3')
```

Plot the boxplots of self interaction ratios (differences between self
interactions and the median of other interactions) of each compartment:

```r
plotSelfInteractionRatios(hic.experiment, chromosome = '3')
```


# References

John C Stansfield, Kellen G Cresswell, Mikhail G Dozmorov, multiHiCcompare:
joint normalization and comparative analysis of complex Hi-C experiments,
_Bioinformatics_, 2019, https://doi.org/10.1093/bioinformatics/btz048

Philip A. Knight, Daniel Ruiz, A fast algorithm for matrix balancing, _IMA
Journal of Numerical Analysis_, Volume 33, Issue 3, July 2013, Pages 1029–1047,
https://doi.org/10.1093/imanum/drs019

Kiri Wagstaff, Claire Cardie, Seth Rogers, Stefan Schrödl, Constrained K-means
Clustering with Background Knowledge, _Proceedings of 18th International
Conference on Machine Learning_, 2001, Pages 577-584,
https://pdfs.semanticscholar.org/0bac/ca0993a3f51649a6bb8dbb093fc8d8481ad4.pdf

[multihiccompare-publication]: https://doi.org/10.1093/bioinformatics/btz048
[knight-ruiz-publication]: https://doi.org/10.1093/imanum/drs019
[constrained-k-means-publication]: https://pdfs.semanticscholar.org/0bac/ca0993a3f51649a6bb8dbb093fc8d8481ad4.pdf

[cooler-documentation]: https://cooler.readthedocs.io/en/latest/
[juicer-documentation]: https://github.com/aidenlab/juicer/wiki/Data
[hicpro-documentation]: https://github.com/nservant/HiC-Pro
