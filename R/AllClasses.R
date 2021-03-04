#' @title
#' \code{HiCDOCDataSet} S4 class.
#'
#' @aliases
#' HiCDOCDataSet HiCDOCDataSet-class defaultHiCDOCParameters chromosomes
#' conditions replicates interactions positions compartments concordances
#' differences centroids show
#'
#' @description
#' Data structure for a Hi-C experiment.
#'
#' @details
#' An instance of \code{HiCDOCDataSet} describes a Hi-C experiment with slots
#' for path(s) to input file(s), interactions, pipeline parameters defaulting to
#' \code{defaultHiCDOCParameters}, and computation results. It can be
#' constructed from 4 different types of data:
#' - Tabular files: see \code{\link{HiCDOCDataSetFromTabular}}
#' - (m)Cool files: see \code{\link{HiCDOCDataSetFromCool}}
#' - HiC files: see \code{\link{HiCDOCDataSetFromHiC}}
#' - HiC-Pro matrices and bed files: see \code{\link{HiCDOCDataSetFromHiCPro}}
#' An example \code{HiCDOCDataSet} is also available, see
#' \code{\link{HiCDOCDataSetExample}}.
#'
#' @slot input
#' A vector of path(s) to input file(s).
#' @slot parameters
#' A list of parameters for filtering and computation.
#' @slot interactions
#' A tibble of interactions.
#' @slot chromosomes
#' A vector of names of chromosomes.
#' @slot conditions
#' A vector of names of conditions repeated along the replicates.
#' @slot replicates
#' A vector of names of replicates repeated along the conditions.
#' @slot positions
#' A tibble of the position of each bin.
#' @slot binSize
#' The resolution (number of bases per bin).
#' @slot totalBins
#' A list of the number of bins in each chromosome.
#' @slot weakBins
#' A list of weak bins that are filtered out in each chromosome.
#' @slot sparseConditions
#' A list of sparse conditions repeated along the sparse replicates in each
#' chromosome.
#' @slot sparseReplicates
#' A list of sparse replicates repeated along the sparse conditions in each
#' chromosome.
#' @slot compartments
#' A tibble of the A or B compartment of each bin in each condition.
#' @slot concordances
#' A tibble of the concordance of each bin in each replicate.
#' @slot differences
#' A tibble of detected compartment differences between conditions.
#' @slot distances
#' A tibble of the distances to centroids of each bin in each replicate.
#' @slot centroids
#' A tibble of centroids in each chromosome and condition.
#' @slot selfInteractionRatios
#' A tibble of differences between self interaction and other interactions for
#' each bin in each replicate.
#'
#' @seealso
#' \code{\link{HiCDOCExample}},
#' \code{\link{HiCDOCDataSetFromTabular}},
#' \code{\link{HiCDOCDataSetFromCool}},
#' \code{\link{HiCDOCDataSetFromHiC}},
#' \code{\link{HiCDOCDataSetFromHiCPro}}
#'
#' @md
#' @export
setClass(
    "HiCDOCDataSet",
    slots = c(
        input = "ANY",
        parameters = "ANY",
        interactions = "ANY",
        chromosomes = "ANY",
        conditions = "ANY",
        replicates = "ANY",
        positions = "ANY",
        binSize = "ANY",
        totalBins = "ANY",
        weakBins = "ANY",
        sparseConditions = "ANY",
        sparseReplicates = "ANY",
        compartments = "ANY",
        concordances = "ANY",
        differences = "ANY",
        distances = "ANY",
        centroids = "ANY",
        selfInteractionRatios = "ANY"
    )
)

#' @describeIn HiCDOCDataSet
#' Provides default parameters, imported into a \code{HiCDOCDataSet} upon
#' creation. The \code{HiCDOCDataSet} parameters are then used for the
#' \code{\link{HiCDOC}} pipeline.
#'
#' @format
#'
#' @usage
#'
#' @export
defaultHiCDOCParameters <- list(
    smallChromosomeThreshold = 100,
    weakPositionThreshold = 1,
    sparseReplicateThreshold = 0.05,
    loessSampleSize = 20000,
    kMeansDelta = 0.0001,
    kMeansIterations = 50,
    kMeansRestarts = 20
)

#' @title
#' \code{\link{HiCDOCDataSet}} constructor from a tabular file.
#'
#' @description
#' Constructs a \code{\link{HiCDOCDataSet}} from a tabular file.
#'
#' @details
#' Accepts a tabular file with \code{chromosome}, \code{posiiton 1},
#' \code{position 2}, and multiple replicate columns listing interaction counts.
#' Null interactions do not have to be listed. Values must be separated by
#' tabulations. The header must be
#' \code{chromosome	position 1	position 2	x.y	x.y	x.y	...} with \code{x}
#' replaced by condition names and \code{y} replaced by replicate names.
#'
#' @param path
#' A path to a tabular file.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @examples
#' path <- system.file("extdata", "sample.tsv", package = "HiCDOC")
#' dataSet <- HiCDOCDataSetFromTabular(path)
#'
#' @usage
#' HiCDOCDataSetFromTabular(path)
#'
#' @export
HiCDOCDataSetFromTabular <- function(path = NULL) {

    if (!is.character(path) || length(path) > 1) {
        stop("'paths' must be a string of characters.", call. = FALSE)
    }
    if (!file.exists(path)) {
        stop("'", path, "' does not exist.", call. = FALSE)
    }

    object <- new("HiCDOCDataSet")
    object@input <- path
    object <- .parseTabular(object)
    object <- .fillHiCDOCDataSet(object)
    return(invisible(object))
}

#' @title
#' \code{\link{HiCDOCDataSet}} constructor from Cool files.
#'
#' @description
#' Constructs a \code{\link{HiCDOCDataSet}} from a set of \code{.cool} or
#' \code{.mcool} files.
#'
#' @param paths
#' A vector of paths to \code{.cool} or \code{.mcool} files.
#' @param replicates
#' A vector of replicate names repeated along the conditions.
#' @param conditions
#' A vector of condition names repeated along the replicates.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @examples
#' # Retrieve Cool files
#' directory <- system.file("extdata", package = "HiCDOC")
#' paths <- list.files(directory, '\\.cool$')[0:5]
#' # Specify replicate and condition for each file
#' # In this case:
#' # - The first file in the paths vector is replicate 1 of condition a
#' # - The second file in the paths vector is replicate 2 of condition a
#' # - ...
#' # - The last file in the paths vector is replicate x of condition 3
#' replicates <- (  1,   2,   1,   2, 'x')
#' conditions <- ('a', 'a',   2,   2,   3)
#' # Create dataSet
#' dataSet <- HiCDOCDataSetFromCool(
#'     paths,
#'     replicates = replicates
#'     conditions = conditions
#' )
#'
#' @usage
#' HiCDOCDataSetFromCool(paths, replicates, conditions)
#'
#' @export
HiCDOCDataSetFromCool <- function(
    paths = NULL,
    replicates = NULL,
    conditions = NULL
) {

    if (is.factor(paths)) paths <- as.vector(paths)
    if (!is.character(paths)) {
        stop("'paths' must be a vector of characters.", call. = FALSE)
    }
    for (path in paths) {
        # Remove trailing URI in case of an mCool file
        path <- strsplit(path, "::/")[[1]][1]
        if (!file.exists(path)) {
            stop("'", path, "' does not exist.", call. = FALSE)
        }
    }

    if (is.factor(replicates)) conditions <- as.vector(replicates)
    if (is.null(replicates)) {
        stop("'replicates' must be a vector of replicates.", call. = FALSE)
    }

    if (is.factor(conditions)) conditions <- as.vector(conditions)
    if (is.null(conditions)) {
        stop("'conditions' must be a vector of conditions.", call. = FALSE)
    }

    object <- new("HiCDOCDataSet")
    object@input <- paths
    object@replicates <- replicates
    object@conditions <- conditions
    object <- .parseCool(object)
    object <- .fillHiCDOCDataSet(object)
    return(invisible(object))
}

#' @title
#' \code{\link{HiCDOCDataSet}} constructor from HiC files.
#'
#' @description
#' Constructs a \code{\link{HiCDOCDataSet}} from a set of
#' \code{.hic} files.
#'
#' @param paths
#' A vector of paths to \code{.hic} files.
#' @param replicates
#' A vector of replicate names repeated along the conditions.
#' @param conditions
#' A vector of condition names repeated along the replicates.
#' @param resolution
#' The resolution (span of each position in number of bases) to select within
#' the \code{.hic} files.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @examples
#' # Retrieve Hi-C files
#' directory <- system.file("extdata", package = "HiCDOC")
#' paths <- list.files(directory, '\\.hic$')[0:5]
#' # Specify replicate and condition for each file
#' # In this case:
#' # - The first file in the paths vector is replicate 1 of condition a
#' # - The second file in the paths vector is replicate 2 of condition a
#' # - ...
#' # - The last file in the paths vector is replicate x of condition 3
#' replicates <- (  1,   2,   1,   2, 'x')
#' conditions <- ('a', 'a',   2,   2,   3)
#' # Create dataSet
#' dataSet <- HiCDOCDataSetFromHiC(
#'     paths,
#'     replicates = replicates
#'     conditions = conditions,
#'     resolution = 500000
#' )
#'
#' @usage
#' HiCDOCDataSetFromHiC(paths, replicates, conditions, resolution)
#'
#' @export
HiCDOCDataSetFromHiC <- function(
    paths = NULL,
    replicates = NULL,
    conditions = NULL,
    resolution = NULL
) {

    if (is.factor(paths)) paths <- as.vector(paths)
    if (!is.character(paths)) {
        stop("'paths' must be a vector of characters.", call. = FALSE)
    }
    for (path in paths) {
        if (!file.exists(path)) {
            stop("'", path, "' does not exist.", call. = FALSE)
        }
    }

    if (is.factor(replicates)) replicates <- as.vector(replicates)
    if (is.null(replicates)) {
        stop("'replicates' must be a vector of replicates.", call. = FALSE)
    }

    if (is.factor(conditions)) conditions <- as.vector(conditions)
    if (is.null(conditions)) {
        stop("'conditions' must be a vector of conditions.", call. = FALSE)
    }

    if (!is.numeric(resolution) || length(resolution) != 1) {
        stop("'resolution' must be an integer.", call. = FALSE)
    }

    object <- new("HiCDOCDataSet")
    object@input <- paths
    object@replicates <- replicates
    object@conditions <- conditions
    object@binSize <- resolution
    object <- .parseHiC(object)
    object <- .fillHiCDOCDataSet(object)
    return(invisible(object))
}

#' @title
#' \code{\link{HiCDOCDataSet}} constructor from HiC-Pro files.
#'
#' @description
#' Constructs a \code{\link{HiCDOCDataSet}} from a set of HiC-Pro generated
#' files.
#'
#' @param matrixPaths
#' A vector of paths to HiC-Pro matrix files.
#' @param bedPaths
#' A vector of paths to HiC-Pro bed files.
#' @param replicates
#' A vector of replicate names repeated along the conditions.
#' @param conditions
#' A vector of condition names repeated along the replicates.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @examples
#' # Retrieve HiC-Pro files
#' directory <- system.file("extdata", package = "HiCDOC")
#' matrixPaths <- list.files(directory, '\\.hicpro$')[0:5]
#' bedPaths <- list.files(directory, '\\.bed$')[0:5]
#' # Specify replicate and condition for each pair of files
#' # In this case:
#' # - The first files in the paths vectors are replicate 1 of condition a
#' # - The second files in the paths vectors are replicate 2 of condition a
#' # - ...
#' # - The last files in the paths vectors are replicate x of condition 3
#' replicates <- (  1,   2,   1,   2, 'x')
#' conditions <- ('a', 'a',   2,   2,   3)
#' # Create dataSet
#' dataSet <- HiCDOCDataSetFromHiCPro(
#'     matrixPaths,
#'     bedPaths,
#'     replicates = replicates
#'     conditions = conditions
#' )
#'
#' @usage
#' HiCDOCDataSetFromHiCPro(matrixPaths, bedPaths, replicates, conditions)
#'
#' @export
HiCDOCDataSetFromHiCPro <- function(
    matrixPaths = NULL,
    bedPaths = NULL,
    replicates = NULL,
    conditions = NULL
) {

    if (is.factor(matrixPaths)) matrixPaths <- as.vector(matrixPaths)
    if (!is.character(matrixPaths)) {
        stop("'matrixPaths' must be a vector of characters.", call. = FALSE)
    }

    if (is.factor(bedPaths)) bedPaths <- as.vector(bedPaths)
    if (!is.character(bedPaths)) {
        stop("'bedPaths' must be a vector of characters.", call. = FALSE)
    }

    if (length(matrixPaths) != length(bedPaths)) {
        stop(
            "'matrixPaths' and 'bedPaths' must have the same length.",
            call. = FALSE
        )
    }

    paths <-
      split(
          cbind(matrixPaths, bedPaths),
          seq(length(matrixPaths))
      )

    for (path in unlist(paths)) {
        if (!file.exists(path)) {
            stop("'", path, "' does not exist.", call. = FALSE)
        }
    }

    if (is.factor(replicates)) replicates <- as.vector(replicates)
    if (is.null(replicates)) {
        stop("'replicates' must be a vector of replicates.", call. = FALSE)
    }

    if (is.factor(conditions)) conditions <- as.vector(conditions)
    if (is.null(conditions)) {
        stop("'conditions' must be a vector of conditions.", call. = FALSE)
    }

    object <- new("HiCDOCDataSet")
    object@input <- paths
    object@replicates <- replicates
    object@conditions <- conditions
    parsed <- .parseHiCPro(object)
    object <- parsed[["matrix"]]
    object@binSize <- parsed[["resolution"]]
    object@positions <- parsed[["positions"]]
    object <- .fillHiCDOCDataSet(object)
    return(invisible(object))
}

#' @title
#' Example \code{\link{HiCDOCDataSet}} constructor.
#'
#' @description
#' Constructs a toy \code{\link{HiCDOCDataSet}} from simulated data.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @examples
#' dataSet <- HiCDOCExample()
#'
#' @export
HiCDOCExample <- function() {
    object <- NULL
    basedir <- system.file("extdata", package = "HiCDOC", mustWork = TRUE)
    path <- file.path(basedir, "sampleMatrix.tsv")
    object <- HiCDOCDataSetFromTabular(path)
    return(invisible(object))
}

#' @title
#' Default pipeline to run on a \code{\link{HiCDOCDataSet}}.
#'
#' @description
#' Runs the default filtering, normalization, and computational steps on a
#' \code{\link{HiCDOCDataSet}}. To learn more about HiCDOC, browse the vignette:
#' \code{browseVignettes(package = "HiCDOC")}.
#'
#' @details
#' \subsection{\code{HiCDOCDataSet} object}{
#'     A \code{\link{HiCDOCDataSet}} object can be created from different files:
#'     \itemize{
#'         \item{Tabular files: }{
#'             see \code{\link{HiCDOCDataSetFromTabular}}
#'         }
#'         \item{(m)Cool files: }{
#'             see \code{\link{HiCDOCDataSetFromCool}}
#'         }
#'         \item{HiC files: }{
#'             see \code{\link{HiCDOCDataSetFromHiC}}
#'         }
#'         \item{HiC-Pro matrices and bed files: }{
#'             see \code{\link{HiCDOCDataSetFromHiCPro}}
#'         }
#'     }
#' }
#' \subsection{\code{HiCDOC} pipeline}{
#'     The HiCDOC pipeline has seven steps:
#'     \describe{
#'         \item{Three filtering steps:}{
#'             \itemize{
#'                 \item{\code{\link{filterSmallChromosomes}}}{
#'                     to filter out small chromosomes
#'                 }
#'                 \item{\code{\link{filterWeakPositions}}}{
#'                     to filter out weak positions with very few interactions
#'                 }
#'                 \item{\code{\link{filterSparseReplicates}}}{
#'                     to filter out sparse replicates with many null
#'                     interactions
#'                 }
#'             }
#'         }
#'         \item{Three normalization steps:}{
#'             \itemize{
#'                 \item{\code{\link{normalizeTechnicalBiases}}}{
#'                     to normalize technical biases in each replicates
#'                 }
#'                 \item{\code{\link{normalizeBiologicalBiases}}}{
#'                     to normalize biological biases in each replicate
#'                 }
#'                 \item{\code{\link{normalizeDistanceEffect}}}{
#'                     to normalize the distance effect in each chromosome
#'                 }
#'             }
#'         }
#'         \item{One computational step:}{
#'             \itemize{
#'                 \item{\code{\link{detectCompartments}}}{
#'                     to detect compartments in each condition and find
#'                     significant changes between conditions.
#'                 }
#'             }
#'         }
#'     }
#' }
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param parallel
#' Whether or not to parallelize each step. Defaults to TRUE.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}} with all slots filled.
#'
#' @seealso
#' \code{\link{HiCDOCDataSet}}
#'
#' @examples
#' object <- HiCDOCExample()
#' object <- HiCDOC(object)
#'
#' # Equivalent to
#' \dontrun{
#' object <- filterSmallChromosomes(object)
#' object <- filterSparseReplicates(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' }
#'
#' @export
HiCDOC <- function(object, parallel = TRUE) {
    object <- filterSmallChromosomes(object)
    object <- filterSparseReplicates(object)
    object <- filterWeakPositions(object)
    object <- normalizeTechnicalBiases(object, parallel = parallel)
    object <- normalizeBiologicalBiases(object)
    object <- normalizeDistanceEffect(object)
    object <- detectCompartments(object, parallel = parallel)
    return(invisible(object))
}
