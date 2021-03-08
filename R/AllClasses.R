#' @title
#' \code{HiCDOCDataSet} S4 class.
#'
#' @name
#' HiCDOCDataSet-class
#'
#' @aliases
#' HiCDOCDataSet
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
#' A list of parameters used for filtering, normalization, and prediction of
#' compartments.
#' @slot interactions
#' A tibble of interactions.
#' @slot chromosomes
#' A vector of names of chromosomes.
#' @slot conditions
#' A vector of names of conditions, repeated along the replicates.
#' @slot replicates
#' A vector of names of replicates, repeated along the conditions.
#' @slot positions
#' A tibble of positions and their corresponding bin.
#' @slot binSize
#' The computed bin size (span of each bin in number of bases).
#' @slot totalBins
#' A list of the number of bins in each chromosome.
#' @slot weakBins
#' A list of weak bins that are filtered out in each chromosome.
#' @slot validConditions
#' A list of non-sparse valid conditions, repeated along the valid replicates
#' in each chromosome.
#' @slot validReplicates
#' A list of non-sparse valid replicates, repeated along the valid conditions
#' in each
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
#' \code{\link{HiCDOCDataSetExample}},
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
        validConditions = "ANY",
        validReplicates = "ANY",
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
#' instantiation. The \code{HiCDOCDataSet} parameters are then used for the
#' \code{\link{HiCDOC}} pipeline.
#'
#' @format
#'
#' @usage
#'
#' @export
defaultHiCDOCParameters <- list(
    smallChromosomeThreshold = 100,
    sparseReplicateThreshold = 0.05,
    weakPositionThreshold = 1,
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
#' path <- system.file("extdata", "example.tsv", package = "HiCDOC")
#' object <- HiCDOCDataSetFromTabular(path)
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
#' @param resolution
#' The resolution (span of each position in number of bases). Optionally
#' provided to select the appropriate resolution in \code{.mcool} files.
#' Defaults to NULL.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @examples
#' \dontrun{
#' # Path to each file
#' paths = c(
#'   'path/to/condition-1.replicate-1.cool',
#'   'path/to/condition-1.replicate-2.cool',
#'   'path/to/condition-2.replicate-1.cool',
#'   'path/to/condition-2.replicate-2.cool',
#'   'path/to/condition-3.replicate-1.cool'
#' )
#'
#' # Replicate and condition of each file. Can be names instead of numbers.
#' replicates <- c(1, 2, 1, 2, 1)
#' conditions <- c(1, 1, 2, 2, 3)
#'
#' # Resolution to select in .mcool files
#' resolution = 500000
#'
#' # Instantiation of data set
#' object <- HiCDOCDataSetFromCool(
#'   paths,
#'   replicates = replicates,
#'   conditions = conditions,
#'   resolution = resolution # Specified for .mcool files.
#' )
#' }
#'
#' @usage
#' HiCDOCDataSetFromCool(paths, replicates, conditions)
#'
#' @export
HiCDOCDataSetFromCool <- function(
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

    if (is.factor(replicates)) conditions <- as.vector(replicates)
    if (is.null(replicates)) {
        stop("'replicates' must be a vector of replicates.", call. = FALSE)
    }

    if (is.factor(conditions)) conditions <- as.vector(conditions)
    if (is.null(conditions)) {
        stop("'conditions' must be a vector of conditions.", call. = FALSE)
    }

    if (
        !is.null(resolution) &&
        (!is.numeric(resolution) || length(resolution) != 1)
    ) {
        stop("'resolution' must be an integer.", call. = FALSE)
    }

    object <- new("HiCDOCDataSet")
    object@input <- paths
    object@replicates <- replicates
    object@conditions <- conditions
    object <- .parseCool(object, resolution)
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
#' \dontrun{
#' # Path to each file
#' paths = c(
#'   'path/to/condition-1.replicate-1.hic',
#'   'path/to/condition-1.replicate-2.hic',
#'   'path/to/condition-2.replicate-1.hic',
#'   'path/to/condition-2.replicate-2.hic',
#'   'path/to/condition-3.replicate-1.hic'
#' )
#'
#' # Replicate and condition of each file. Can be names instead of numbers.
#' replicates <- c(1, 2, 1, 2, 1)
#' conditions <- c(1, 1, 2, 2, 3)
#'
#' # Resolution to select
#' resolution <- 500000
#'
#' # Instantiation of data set
#' hic.experiment <- HiCDOCDataSetFromHiC(
#'   paths,
#'   replicates = replicates,
#'   conditions = conditions,
#'   resolution = resolution
#' )
#' }
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
    object <- .parseHiC(object, resolution)
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
#' \dontrun{
#' # Path to each matrix file
#' matrixPaths = c(
#'   'path/to/condition-1.replicate-1.matrix',
#'   'path/to/condition-1.replicate-2.matrix',
#'   'path/to/condition-2.replicate-1.matrix',
#'   'path/to/condition-2.replicate-2.matrix',
#'   'path/to/condition-3.replicate-1.matrix'
#' )
#'
#' # Path to each bed file
#' bedPaths = c(
#'   'path/to/condition-1.replicate-1.bed',
#'   'path/to/condition-1.replicate-2.bed',
#'   'path/to/condition-2.replicate-1.bed',
#'   'path/to/condition-2.replicate-2.bed',
#'   'path/to/condition-3.replicate-1.bed'
#' )
#'
#' # Replicate and condition of each file. Can be names instead of numbers.
#' replicates <- c(1, 2, 1, 2, 1)
#' conditions <- c(1, 1, 2, 2, 3)
#'
#' # Instantiation of data set
#' hic.experiment <- HiCDOCDataSetFromHiCPro(
#'   matrixPaths = matrixPaths,
#'   bedPaths = bedPaths,
#'   replicates = replicates,
#'   conditions = conditions
#' )
#' }
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
    object <- .parseHiCPro(object)
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
#' object <- HiCDOCDataSetExample()
#'
#' @export
HiCDOCDataSetExample <- function() {
    object <- NULL
    directory <- system.file("extdata", package = "HiCDOC", mustWork = TRUE)
    path <- file.path(directory, "example.tsv")
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
#' object <- HiCDOCDataSetExample()
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
