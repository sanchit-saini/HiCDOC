#### - HiCDOCDefaultParameters -------------------------------------------####
# @describeIn parameters Default parameters for the HiCDOC pipeline
#' @rdname parameters
#' @examples
#' # All the default parameters
#' HiCDOCDefaultParameters
#'
#' @export
HiCDOCDefaultParameters <- list(
    smallChromosomeThreshold = 100,
    weakPositionThreshold = 0,
    sparseReplicateThreshold = 0.95,
    loessSampleSize = 20000,
    kMeansIterations = 50,
    kMeansDelta = 0.0001,
    kMeansRestarts = 20
)


## - HiCDOCDataSet S4 class constructor from sparse matrix ------------------#
## --------------------------------------------------------------------------#
#' Reads a sparse matrix and fill a \code{HiCDOCDataSet} with its content.
#'
#' @rdname HiCDOCDataSetFromTabular
#' @docType class
#' @param matrix A matrix with the data.
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
#' @examples
#' link <- system.file("extdata", "sampleMatrix.tsv", package = "HiCDOC")
#' dataSet <- HiCDOCDataSetFromTabular(link)
#' dataSet
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
    object <- parseTabular(object)
    object <- HiCDOCDataSet(object)
    return(invisible(object))
}


## - HiCDOCDataSet S4 class constructor from cool files ---------------------#
## --------------------------------------------------------------------------#
#' Construct a \code{HiCDOCDataSet} from a cool file.
#' @rdname HiCDOCDataSetFromCool
#' @docType class
#'
#' @param paths A vector of file names, each one being a .cool file.
#' @param replicates   character vector. The list of the replicate names.
#' This vector can have repeated values if the replicates are repeated
#' over the conditions.
#' @param conditions   character vector. The names of the two conditions.
#' This should be repeated values of the names, the vector should have the
#' same size as \code{replicates}.
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
#'
#' @examples
#' basedir <- system.file("extdata", package = "HiCDOC")
#' data <- read.csv(file.path(basedir, "coolData.csv"))
#' data
#' dataSet <- HiCDOCDataSetFromCool(
#'     file.path(basedir, data$FileName),
#'     data$Replicate,
#'     data$Condition
#' )
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

    if (is.factor(conditions)) conditions <- as.vector(conditions)
    if (is.null(conditions)) {
        stop("'conditions' must be a vector of conditions.", call. = FALSE)
    }

    object <- new("HiCDOCDataSet")
    object@input <- paths
    object@replicates <- replicates
    object@conditions <- conditions
    object <- parseCool(object)
    object <- HiCDOCDataSet(object)

    return(invisible(object))
}


## - HiCDOCDataSet S4 class constructor from hic files ----------------------#
## --------------------------------------------------------------------------#
#' Construct a \code{HiCDOCDataSet} from a hic file.
#' @rdname HiCDOCDataSetFromHiC
#' @docType class
#'
#' @param paths A vector of file names, each one being a .hic file.
#' @param replicates   character vector. The list of the replicate names.
#' This vector can have repeated values if the replicates are repeated
#' over the conditions.
#' @param conditions   character vector. The names of the two conditions.
#' This should be repeated values of the names, the vector should have the
#' same size as \code{replicates}.
#' @param resolution   integer. The resolution of the matrices.
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
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
    object <- parseHiC(object)
    object <- HiCDOCDataSet(object)

    return(invisible(object))
}

## - HiCDOCDataSet S4 class constructor from HiC-Pro files ------------------#
## --------------------------------------------------------------------------#
#' Construct a \code{HiCDOCDataSet} from hHiC-Pro files
#' @rdname HiCDOCDataSetFromHiCPro
#' @docType class
#'
#' @param matrixPaths A vector of file names, each one being a .matrix file.
#' @param bedPaths A vector of file names, each one being a .bed file.
#' @param replicates   character vector. The list of the replicate names.
#' This vector can have repeated values if the replicates are repeated
#' over the conditions.
#' @param conditions   character vector. The names of the two conditions.
#' This should be repeated values of the names, the vector should have the
#' same size as \code{replicates}.
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
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

    if (is.factor(conditions)) conditions <- as.vector(conditions)
    if (is.null(conditions)) {
        stop("'conditions' must be a vector of conditions.", call. = FALSE)
    }

    object <- new("HiCDOCDataSet")
    object@input <- paths
    object@replicates <- replicates
    object@conditions <- conditions
    parsed <- parseHiCPro(object)
    object <- parsed[["matrix"]]
    object@binSize <- parsed[["resolution"]]
    object@positions <- parsed[["positions"]]
    object <- HiCDOCDataSet(object)

    return(invisible(object))
}

#### - HiCDOCDataSet constructors-----------------------------------------####

## - Example constructor ----------------------------------------------------#
## --------------------------------------------------------------------------#
#' Example constructor
#'
#' This function provides an example of a \code{HiCDOCDataSet} object from
#' Marti-Marimon et al.
#'
#' @return An \code{HiCDOCDataSet} object called '\code{exp}'.
#'
#' @examples
#' ## The 'HiCDOCDataSet' object in this example was constructed by:
#' exp <- HiCDOCExample()
#' exp
#' @export
HiCDOCExample <- function() {
    object <- NULL
    basedir <- system.file("extdata", package = "HiCDOC", mustWork = TRUE)
    path <- file.path(basedir, "sampleMatrix.tsv")
    object <- HiCDOCDataSetFromTabular(path)
    return(invisible(object))
}

#### - class definition ####
## HiCDOCDataSet S4 class definition ----------------------------------------#
## --------------------------------------------------------------------------#
#' Infrastructure for HiCDOC experiment and differential interaction.
#'
#' \code{HiCDOCDataSet} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @aliases HiCDOCDataSet HiCDOCDataSet-class chromosomes conditions
#' replicates interactions positions differences concordances compartments
#' centroids show
#' @details
#' A HiCDOCDataSet can be constructed from 4 different types of data :
#' * Matrices: see \code{\link{HiCDOCDataSetFromTabular}}
#' * Cool or mCool data: see \code{\link{HiCDOCDataSetFromCool}}
#' * Hic files: see \code{\link{HiCDOCDataSetFromHiC}}
#' * Hic Pro matrices and bed files: see \code{\link{HiCDOCDataSetFromHiCPro}}
#' @slot input The names of the matrix input files.
#' @slot interactions The interaction matrices.
#' @slot weakBins The empty bins.
#' @slot chromosomes The list of chromosomes.
#' @slot replicates The names of the replicates.
#' @slot conditions The names of the conditions (exactly two different).
#' @slot totalBins The number of bins per chromosome.
#' @slot binSize The resolution.
#' @slot distances The distribution of distances to the centroids.
#' @slot selfInteractionRatios TODO
#' @slot compartments The A/B compartments distribution, along the chromosomes.
#' @slot concordances The concordance distribution, along the chromosomes.
#' @slot differences The distribution of the difference of the concordance.
#' @slot centroids The position of the centroids.
#' @slot parameters An named \code{list}. The parameters for the
#' segmentation methods. See \code{\link{parameters}}.
#' @slot positions The position of the bins.
#' @seealso \code{\link{HiCDOCExample}},
#' \code{\link{parameters}},
#' \code{\link{HiCDOCDataSetFromTabular}},
#' \code{\link{HiCDOCDataSetFromCool}},
#' \code{\link{HiCDOCDataSetFromHiC}},
#' \code{\link{HiCDOCDataSetFromHiCPro}}
#' @export
#' @md
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
        weakBins = "ANY",
        totalBins = "ANY",
        compartments = "ANY",
        concordances = "ANY",
        differences = "ANY",
        distances = "ANY",
        centroids = "ANY",
        selfInteractionRatios = "ANY"
    )
)


## - HiCDOCDataSet S4 class constructor -------------------------------------#
## --------------------------------------------------------------------------#
#' Constructor for the HiCDOCDataSet class.
#'
#' This function should not be called directly. It is called by the functions
#' \code{\link{HiCDOCDataSetFromTabular}},
#' \code{\link{HiCDOCDataSetFromCool}},
#' \code{\link{HiCDOCDataSetFromHiC}} and
#' \code{\link{HiCDOCDataSetFromHiCPro}}.
#' @name HiCDOCDataSet-constructor
#' @rdname HiCDOCDataSet-constructor
#' @aliases HiCDOCDataSet-constructor
#' @param object A prefilled \code{HiCDOCDataSet} object.
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
#' @keywords internal
#' @noRd
HiCDOCDataSet <- function(object = NULL) {
    # Create new empty object
    if (is.null(object)) object <- new("HiCDOCDataSet")

    # Fill binSize slot
    if (is.null(object@binSize)) object@binSize <- computeBinSize(object)

    if (!is.null(object@interactions)) {
        chromosomes <-
            gtools::mixedsort(
                unique(as.character(object@interactions$chromosome))
            )
        object@chromosomes <- chromosomes

        # Determine positions
        if (is.null(object@positions)) {
            object@positions <- determinePositions(object)
        }

        # Replace positions by bins
        object@interactions <- replacePositionsByBins(object)

        object@weakBins <- vector("list", length(chromosomes))
        names(object@weakBins) <- chromosomes

        # Make upper triangular matrix, values as numeric, average duplicates
        object@interactions <- reformatInteractions(object@interactions)
    }

    # Fill totalBins slot
    if (!is.null(object@binSize) & !is.null(object@interactions)) {
        object@totalBins <-
            vapply(
                as.character(object@chromosomes),
                function(x) {
                    nrow(object@positions[object@positions$chromosome == x, ])
                },
                FUN.VALUE = 0
            )
    }

    object@parameters <- HiCDOCDefaultParameters
    return(invisible(object))
}

#### - Pipeline method ####
## - Main method ------------------------------------------------------------#
## --------------------------------------------------------------------------#
#' @description
#' To learn more about HiCDOC, start with the vignette:
#' \code{browseVignettes(package = "HiCDOC")}.
#'
#' The HiCDOC function allow to run the HiCDOC pipeline, on a HiCDOCDataSet
#' object.
#'
#' @usage
#' HiCDOC(object, parallel = FALSE) # to run the entire HiCDOC pipeline
#' @details
#' \subsection{HiCDOCDataSet object}{
#'     A \code{HiCDOCDataSet} object can be created from different object :
#'     \itemize{
#'         \item{matrices: }{see \code{\link{HiCDOCDataSetFromTabular}}},
#'         \item{cool or mcool data: }{see \code{\link{HiCDOCDataSetFromCool}}},
#'         \item{hic data: }{see \code{\link{HiCDOCDataSetFromHiC}}} and
#'         \item{HicPro data: }{see \code{\link{HiCDOCDataSetFromHiCPro}}}
#'     }
#' }
#' \subsection{The HiCDOCDataSet pipeline have 7 steps :}{
#'    \describe{
#'        \item{3 filters to clean the data:}{
#'        \itemize{
#'            \item{\code{\link{filterSmallChromosomes}}}{
#'            to filter the too small chromosomes}
#'            \item{\code{\link{filterWeakPositions}}}{
#'            to filter the "weak" positions of a chromosome}
#'            \item{\code{\link{filterSparseChromosomes}}}{
#'            to filter the too sparse interactions matrix
#'            by chromosome}
#'        }
#'      }
#'    \item{3 normalization steps to normalize the data:}{
#'    \itemize{
#'      \item{\code{\link{normalizeTechnicalBiases}}}{
#'      to normalize the data in case of tehnical bias}
#'      \item{\code{\link{normalizeBiologicalBiases}}}{
#'      to normalize the data in case of biological bias}
#'      \item{\code{\link{normalizeDistanceEffect}}}{
#'      to normalize the distance effect}
#'    }
#'    }
#'    \item{the detection of compartments}{
#' \itemize{
#'      \item{\code{\link{detectCompartments}}}{
#'      to run the detection of compartments and compute
#'      significant differences}
#'    }
#'  }
#'    }
#'}
#' @param object A \code{HiCDOCDataSet} object.
#' @param parallel Logical, default to FALSE. Should parallel computation
#' should be used ?
#' @return Returns an \code{HiCDOCDataSet} object.
#' @seealso \code{\link{HiCDOCDataSet}}
#' @examples
#' object <- HiCDOCExample()
#' object <- HiCDOC(object)
#'
#' # This is equivalent of
#' \dontrun{
#' parallel = FALSE
#' object <- filterSmallChromosomes(object)
#' object <- filterSparseChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object, parallel = parallel)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object, parallel = parallel)
#' }
#' @export
HiCDOC <- function(object, parallel = FALSE) {
    object <- filterSmallChromosomes(object)
    object <- filterSparseChromosomes(object)
    object <- filterWeakPositions(object)
    object <- normalizeTechnicalBiases(object, parallel = parallel)
    object <- normalizeBiologicalBiases(object)
    object <- normalizeDistanceEffect(object)
    object <- detectCompartments(object, parallel = parallel)
    return(invisible(object))
}
