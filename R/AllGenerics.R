#' @title
#' Methods to access a \code{\link{HiCDOCDataSet}} components.
#'
#' @docType methods
#'
#' @name
#' HiCDOCDataSet-methods
#'
#' @aliases
#' chromosomes positions conditions replicates binSize interactions
#' compartments concordances differences show
#'
#' @description
#' Retrieve information and results from a \code{\link{HiCDOCDataSet}}.
#' @examples
#' data(exampleHiCDOCDataSet)
#' chromosomes(exampleHiCDOCDataSet)
#' conditions(exampleHiCDOCDataSet)
#' binSize(exampleHiCDOCDataSet)
#' @return
#' A character vector (for \code{chromosomes}, \code{conditions},
#' \code{replicates}), an integer(for \code{binSize}), a tibble
#' (for \code{interactions} and \code{positions}), or a GRanges object
#' (for \code{compartments}, \code{concordances}, \code{differences}).
NULL

#### chromosomes ####
#' @describeIn HiCDOCDataSet-methods
#' Retrieves the vector of chromosome names.
#' @export
setGeneric(
    name = "chromosomes",
    def = function(object) {
        standardGeneric("chromosomes")
    }
)

#### positions ####
#' @describeIn HiCDOCDataSet-methods
#' Retrieves the genomic positions corresponding to bins for each chromosome.
#' @export
setGeneric(
    name = "positions",
    def = function(object) {
        standardGeneric("positions")
    }
)

#### conditions ####
#' @describeIn HiCDOCDataSet-methods
#' Retrieves the vector of condition names.
#' @export
setGeneric(
    name = "conditions",
    def = function(object) {
        standardGeneric("conditions")
    }
)

#### replicates ####
#' @describeIn HiCDOCDataSet-methods
#' Retrieves the vector of replicate names.
#' @export
setGeneric(
    name = "replicates",
    def = function(object) {
        standardGeneric("replicates")
    }
)

#### interactions ####
#' @describeIn HiCDOCDataSet-methods
#' Retrieves a tibble of the interactions.
#' @export
setGeneric(
    name = "interactions",
    def = function(object) {
        standardGeneric("interactions")
    }
)

#### binSize ####
#' @describeIn HiCDOCDataSet-methods
#' Retrieves the resolution (span of each position in number of bases).
#' @export
setGeneric(
    name = "binSize",
    def = function(object) {
        standardGeneric("binSize")
    }
)

#### compartments ####
#' @describeIn HiCDOCDataSet-methods
#' Retrieves a \code{GenomicRange} of the compartment of every position
#' in every condition.
#' @export
setGeneric(
    name = "compartments",
    def = function(object) {
        standardGeneric("compartments")
    }
)

#### differences ####
#' @describeIn HiCDOCDataSet-methods
#' Retrieves a \code{GenomicRange} of the significant compartment differences
#' between conditions, and their p-values.
#' @usage
#' differences(object, threshold = 0.05)
#' @param object
#' a HiCDOCDataSet object
#' @param threshold
#' a numeric value between 0 and 1. If no threshold, all the differences will
#' be printed even the non significant ones. Otherwise the differences printed
#' are filtered to show the ones with an adjusted p-value <= \code{threshold}.
#' @export
setGeneric(
    name = "differences",
    def = function(object, threshold = 0.05) {
        standardGeneric("differences")
    }
)

#### concordances ####
#' @describeIn HiCDOCDataSet-methods
#' Retrieves a \code{GenomicRange} of the concordance (confidence in assigned
#' compartment) of every position in every replicate.
#' @export
setGeneric(
    name = "concordances",
    def = function(object) {
        standardGeneric("concordances")
    }
)

#' @title
#' Access the parameters of a \code{\link{HiCDOCDataSet}}.
#'
#' @docType methods
#'
#' @name
#' HiCDOCDataSet-parameters
#'
#' @aliases
#' parameters parameters<- defaultHiCDOCParameters
#'
#' @description
#' Retrieves or sets parameters used for filtering, normalization, and
#' prediciton of compartments.
#'
#' @details
#' A \code{\link{HiCDOCDataSet}}'s parameters are automatically set to default
#' values retrieved from \code{\link{defaultHiCDOCParameters}}. They are
#' accessed by filtering, normalization, and compartment detection functions.
#' If those functions are called with custom arguments, the object's
#' parameters are updated to record the actual parameters used. If the
#' object's parameters are customized before calling the functions, the
#' custom parameters will be used.
#'
#' \subsection{All parameters are listed here:}{
#'     \describe{
#'         \item{\code{smallChromosomeThreshold}}{
#'             The minimum length (number of positions) for a chromosome to be
#'             kept when filtering with \code{\link{filterSmallChromosomes}}.
#'             Defaults to
#'             \code{defaultHiCDOCParameters$smallChromosomeThreshold} = 100.
#'         }
#'         \item{\code{sparseReplicateThreshold}}{
#'             The minimum percentage of non-zero interactions for a chromosome
#'             replicate to be kept when filtering with
#'             \code{\link{filterSparseReplicates}}. If a chromosome replicate's
#'             percentage of non-zero interactions is lower than this value, it
#'             is removed. Defaults to
#'             \code{defaultHiCDOCParameters$smallChromosomeThreshold} = 0.05.
#'         }
#'         \item{\code{weakPositionThreshold}}{
#'             The minimum average interaction for a position to be kept when
#'             filtering with \code{\link{filterWeakPositions}}. If a position's
#'             average interaction with the entire chromosome is lower than this
#'             value in any of the replicates, it is removed from all replicates
#'             and conditions. Defaults to
#'             \code{defaultHiCDOCParameters$smallChromosomeThreshold} = 1.
#'         }
#'         \item{\code{loessSampleSize}}{
#'             The number of positions used as a sample to estimate the effect
#'             of distance on proportion of interactions when normalizing with
#'             \code{\link{normalizeDistanceEffect}} Defaults to
#'             \code{defaultHiCDOCParameters$loessSampleSize} = 20000.
#'         }
#'         \item{\code{kMeansDelta}}{
#'             The convergence stop criterion for the clustering when detecting
#'             compartments with \code{\link{detectCompartments}}. When the
#'             centroids' distances between two iterations is lower than this
#'             value, the clustering stops. Defaults to
#'             \code{defaultHiCDOCParameters$kMeansDelta} = 0.0001.
#'         }
#'         \item{\code{kMeansIterations}}{
#'             The maximum number of iterations during clustering when detecting
#'             compartments with \code{\link{detectCompartments}}. Defaults to
#'             \code{defaultHiCDOCParameters$kMeansIterations} = 50.
#'         }
#'         \item{\code{kMeansRestarts}}{
#'             The amount of times the clustering is restarted when detecting
#'             compartments with \code{\link{detectCompartments}}. For each
#'             restart, the clustering iterates until convergence or reaching
#'             the maximum number of iterations. The clustering that minimizes
#'             inner-cluster variance is selected. Defaults to
#'             \code{defaultHiCDOCParameters$kMeansRestarts} = 20.
#'         }
#'     }
#' }
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @examples
#' data(exampleHiCDOCDataSet)
#'
#' # Retrieve parameters
#' parameters(exampleHiCDOCDataSet)
#'
#' # Set parameters
#' parameters(exampleHiCDOCDataSet) <- list("smallChromosomeThreshold" = 50)
#' parameters(exampleHiCDOCDataSet) <- list(
#'     "weakPositionThreshold" = 10,
#'     "kMeansRestarts" = 30
#' )
NULL

#### parameters ####
#' @describeIn HiCDOCDataSet-parameters
#' Retrieves the parameters used for filtering, normalization, and prediction of
#' compartments. See
#' \code{\link{filterSmallChromosomes}},
#' \code{\link{filterSparseReplicates}},
#' \code{\link{filterWeakPositions}},
#' \code{\link{normalizeDistanceEffect}}, and
#' \code{\link{detectCompartments}},
#' for details on how these parameters are used.
#' @export
setGeneric(
    name = "parameters",
    def = function(object) {
        standardGeneric("parameters")
    }
)

#### parameters <- ####
#' @describeIn HiCDOCDataSet-parameters
#' Sets the parameters used for filtering, normalization, and prediciton of
#' compartments. See
#' \code{\link{filterSmallChromosomes}},
#' \code{\link{filterSparseReplicates}},
#' \code{\link{filterWeakPositions}},
#' \code{\link{normalizeDistanceEffect}}, and
#' \code{\link{detectCompartments}},
#' for details on how these parameters are used.
#' @param value a named list containing the names and valued of the
#' parameters to change (see Details).
#' @export
setGeneric(
    name = "parameters<-",
    def = function(object, value) {
        standardGeneric("parameters<-")
    }
)
