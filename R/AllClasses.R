#' List of HiCDOC default parameters.
#'
#' @name HiCDOCDefaultParameters
#' @docType data
#' @keywords data
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp, useParameters = HiCDOCDefaultParameters)
#'
#' @export
HiCDOCDefaultParameters <- list(kMeansIterations = 10,
                                kMeansDistance   = 0.0001)


###############################################################################
### HiCDOC S4 class definition
###############################################################################
#' Infrastructure for HiCDOC experiment and differential interaction
#'
#' \code{HiCDOC} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details \code{HiCDOCExp} does this and that...
#' TODO
#'
#' @name HiCDOCExp
#' @rdname HiCDOCExp
#' @docType class
#' @aliases HiCDOCExp HiCDOCExp-class
#'
#' @slot inputMatrix  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("HiCDOCExp", slots = c(inputMatrixPath    = "ANY",
                                interactionMatrix  = "ANY",
                                chromosomes        = "ANY",
                                replicates         = "ANY",
                                nReplicates        = "ANY",
                                conditions         = "ANY",
                                binSize            = "ANY",
                                sampleSize         = "ANY",
                                distances          = "ANY",
                                compartments       = "ANY",
                                concordances       = "ANY",
                                DIR                = "ANY",
                                loessSpan          = "ANY",
                                parameters         = "ANY")
)


##- HiCDOCExp S4 class constructor -------------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname HiCDOCExp
#' @docType class
#'
#' @param inputMatrix A matrix with the data.
#'
#' @return \code{HiCDOCExp} constructor returns an \code{HiCDOCExp}
#'         object of class S4.
#'
#' @examples
#' basedir    <- system.file("extdata", package="HiCDOC", mustWork = TRUE)
#' matrix     <- read.csv(file.path(basedir, "sampleMatrix.tsv"))
#'
#' srnaExp <- HiCDOCExp(matrix)
#' srnaExp
#'
#' @export
HiCDOCExp <- function(matrix     = NULL,
                      sampleSize = 20000,
                      loessSpan  = 0.75) {

    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#

    ##- bamFiles
    if (is.null(matrix)) {
        stop("'matrix' must be the path to a matrix", call. = FALSE)
    }

    if (!is.character(matrix)) {
        stop("'matrix' must be a character.", call. = FALSE)
    }

    if (!file.exists(matrix)) {
        stop("'matrix' must be a valid file.", call. = FALSE)
    }

    ##- end checking ---------------------------------------------------------#

    object <- new("HiCDOCExp")

    object@inputMatrixPath <- matrix
    object@sampleSize      <- sampleSize
    object@loessSpan       <- loessSpan

    return(invisible(object))
}


##- Example constructor ------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Example constructor
#'
#' This function provides an example of a \code{HiCDOCExp} object from
#' TODO
#'
#' @return An \code{HiCDOCExp} object called '\code{exp}'.
#'
#' @examples
#' ## The 'HiCDOCExp' object in this example was constructed by:
#'
#' \dontrun{
#'
#' TODO
#' }
#'
#' exp <- HiCDOCExample()
#' exp
#'
#' @export
HiCDOCExample <- function() {
    object <- NULL
    basedir <- system.file("extdata", package = "HiCDOC", mustWork = TRUE)
    matrix <- file.path(basedir, "sample.tsv")
    object <- HiCDOCExp(matrix)
    return(invisible(object))
}


##- Main method --------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Main method
#'
#' TODO
#'
#' @export
HiCDOC <- function(object) {
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
    DIR(object, pvalue = 0.05)
    concordances(object)
    compartments(object)
    return(invisible(object))
}
