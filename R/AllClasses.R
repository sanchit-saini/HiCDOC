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
HiCDOCDefaultParameters <- list(kMeansIterations = 50,
                                kMeansDistance   = 0.0001,
                                kMeansRestarts   = 10,
                                sampleSize       = 20000,
                                loessSpan        = 0.75)

###############################################################################
### HiCDOCDataSet S4 class definition
###############################################################################
#' Infrastructure for HiCDOC data set
#'
#' \code{HiCDOCDataSet} is an S4 class providing the infrastructure (slots)
#' to store the input data.
#'
#' @details \code{HiCDOCDataSet} does this and that...
#' TODO
#'
#' @name HiCDOCDataSet
#' @rdname HiCDOCDataSet
#' @docType class
#' @aliases HiCDOCDatSet HiCDOCDataSet-class
#'
#' @slot inputMatrix  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("HiCDOCDataSet", slots = c(inputMatrixPath    = "ANY",
                                    interactionMatrix  = "ANY",
                                    replicates         = "ANY",
                                    conditions         = "ANY")
)


##- HiCDOCDataSet S4 class constructor from sparse matrix --------------------#
##----------------------------------------------------------------------------#
#' Reads a sparse matrix a fill a \code{HiCDOCDataSet} with its content.
#'
#' @rdname HiCDOCDataSetFromSparseMatrix
#' @docType class
#'
#' @param matrix A matrix with the data.
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
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
HiCDOCDataSetFromSparseMatrix <- function(matrix = NULL) {

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

    object <- new("HiCDOCDataSet")

    object@inputMatrixPath <- matrix
    object <- parseInteractionMatrix3Columns(object)

    return(invisible(object))
}


##- HiCDOCDataSet S4 class constructor from cool files -----------------------#
##----------------------------------------------------------------------------#
#' @rdname HiCDOCDataSetFromCool
#' @docType class
#'
#' @param matrices s
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
#'
#' @examples
#' basedir <- system.file("extdata", package="HiCDOC", mustWork = TRUE)
#' data    <- read.csv(file.path(basedir, "coolData.csv"))
#' dataSet <- HiCDOCDataSetFromCool(file.path(basedir, data$FileName),
#'                                  data$Replicate,
#'                                  data$Condition)
#' dataSet
#'
#' @export
HiCDOCDataSetFromCool <- function(coolFileNames,
                                  replicates,
                                  conditions) {

    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#

    ##- coolFileNames
    if (is.null(coolFileNames)) {
        stop("'coolFileNames' must be paths to cool files", call. = FALSE)
    }
    if (is.factor(coolFileNames)) {
        coolFileNames <- as.vector(coolFileNames)
    }
    if (!is.vector(coolFileNames)) {
        stop("'coolFileNames' must be a vector.", call. = FALSE)
    }
    if (!is.character(coolFileNames)) {
        stop("'coolFileNames' must be a vector characters.", call. = FALSE)
    }
    for (coolFileName in coolFileNames)
        if (!file.exists(coolFileName)) {
            stop(paste0("Cool file name ",
                        coolFileName,
                        " is not a valid file."), call. = FALSE)
        }

    ##- conditions
    if (is.null(conditions)) {
        stop("'conditions' should not be null.", call. = FALSE)
    }
    if (is.factor(conditions)) {
        conditions <- as.vector(conditions)
    }
    if (!is.integer(conditions)) {
        stop("'replicates' should be a vector of integers.", call. = FALSE)
    }
    if (any(sort(unique(conditions)) != c(1, 2))) {
        stop("'conditions' should be '1's or '2's.", call. = FALSE)
    }
    for (i in c(1, 2)) {
        if (sum(conditions == i) < 2) {
            stop(paste0("'conditions' should contain at least two '", i, "'s."),
                 call. = FALSE)
        }
    }
    ##- end checking ---------------------------------------------------------#

    object <- new("HiCDOCDataSet")

    object@inputMatrixPath <- coolFileNames
    object@replicates      <- replicates
    object@conditions      <- conditions

    object <- parseInteractionMatrixCool(object)

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
    object  <- NULL
    basedir <- system.file("extdata", package = "HiCDOC", mustWork = TRUE)
    matrix  <- file.path(basedir, "sample.tsv")
    dataSet <- HiCDOCDataSetFromSparseMatrix(matrix)
    object  <- HiCDOCExp(dataSet)
    return(invisible(object))
}

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
                                nReplicatesPerCond = "ANY",
                                conditions         = "ANY",
                                binSize            = "ANY",
                                sampleSize         = "ANY",
                                distances          = "ANY",
                                compartments       = "ANY",
                                concordances       = "ANY",
                                DIR                = "ANY",
                                loessSpan          = "ANY",
                                kMeansIterations   = "ANY",
                                kMeansDistance     = "ANY",
                                kMeansRestarts     = "ANY",
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
HiCDOCExp <- function(dataSet    = NULL,
                      parameters = NULL,
                      binSize    = NULL) {

    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#

    ##- dataSet
    if (is.null(dataSet)) {
        stop("'datSet' must be specified", call. = FALSE)
    }

    ##- parameters


    ##- end checking ---------------------------------------------------------#

    object <- new("HiCDOCExp")

    object@interactionMatrix <- dataSet@interactionMatrix
    object@replicates        <- dataSet@replicates
    object@conditions        <- dataSet@conditions

    object@chromosomes <- as.vector(unique(object@interactionMatrix$chromosome))
    object@nReplicates <- length(object@replicates)
    object@nReplicatesPerCond <- sapply(c(1, 2),
                                        function(x) {
                                            length(which(object@conditions ==
                                                             x))
                                        })
    if (is.null(binSize)) {
        object@binSize <- min(object@interactionMatrix$`position 1`
                              [object@interactionMatrix$`position 1` > 0])
    }

    object@kMeansIterations <- HiCDOCDefaultParameters$kMeansIterations
    object@kMeansDistance   <- HiCDOCDefaultParameters$kMeansDistance
    object@kMeansRestarts   <- HiCDOCDefaultParameters$kMeansRestarts
    object@sampleSize       <- HiCDOCDefaultParameters$sampleSize
    object@loessSpan        <- HiCDOCDefaultParameters$loessSpan

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
    DIR(object, pvalue = 0.1)
    concordances(object)
    compartments(object)
    return(invisible(object))
}
