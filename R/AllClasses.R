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
                                kMeansDelta      = 0.0001,
                                kMeansRestarts   = 20,
                                sampleSize       = 20000,
                                loessSpan        = 0.75,
                                minLength        = 100)

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
setClass("HiCDOCDataSet", slots = c(
  inputPath    = "ANY",
  interactions = "ANY",
  replicates   = "ANY",
  conditions   = "ANY",
  binSize      = "ANY"
))


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

  ##- matrix
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

  object@inputPath <- matrix
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
    stop("'coolFileNames' must be a vector of characters.", call. = FALSE)
  }
  for (coolFileName in coolFileNames) {
    if (!file.exists(coolFileName)) {
      stop(
        paste("Cool file name", coolFileName, "is not a valid file."),
        call. = FALSE
      )
    }
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

  object@inputPath  <- coolFileNames
  object@replicates <- replicates
  object@conditions <- conditions

  object <- parseInteractionMatrixCool(object)

  return(invisible(object))
}


##- HiCDOCDataSet S4 class constructor from hic files ------------------------#
##----------------------------------------------------------------------------#
#' @rdname HiCDOCDataSetFromHic
#' @docType class
#'
#' @param matrices s
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
#'
#' @examples
#' basedir <- system.file("extdata", package="HiCDOC", mustWork = TRUE)
#' data    <- read.csv(file.path(basedir, "hicData.csv"))
#' dataSet <- HiCDOCDataSetFromHic(file.path(basedir, data$FileName),
#'                                  data$Replicate,
#'                                  data$Condition,
#'                                  100000)
#' dataSet
#'
#' @export
HiCDOCDataSetFromHic <- function(hicFileNames,
                                 replicates,
                                 conditions,
                                 resolution) {

  ##- checking general input arguments -------------------------------------#
  ##------------------------------------------------------------------------#

  ##- hicFileNames
  if (is.null(hicFileNames)) {
    stop("'hicFileNames' must be paths to hic files", call. = FALSE)
  }
  if (is.factor(hicFileNames)) {
    hicFileNames <- as.vector(hicFileNames)
  }
  if (!is.vector(hicFileNames)) {
    stop("'hicFileNames' must be a vector.", call. = FALSE)
  }
  if (!is.character(hicFileNames)) {
    stop("'hicFileNames' must be a vector of characters.", call. = FALSE)
  }
  for (hicFileName in hicFileNames) {
    if (!file.exists(hicFileName)) {
      stop(
        paste("hic file name", hicFileName, "is not a valid file."),
        call. = FALSE
      )
    }
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
      stop(
        paste0("'conditions' should contain at least two '", i, "'s."),
        call. = FALSE
      )
    }
  }

  ##- resolution
  if (is.null(resolution)) {
    stop("'resolution' should not be null.", call. = FALSE)
  }
  if (is.factor(resolution)) {
    resolution <- as.vector(resolution)
  }
  if (!is.numeric(resolution)) {
    stop("'resolution' should be an integer.", call. = FALSE)
  }
  if (length(resolution) != 1) {
    stop("'resolution' should not be a vector.", call. = FALSE)
  }

  ##- end checking ---------------------------------------------------------#

  object <- new("HiCDOCDataSet")

  object@inputPath  <- hicFileNames
  object@replicates <- replicates
  object@conditions <- conditions
  object@binSize    <- resolution

  object <- parseInteractionMatrixHic(object)

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
setClass("HiCDOCExp", slots = c(
  inputPath                   = "ANY",
  interactions                = "ANY",
  chromosomes                 = "ANY",
  replicates                  = "ANY",
  totalReplicates             = "ANY",
  totalReplicatesPerCondition = "ANY",
  conditions                  = "ANY",
  binSize                     = "ANY",
  sampleSize                  = "ANY",
  distances                   = "ANY",
  compartments                = "ANY",
  concordances                = "ANY",
  differences                 = "ANY",
  centroids                   = "ANY",
  loessSpan                   = "ANY",
  kMeansIterations            = "ANY",
  kMeansDelta                 = "ANY",
  kMeansRestarts              = "ANY",
  minLength                   = "ANY",
  parameters                  = "ANY"
))


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

  object@interactions <- dataSet@interactions
  object@replicates   <- dataSet@replicates
  object@conditions   <- dataSet@conditions

  object@chromosomes     <- as.vector(unique(object@interactions$chromosome))
  object@totalReplicates <- sum(sapply(object@replicates, length))

  object@totalReplicatesPerCondition <- sapply(c(1, 2), function(x) {
    length(which(object@conditions == x))
  })

  if (is.null(binSize)) {
    object@binSize <- min(
      object@interactions$`position 1`[object@interactions$`position 1` > 0]
    )
  }

  object@kMeansIterations <- HiCDOCDefaultParameters$kMeansIterations
  object@kMeansDelta      <- HiCDOCDefaultParameters$kMeansDelta
  object@kMeansRestarts   <- HiCDOCDefaultParameters$kMeansRestarts
  object@sampleSize       <- HiCDOCDefaultParameters$sampleSize
  object@loessSpan        <- HiCDOCDefaultParameters$loessSpan
  object@minLength        <- HiCDOCDefaultParameters$minLength

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
  differences(object, pvalue = 0.1)
  concordances(object)
  compartments(object)
  return(invisible(object))
}
