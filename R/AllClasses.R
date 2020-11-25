#' List of HiCDOC default parameters.
#'
#' @name HiCDOCDefaultParameters
#' @docType data
#' @keywords data
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp)
#'
#' @export
HiCDOCDefaultParameters <- list(
    kMeansIterations = 50,
    kMeansDelta      = 0.0001,
    kMeansRestarts   = 20,
    sampleSize       = 20000,
    loessSpan        = 0.75,
    minLengthChr     = 100,
    weakPosThreshold = 0
)

##- makeHiCDOCDataSet --------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Constructor function for the \code{HiCDOCDataSet} class.
#'
#' @param inputPath    One or several file(s) that contain the input matrix/ces.
#' @param interactions The interaction matrices.
#' @param replicates   The list of the replicate names.
#' @param conditions   The names of the two conditions.
#' @param binSize      The resolution.
#' @return A \code{HiCDOCDataSet} object.
#' @examples
#' basedir <- system.file("extdata", package="HiCDOC", mustWork = TRUE)
#' matrix  <- file.path(basedir, "sampleMatrix.tsv")
#' data    <- makeHiCDOCDataSet(inputPath = matrix)
makeHiCDOCDataSet <- function(
    inputPath    = NULL,
    interactions = NULL,
    replicates   = NULL,
    conditions   = NULL,
    binSize      = NULL) {
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#

    ##- matrix
    if (is.null(inputPath)) {
        stop("'inputPath' must be the path(es) to a matrix(ces)",
                 call. = FALSE)
    }

    if (!is.character(inputPath)) {
        stop("'inputPath' must be a character.", call. = FALSE)
    }

    for (fileName in inputPath) {
        if (!file.exists(fileName)) {
            stop(paste("File name", fileName, "is not a valid file."),
                    call. = FALSE)
        }
    }

    ##- end checking ---------------------------------------------------------#

    object <- new("HiCDOCDataSet")
    object@inputPath <- inputPath
    if (!is.null(interactions)) {
        object@interactions <- interactions
    }
    if (!is.null(replicates)) {
        object@replicates <- replicates
    }
    if (!is.null(conditions)) {
        object@conditions <- conditions
    }
    if (!is.null(binSize)) {
        object@binSize <- binSize
    }
    return(object)
}


##- HiCDOCDataSet S4 class constructor from sparse matrix --------------------#
##----------------------------------------------------------------------------#
#' Reads a sparse matrix and fill a \code{HiCDOCDataSet} with its content.
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
#' linkToMatrix <- system.file("extdata", "sampleMatrix.tsv",  package="HiCDOC")
#' dataSet <- HiCDOCDataSetFromSparseMatrix(linkToMatrix)
#' dataSet
#'
#' @export
HiCDOCDataSetFromSparseMatrix <- function(matrix = NULL) {
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#

    ##- end checking ---------------------------------------------------------#

    data <- makeHiCDOCDataSet(inputPath = matrix)
    object <- parseInteractionMatrix3Columns(data)
    object <- HiCDOCDataSet(object)
    return(invisible(object))
}


##- HiCDOCDataSet S4 class constructor from cool files -----------------------#
##----------------------------------------------------------------------------#
#' Construct a \code{HiCDOCDataSet} from a cool file.
#' @rdname HiCDOCDataSetFromCool
#' @docType class
#'
#' @param coolFileNames A vector of file names, each one being a .cool file.
#' @param replicates   The list of the replicate names.
#' @param conditions   The names of the two conditions.
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
#'
#' @examples
#' basedir <- system.file("extdata", package="HiCDOC")
#' data    <- read.csv(file.path(basedir, "coolData.csv"))
#' data
#' dataSet <- HiCDOCDataSetFromCool(file.path(basedir, data$FileName),
#'                                  data$Replicate,
#'                                  data$Condition)
#'
#' @export
HiCDOCDataSetFromCool <- function(
  coolFileNames,
  replicates,
  conditions
) {

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
    # Remove trailing URI in case of an mcool file
    coolFilePath <- strsplit(coolFileName, '::/')[[1]][1]
    if (!file.exists(coolFilePath)) {
      stop(
        paste("Cool file name", coolFilePath, "is not a valid file."),
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
  ##- end checking ---------------------------------------------------------#

  data   <- makeHiCDOCDataSet(inputPath  = coolFileNames,
                              replicates = replicates,
                              conditions = conditions)
  object <- parseInteractionMatrixCool(data)
  object <- HiCDOCDataSet(object)

  return(invisible(object))
}


##- HiCDOCDataSet S4 class constructor from hic files ------------------------#
##----------------------------------------------------------------------------#
#' Construct a \code{HiCDOCDataSet} from a hic file.
#' @rdname HiCDOCDataSetFromHic
#' @docType class
#'
#' @param hicFileNames A vector of file names, each one being a .hic file.
#' @param replicates   The list of the replicate names.
#' @param conditions   The names of the two conditions.
#' @param resolution   The resolution.
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
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
            stop(paste("hic file name", hicFileName, "is not a valid file."),
                 call. = FALSE)
        }
    }

    ##- conditions
    if (is.null(conditions)) {
        stop("'conditions' should not be null.", call. = FALSE)
    }
    if (is.factor(conditions)) {
        conditions <- as.vector(conditions)
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

    data   <- makeHiCDOCDataSet(
        inputPath  = hicFileNames,
        replicates = replicates,
        conditions = conditions,
        binSize    = resolution
    )
    object <- parseInteractionMatrixHic(data)
    object <- HiCDOCDataSet(object)

    return(invisible(object))
}

##- Example constructor ------------------------------------------------------#
##----------------------------------------------------------------------------#
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
#'
#' @export
HiCDOCExample <- function() {
    object  <- NULL
    basedir <-
        system.file("extdata", package = "HiCDOC", mustWork = TRUE)
    matrix  <- file.path(basedir, "sampleMatrix.tsv")
    dataSet <- HiCDOCDataSetFromSparseMatrix(matrix)
    object  <- HiCDOCDataSet(dataSet)
    return(invisible(object))
}


## HiCDOCDataSet S4 class definition ----------------------------------------#
##---------------------------------------------------------------------------#
#' Infrastructure for HiCDOC experiment and differential interaction
#'
#' \code{HiCDOCDataSet} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details HiCDOCDataSet does this
#'
#' @name HiCDOCDataSet
#' @rdname HiCDOCDataSet
#' @docType class
#' @aliases HiCDOCDataSet HiCDOCDataSet-class
#'
#' @slot inputPath The names of the matrix input files.
#' @slot interactions The interaction matrices.
#' @slot weakBins The empty bins.
#' @slot chromosomes The list of chromosomes.
#' @slot replicates The names of the replicates.
#' @slot totalReplicates The names of the replicates, glued with the
#' name of the conditions.
#' @slot totalReplicatesPerCondition A 2-element list, one for each condition,
#' where the union is the names of the replicates.
#' @slot conditions The names of the conditions (exactly two different).
#' @slot totalBins The number of bins per chromosome.
#' @slot binSize The resolution.
#' @slot distances The distribution of distances to the centroids.
#' @slot diagonalRatios TODO
#' @slot compartments The A/B compartments distribution, along the chromosomes.
#' @slot concordances The concordance distribution, along the chromosomes.
#' @slot differences The distribution of the difference of the concordance.
#' @slot centroids The position of the centroids.
#' @slot parameters An named \code{list}. The parameters for the
#' segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass(
    "HiCDOCDataSet",
    slots = c(
        inputPath                   = "ANY",
        interactions                = "ANY",
        weakBins                    = "ANY",
        chromosomes                 = "ANY",
        replicates                  = "ANY",
        totalReplicates             = "ANY",
        totalReplicatesPerCondition = "ANY",
        conditions                  = "ANY",
        totalBins                   = "ANY",
        binSize                     = "ANY",
        distances                   = "ANY",
        diagonalRatios              = "ANY",
        compartments                = "ANY",
        concordances                = "ANY",
        differences                 = "ANY",
        centroids                   = "ANY",
        parameters                  = "ANY"
    )
)


#' Number of bins of a chromosome, from the interactions matrix
#'
#' @param chromosomeId Name of the chromosome
#' @param interactions Interactions matrix
#' @param binSize Numeric value, resolution
#'
#' @return a numeric value
nbBinsChromosome <- function(chromosomeId, interactions, binSize){
  chromosomeInteractions <- interactions %>%
    dplyr::filter(chromosome == chromosomeId)
  nbBins = max(chromosomeInteractions$position.1,
               chromosomeInteractions$position.2) / binSize + 1
  return(nbBins)
}



##- HiCDOCDataSet S4 class constructor -------------------------------------#
##--------------------------------------------------------------------------#
#' @rdname HiCDOCDataSet
#' @docType class
#'
#' @param object     A \code{HiCDOCDataSet} object.
#' @param parameters A named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#' @param binSize    The resolution.
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
#'
#' @examples
#' linkToMatrix <- system.file("extdata", "sampleMatrix.tsv", package="HiCDOC")
#' srnaDataSet <- HiCDOCDataSetFromSparseMatrix(linkToMatrix)
#' srnaExp <- HiCDOCDataSet(srnaDataSet)
#' @export
HiCDOCDataSet <- function(object = NULL,
                      parameters = NULL,
                      binSize = NULL) {

      if ( is.null(object) ) {
        # stop("'dataSet' must be specified", call. = FALSE)
        object <- new("HiCDOCDataSet")
    }

    if( !is.null(object@interactions) ) {
        object@chromosomes <-
            mixedsort(as.vector(unique(object@interactions$chromosome)))

        object@weakBins <- vector("list", length(object@chromosomes))
        names(object@weakBins) <- object@chromosomes
    }
    if( !is.null(object@replicates) ) {
        object@totalReplicates <- length(object@replicates)
    }
    if( !is.null(object@conditions) ) {
        object@totalReplicatesPerCondition <-
            vapply(unique(object@conditions), function(x) {
                length(which(object@conditions == x))
            }, FUN.VALUE = 0)
    }

    if ( !is.null(binSize) ) {
        object@binSize <- binSize
    } else {
      if( !is.null(object@interactions) ) {
          object@binSize <- min(abs(
            object@interactions$position.1[object@interactions$position.1
                                           != object@interactions$position.2]
            - object@interactions$position.2[object@interactions$position.1
                                             != object@interactions$position.2]
          ))
      }
    }

    if( !is.null(object@binSize) & !is.null(object@interactions) ){
        object@totalBins <- vapply(object@chromosomes,
            function(x) nbBinsChromosome(x,
                                         object@interactions,
                                         object@binSize),
            FUN.VALUE = 0
            )
    }

    acceptedParam <- names(HiCDOCDefaultParameters)
    object@parameters <- vector("list", length(acceptedParam))
    names(object@parameters) <- acceptedParam

    if( !is.null(parameters) ){
      paramnames <- names(parameters)
      paramOK <- match(paramnames, acceptedParam)
      if(length(paramOK[is.na(paramOK)])>0) {
          warning(paste0("Unknown parameters '",
                        paste(paramnames[is.na(paramOK)], collapse="', '"),
                        "'.\nAccepted parameters: '",
                        paste(acceptedParam, collapse = "', '"), "'"),
                  call. = FALSE)
      }
      paramnames <- paramnames[ !is.na(paramOK) ]
      object@parameters[paramnames] <- parameters[paramnames]
    }

    return(invisible(object))
}


##- Main method --------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Main method.  Start the pipeline with default parameters.
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @return Returns an \code{HiCDOCDataSet} object.
#' @export
HiCDOC <- function(object) {
    object <- filterSmallChromosomes(object)
    object <- filterWeakPositions(object)
    object <- normalizeTechnicalBiases(object)
    object <- normalizeBiologicalBiases(object)
    object <- normalizeDistanceEffect(object)
    object <- detectCompartments(object)
    return(invisible(object))
}
