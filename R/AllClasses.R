#### - HiCDOCDefaultParameters -------------------------------------------####
# @describeIn parameters Default parameters for the HiCDOC pipeline
#' @rdname parameters
#' @examples
#' # All the default parameters
#' HiCDOCDefaultParameters
#' 
#' @export
HiCDOCDefaultParameters <- list(
    minLengthChr = 100,
    weakPosThreshold = 0,
    sparseThreshold = 0.95,
    sampleSize = 20000,
    kMeansIterations = 50,
    kMeansDelta = 0.0001,
    kMeansRestarts = 20
)

#### - HiCDOCDataSet parsers ---------------------------------------------####

## - makeHiCDOCDataSet ------------------------------------------------------#
## --------------------------------------------------------------------------#
#' Constructor function for the \code{HiCDOCDataSet} class.
#'
#' @param inputPath    One or several file(s) that contain the input 
#' matrix/ces.
#' @param interactions The interaction matrices.
#' @param replicates   The list of the replicate names.
#' @param conditions   The names of the two conditions.
#' @param binSize      The resolution.
#' @param hicPro       Logical. Default to FALSE. Does the data sent are in
#' HiC-Pro format ?
#' @return A \code{HiCDOCDataSet} object.
#' @keywords internal
#' @noRd
makeHiCDOCDataSet <- function(
    inputPath = NULL,
    interactions = NULL,
    replicates = NULL,
    conditions = NULL,
    binSize = NULL,
    hicPro = FALSE) {
    ## - checking general input arguments -----------------------------------#
    ## ----------------------------------------------------------------------#

    ## - matrix
    if (is.null(inputPath)) {
        stop("'inputPath' must be the path(es) to a matrix(ces)",
            call. = FALSE
        )
    }
    if (!hicPro) {
        if (!is.character(inputPath)) {
            stop("'inputPath' must be a character.", call. = FALSE)
        }

        for (fileName in inputPath) {
            if (!file.exists(fileName)) {
                stop(paste("File name", fileName, "is not a valid file."),
                    call. = FALSE
                )
            }
        }
    } else {
        if (!is.list(inputPath)) {
            stop("'inputPath' must be a list", call. = FALSE)
        }

        for (fileName in unlist(inputPath)) {
            if (!file.exists(fileName)) {
                stop(paste("File name", fileName, "is not a valid file."),
                    call. = FALSE
                )
            }
        }
    }
    ## - end checking -------------------------------------------------------#

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


## - HiCDOCDataSet S4 class constructor from sparse matrix ------------------#
## --------------------------------------------------------------------------#
#' Reads a sparse matrix and fill a \code{HiCDOCDataSet} with its content.
#'
#' @rdname HiCDOCDataSetFromSparseMatrix
#' @docType class
#' @param matrix A matrix with the data.
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
#' @examples
#' link <- system.file("extdata", "sampleMatrix.tsv", package = "HiCDOC")
#' dataSet <- HiCDOCDataSetFromSparseMatrix(link)
#' dataSet
#' @export
HiCDOCDataSetFromSparseMatrix <- function(matrix = NULL) {
    ## - checking general input arguments -----------------------------------#
    ## ----------------------------------------------------------------------#

    ## - end checking -------------------------------------------------------#

    data <- makeHiCDOCDataSet(inputPath = matrix)
    object <- parseInteractionMatrix3Columns(data)
    object <- HiCDOCDataSet(object)
    return(invisible(object))
}


## - HiCDOCDataSet S4 class constructor from cool files ---------------------#
## --------------------------------------------------------------------------#
#' Construct a \code{HiCDOCDataSet} from a cool file.
#' @rdname HiCDOCDataSetFromCool
#' @docType class
#'
#' @param coolFileNames A vector of file names, each one being a .cool file.
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
    coolFileNames,
    replicates,
    conditions) {

    ## - checking general input arguments -----------------------------------#
    ## ----------------------------------------------------------------------#

    ## - coolFileNames
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
        coolFilePath <- strsplit(coolFileName, "::/")[[1]][1]
        if (!file.exists(coolFilePath)) {
            stop(
                paste("Cool file name", coolFilePath, "is not a valid file."),
                call. = FALSE
            )
        }
    }

    ## - conditions
    if (is.null(conditions)) {
        stop("'conditions' should not be null.", call. = FALSE)
    }
    if (is.factor(conditions)) {
        conditions <- as.vector(conditions)
    }
    ## - end checking -------------------------------------------------------#

    data <- makeHiCDOCDataSet(
        inputPath = coolFileNames,
        replicates = replicates,
        conditions = conditions
    )
    object <- parseInteractionMatrixCool(data)
    object <- HiCDOCDataSet(object)

    return(invisible(object))
}


## - HiCDOCDataSet S4 class constructor from hic files ----------------------#
## --------------------------------------------------------------------------#
#' Construct a \code{HiCDOCDataSet} from a hic file.
#' @rdname HiCDOCDataSetFromHic
#' @docType class
#'
#' @param hicFileNames A vector of file names, each one being a .hic file.
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
HiCDOCDataSetFromHic <- function(hicFileNames,
    replicates,
    conditions,
    resolution) {
    ## - checking general input arguments -----------------------------------#
    ## ----------------------------------------------------------------------#

    ## - hicFileNames
    if (is.null(hicFileNames)) {
        stop("'hicFileNames' must be paths to hic files", call. = FALSE)
    }
    if (is.factor(hicFileNames)) {
        hicFileNames <- as.vector(hicFileNames)
    }
    if (!is.vector(hicFileNames) || !is.character(hicFileNames)) {
        stop("'hicFileNames' must be a character vector.", call. = FALSE)
    }

    ## - conditions
    if (is.null(conditions)) {
        stop("'conditions' should not be null.", call. = FALSE)
    }
    if (is.factor(conditions)) {
        conditions <- as.vector(conditions)
    }

    ## - resolution
    if (is.null(resolution)) {
        stop("'resolution' should not be null.", call. = FALSE)
    }
    if (!is.numeric(resolution) || length(resolution) != 1) {
        stop("'resolution' should be an integer.", call. = FALSE)
    }

    ## - end checking -------------------------------------------------------#

    data <- makeHiCDOCDataSet(
        inputPath = hicFileNames,
        replicates = replicates,
        conditions = conditions,
        binSize = resolution
    )
    object <- parseInteractionMatrixHic(data)
    object <- HiCDOCDataSet(object)

    return(invisible(object))
}

## - HiCDOCDataSet S4 class constructor from HiC-Pro files ------------------#
## --------------------------------------------------------------------------#
#' Construct a \code{HiCDOCDataSet} from hHiC-Pro files
#' @rdname HiCDOCDataSetFromHicPro
#' @docType class
#'
#' @param matrixFileNames A vector of file names, each one being a .matrix file.
#' @param bedFileNames A vector of file names, each one being a .bed file.
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
HiCDOCDataSetFromHicPro <- function(matrixFileNames,
    bedFileNames,
    replicates,
    conditions) {
    ## - checking general input arguments -----------------------------------#
    ## ----------------------------------------------------------------------#

    ## - Matrix files
    if (!is.vector(matrixFileNames) || !is.character(matrixFileNames)) {
        stop(paste("'matrixFileNames' must be a character vector,",
             "giving the links to matrix files."), call. = FALSE)
    }

    # bedFiles
    if (!is.vector(bedFileNames) || !is.character(bedFileNames)) {
        stop("'bedFileNames' must be a character vector, giving the links to
           .bed files.", call. = FALSE)
    }

    if (length(matrixFileNames) != length(bedFileNames)) {
          stop(paste("'matrixFileNames' and 'bedFileNames'",
                     "should have the same length."),
              call. = FALSE
          )
      }
    
    ## - conditions
    if (is.null(conditions)) {
        stop("'conditions' should not be null.", call. = FALSE)
    }
    if (is.factor(conditions)) {
        conditions <- as.vector(conditions)
    }
    
    ## - end checking -------------------------------------------------------#
    hicProFiles <- split(
        cbind(matrixFileNames, bedFileNames),
        seq(length(matrixFileNames))
    )
    data <- makeHiCDOCDataSet(
        inputPath = hicProFiles,
        replicates = replicates,
        conditions = conditions,
        hicPro = TRUE
    )
    temp <- parseInteractionMatrixHicPro(data)
    object <- temp[["matrix"]]
    object@binSize <- temp[["resolution"]]
    object@positions <- temp[["positions"]]
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
    basedir <-
        system.file("extdata", package = "HiCDOC", mustWork = TRUE)
    matrix <- file.path(basedir, "sampleMatrix.tsv")
    object <- HiCDOCDataSetFromSparseMatrix(matrix)
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
#' * Matrices: see \code{\link{HiCDOCDataSetFromSparseMatrix}}
#' * Cool or mCool data: see \code{\link{HiCDOCDataSetFromCool}}
#' * Hic files: see \code{\link{HiCDOCDataSetFromHic}}
#' * Hic Pro matrices and bed files: see \code{\link{HiCDOCDataSetFromHicPro}}
#' @slot inputPath The names of the matrix input files.
#' @slot interactions The interaction matrices.
#' @slot weakBins The empty bins.
#' @slot chromosomes The list of chromosomes.
#' @slot replicates The names of the replicates, with repetition along the 
#' conditions.
#' @slot totalReplicates The names of the replicates, glued with the
#' name of the conditions.
#' @slot totalReplicatesPerCondition A 2-element list, one for each condition,
#' where the union is the names of the replicates.
#' @slot conditions The names of the conditions (exactly two different), with 
#' repetition along the replicates.
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
#' @slot positions The position of the bins.
#' @seealso \code{\link{HiCDOCExample}}, 
#' \code{\link{parameters}}, 
#' \code{\link{HiCDOCDataSetFromSparseMatrix}}, 
#' \code{\link{HiCDOCDataSetFromCool}},
#' \code{\link{HiCDOCDataSetFromHic}},
#' \code{\link{HiCDOCDataSetFromHicPro}}
#' @export
#' @md
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
        parameters                  = "ANY",
        positions                   = "ANY"
    )
)


## - HiCDOCDataSet S4 class constructor -------------------------------------#
## --------------------------------------------------------------------------#
#' Constructor for the HiCDOCDataSet class.
#' 
#' This function should not be called directly. It is called by the functions
#' \code{\link{HiCDOCDataSetFromSparseMatrix}}, 
#' \code{\link{HiCDOCDataSetFromCool}}, 
#' \code{\link{HiCDOCDataSetFromHic}} and 
#' \code{\link{HiCDOCDataSetFromHicPro}}.
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
    if (is.null(object)) {
        # stop("'dataSet' must be specified", call. = FALSE)
        object <- new("HiCDOCDataSet")
    }

    # Fill binSize slot
    if (is.null(object@binSize)) {
        object@binSize <- computeBinSize(object)
    }

    if (!is.null(object@interactions)) {
        chromosomes <- gtools::mixedsort(
            unique(as.character(object@interactions$chromosome))
        )
        object@chromosomes <- chromosomes

        # Determine positions
        if (is.null(object@positions)) {
            object@positions <- determinePositions(object)
        }

        # Replace positions by bins in interactions
        object@interactions <- replacePositionsByBins(object)

        object@weakBins <- vector("list", length(chromosomes))
        names(object@weakBins) <- chromosomes

        # Put matrix in upper only, value as numeric and aggregate duplicates
        object@interactions <- reformatInteractions(object@interactions)
    }

    # Fill totalReplicates slot
    if (!is.null(object@replicates)) {
        object@totalReplicates <- length(object@replicates)
    }
    # Fill totalReplicatesPerCondition slot
    if (!is.null(object@conditions)) {
        object@totalReplicatesPerCondition <-
            vapply(unique(object@conditions), function(x) {
                length(which(object@conditions == x))
            }, FUN.VALUE = 0)
    }
    # Fill totalBins slot
    if (!is.null(object@binSize) & !is.null(object@interactions)) {
        object@totalBins <- vapply(as.character(object@chromosomes),
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
#'         \item{matrices: }{see \code{\link{HiCDOCDataSetFromSparseMatrix}}},
#'         \item{cool or mcool data: }{see \code{\link{HiCDOCDataSetFromCool}}},
#'         \item{hic data: }{see \code{\link{HiCDOCDataSetFromHic}}} and 
#'         \item{HicPro data: }{see \code{\link{HiCDOCDataSetFromHicPro}}}
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
