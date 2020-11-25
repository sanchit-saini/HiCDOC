#' Check object@parameters, return default if NULL
#'
#' @param parameters list of parameters
#' @param paramnames character vector, names of parameters to check
#'
#' @return list of updated parameters, default from HiCDOCDefaultParameters()
#' if null.
checkParameters <- function(parameters, paramnames) {
    nullparam <- paramnames[ vapply(parameters[paramnames], is.null, TRUE) ]
    parameters[nullparam] <- HiCDOCDefaultParameters[nullparam]
    return(parameters)
}


##- sparseInteractionsToMatrix -----------------------------------------------#
##----------------------------------------------------------------------------#
#' Build the interaction matrix for a chromosome in a condition and replicate.
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param chromosomeId A chromosome.
#' @param conditionId A condition.
#' @param replicateId A replicate.
#' @param filter  Shrink the matrix by removing weak rows/columns.
#'
#' @return A matrix.
sparseInteractionsToMatrix <- function(object,
                                       chromosomeId,
                                       conditionId,
                                       replicateId,
                                       filter = FALSE) {
    totalBins <- object@totalBins[[chromosomeId]]

    interactions <- object@interactions %>%
        dplyr::filter(chromosome == chromosomeId) %>%
        dplyr::filter(condition == conditionId) %>%
        dplyr::filter(replicate == replicateId) %>%
        dplyr::mutate(bin.1 = position.1 / object@binSize + 1,
                   bin.2 = position.2 / object@binSize + 1) %>%
        dplyr::select(bin.1, bin.2, value) %>%
        dplyr::filter(value != 0) %>%
        as.matrix()

    if (nrow(interactions) == 0) {
        return(matrix(0, nrow = 0, ncol = 0))
    }
    if (!is.numeric(interactions)) {
        stop("Error: non numeric matrix of interactions.", call. = TRUE)
    }
    result <- matrix(0, nrow = totalBins, ncol = totalBins)
    result[interactions[, c(1, 2)]] <- interactions[, 3]
    result <- result + t(result) - diag(diag(result))
    if (!isSymmetric(result)) {
        stop("Matrix is not symmetric.")
    }

    if (filter && length(object@weakBins[[chromosomeId]]) > 0) {
        result <- result[-object@weakBins[[chromosomeId]],
                         -object@weakBins[[chromosomeId]]]
    }

    return (result)
}

##- matrixToSparseInteractions -----------------------------------------------#
##----------------------------------------------------------------------------#
#' Build the interactions tibble for a chromosome in a condition and replicate.
#'
#' @param m                         A matrix.
#' @param object                A \code{HiCDOCDataSet} object.
#' @param chromosomeId    A chromosome.
#' @param conditionId     A condition.
#' @param replicateId     A replicate.
#'
#' @return An interactions tibble.
matrixToSparseInteractions <- function(m,
                                       object,
                                       chromosomeId,
                                       conditionId,
                                       replicateId) {
    totalBins <- object@totalBins[[chromosomeId]]
    if (nrow(m) < totalBins) {
        refilled <- matrix(0, nrow = totalBins, ncol = totalBins)
        refilled[-object@weakBins[[chromosomeId]],
                 -object@weakBins[[chromosomeId]]] <- m
        m <- refilled
    }
    return (
        tibble(
            chromosome = chromosomeId,
            position.1 = (rep(seq(totalBins), each = totalBins) - 1) *
                object@binSize,
            position.2 = (rep(seq(totalBins), times = totalBins) - 1) *
                object@binSize,
            condition = conditionId,
            replicate = replicateId,
            value = as.vector(t(m))
        ) %>% dplyr::filter(position.1 <= position.2 & value != 0)
    )

}





#' reduce an HiCDOCDataSet
#'
#' @param object and HiCDOCDataSet object
#' @param chromosome default to NULL, chromohomes ID to keep
#' @param condition default to NULL, conditions ID to keep
#' @param replicate default to NULL, replicate ID to keep
#'
#' @return a HiCDOCDataSET object, reduced by keeping only the chromosomes, 
#' conditions or replicates given in parameters
#' @export
#'
#' @examples
#' object <- HiCDOCExample()
#' objectReduced <-reduceHiCDOCDataSet(object, chromosome = "17", replicate = "1")
reduceHiCDOCDataSet <- function(object, 
                                chromosomes = NULL, 
                                conditions = NULL, 
                                replicates=NULL){
    if(!is.null(chromosomes)){
        numchr <- which(object@chromosomes %in% chromosomes)
        object@chromosomes <- 
            object@chromosomes[numchr]
        object@weakBins <- object@weakBins[numchr]
        object@totalBins <- object@totalBins[numchr]
        object@interactions <- 
            object@interactions[object@interactions$chromosome %in% chromosomes,]
        object@distances <- 
            object@distances[object@distances$chromosome %in% chromosomes,]
        object@diagonalRatios <- 
            object@diagonalRatios[object@diagonalRatios$chromosome %in% chromosomes,]
        object@compartments <- 
            object@compartments[object@compartments$chromosome %in% chromosomes,]
        object@concordances <- 
            object@concordances[object@concordances$chromosome %in% chromosomes,]
        object@differences <- 
            object@differences[object@differences$chromosome %in% chromosomes,]
        object@centroids <- 
            object@centroids[object@centroids$chromosome %in% chromosomes,]
    }
    
    if(!is.null(conditions)){
        numcond <- which(object@conditions %in% conditions)
        object@conditions <- 
            object@conditions[numcond]
        object@replicates <- 
            object@replicates[numcond]
        object@interactions <- 
            object@interactions[object@interactions$condition %in% conditions,]
        object@distances <- 
            object@distances[object@distances$condition %in% conditions,]
        object@diagonalRatios <- 
            object@diagonalRatios[object@diagonalRatios$condition %in% conditions,]
        object@compartments <- 
            object@compartments[object@compartments$condition %in% conditions,]
        object@concordances <- 
            object@concordances[object@concordances$condition %in% conditions,]
        object@differences <- 
            object@differences[object@differences$condition %in% conditions,]
        object@centroids <- 
            object@centroids[object@centroids$condition %in% conditions,]
    }
    
    if(!is.null(replicates)){
        numrep <- which(object@replicates %in% replicates)
        object@conditions <- 
            object@conditions[numrep]
        object@replicates <- 
            object@replicates[numrep]
        object@interactions <- 
            object@interactions[object@interactions$replicate %in% replicates,]
        object@distances <- 
            object@distances[object@distances$replicate %in% replicates,]
        object@diagonalRatios <- 
            object@diagonalRatios[object@diagonalRatios$chromosome %in% chromosomes,]
        object@compartments <- 
            object@compartments[object@compartments$chromosome %in% chromosomes,]
        object@concordances <- 
            object@concordances[object@concordances$chromosome %in% chromosomes,]
        object@differences <- 
            object@differences[object@differences$chromosome %in% chromosomes,]
        object@centroids <- 
            object@centroids[object@centroids$chromosome %in% chromosomes,]
    }
    
    object@totalReplicates <- length(object@replicates)
    object@totalReplicatesPerCondition <-
        vapply(c(1, 2), function(x) {
            length(which(object@conditions == x))
        }, FUN.VALUE = 0)
    
    return(object)
    
}
