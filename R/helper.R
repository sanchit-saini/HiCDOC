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
            bin.1 = rep(seq(totalBins), each = totalBins),
            bin.2 = rep(seq(totalBins), times = totalBins),
            condition = conditionId,
            replicate = replicateId,
            value = as.vector(t(m))
        ) %>% dplyr::filter(bin.1 <= bin.2 & value > 0)
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
        for(slotName in c("interactions", 
                          "distances", 
                          "diagonalRatios", 
                          "compartments", 
                          "concordances", 
                          "differences", 
                          "centroids", 
                          "positions")){
            if(!is.null(slot(object, slotName))){
                slot(object, slotName) %<>% filter(chromosome %in% chromosomes)
            }
        }
    }
    
    if(!is.null(conditions)){
        numcond <- which(object@conditions %in% conditions)
        object@conditions <- 
            object@conditions[numcond]
        object@replicates <- 
            object@replicates[numcond]
        for(slotName in c("interactions", 
                          "distances", 
                          "diagonalRatios", 
                          "compartments", 
                          "concordances", 
                          "centroids")){
            if(!is.null(slot(object, slotName))){
                slot(object, slotName) %<>% filter(condition %in% conditions)
            }
        }
    }
    
    if(!is.null(replicates)){
        numrep <- which(object@replicates %in% replicates)
        object@conditions <- 
            object@conditions[numrep]
        object@replicates <- 
            object@replicates[numrep]
        for(slotName in c("interactions", 
                          "distances", 
                          "diagonalRatios",
                          "concordances")){
            if(!is.null(slot(object, slotName))){
                slot(object, slotName) %<>% filter(replicate %in% replicates)
            }
        }
    }
    
    object@totalReplicates <- length(object@replicates)
    object@totalReplicatesPerCondition <-
        vapply(c(1, 2), function(x) {
            length(which(object@conditions == x))
        }, FUN.VALUE = 0)
    
    return(object)
    
}
