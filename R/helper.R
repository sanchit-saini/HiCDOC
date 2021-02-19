## - sparseInteractionsToMatrix ---------------------------------------------#
## --------------------------------------------------------------------------#
#' Build the interaction matrix for a chromosome in a condition and replicate.
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param chromosomeId A chromosome.
#' @param conditionId A condition.
#' @param replicateId A replicate.
#' @param filter  Shrink the matrix by removing weak rows/columns.
#'
#' @return A matrix.
#' @keywords internal
#' @noRd
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
        result <- result[
            -object@weakBins[[chromosomeId]],
            -object@weakBins[[chromosomeId]]
        ]
    }

    return(result)
}

## - matrixToSparseInteractions ---------------------------------------------#
## --------------------------------------------------------------------------#
#' Build the interactions tibble for a chromosome in a condition and replicate.
#'
#' @param m                         A matrix.
#' @param object                A \code{HiCDOCDataSet} object.
#' @param chromosomeId    A chromosome.
#' @param conditionId     A condition.
#' @param replicateId     A replicate.
#'
#' @return An interactions tibble.
#' @keywords internal
#' @noRd
matrixToSparseInteractions <- function(m,
    object,
    chromosomeId,
    conditionId,
    replicateId) {
    totalBins <- object@totalBins[[chromosomeId]]
    if (nrow(m) < totalBins) {
        refilled <- matrix(0, nrow = totalBins, ncol = totalBins)
        refilled[
            -object@weakBins[[chromosomeId]],
            -object@weakBins[[chromosomeId]]
        ] <- m
        m <- refilled
    }
    return(
        dplyr::tibble(
            chromosome = chromosomeId,
            bin.1 = rep(seq(totalBins), each = totalBins),
            bin.2 = rep(seq(totalBins), times = totalBins),
            condition = conditionId,
            replicate = replicateId,
            value = as.vector(t(m))
        ) %>% dplyr::filter(bin.1 <= bin.2 & value > 0)
    )
}

reduceHiCDOCChromosomes <- function(object, chromosomes, dropLevels) {
    numchr <- which(object@chromosomes %in% chromosomes)
    object@chromosomes <-
    object@chromosomes[numchr]
    if (dropLevels == TRUE) {
        object@chromosomes <-
          gtools::mixedsort(as.character(object@chromosomes))
    }
  
    object@weakBins <- object@weakBins[numchr]
    object@totalBins <- object@totalBins[numchr]
    for (slotName in c(
        "interactions",
        "distances",
        "diagonalRatios",
        "compartments",
        "concordances",
        "differences",
        "centroids",
        "positions"
  )) {
        if (!is.null(slot(object, slotName))) {
            slot(object, slotName) %<>%
                dplyr::filter(chromosome %in% chromosomes)
            if (dropLevels == TRUE) {
                slot(object, slotName) %<>%
                    dplyr::mutate(chromosome = droplevels(chromosome))
            }
        }
    }
    return(object)
}

reduceHiCDOCConditions <- function(object, conditions, dropLevels) {
    numcond <- which(object@conditions %in% conditions)
    object@conditions <-
        object@conditions[numcond]
    object@replicates <-
        object@replicates[numcond]
    for (slotName in c(
        "interactions",
        "distances",
        "diagonalRatios",
        "compartments",
        "concordances",
        "centroids"
    )) {
        if (!is.null(slot(object, slotName))) {
            slot(object, slotName) %<>%
                dplyr::filter(condition %in% conditions)
            if (dropLevels == TRUE) {
                slot(object, slotName) %<>%
                    dplyr::mutate(condition = droplevels(condition))
            }
        }
    }
    return(object)
}

reduceHiCDOCReplicates <- function(object, replicates, dropLevels) {
    numrep <- which(object@replicates %in% replicates)
    object@conditions <-
        object@conditions[numrep]
    object@replicates <-
        object@replicates[numrep]
    for (slotName in c(
        "interactions",
        "distances",
        "diagonalRatios",
        "concordances"
    )) {
        if (!is.null(slot(object, slotName))) {
            slot(object, slotName) %<>%
                dplyr::filter(replicate %in% replicates)
            if (dropLevels == TRUE) {
                slot(object, slotName) %<>%
                    dplyr::mutate(replicate = droplevels(replicate))
            }
        }
    }
    return(object)
}
#' reduce an HiCDOCDataSet
#'
#' @param object and HiCDOCDataSet object
#' @param chromosomes default to NULL, chromosomes ID to keep
#' @param conditions default to NULL, conditions ID to keep
#' @param replicates default to NULL, replicates ID to keep
#' @param dropLevels Logical, default to TRUE. Chromosomes, conditions and
#' replicates are in factor format, should the unused levels be removed ?
#' It should be set to FALSE if reduced objects are meant to be re-combined
#' later.
#'
#' @return a HiCDOCDataSET object, reduced by keeping only the chromosomes,
#' conditions or replicates given in parameters
#' @export
#'
#' @examples
#' object <- HiCDOCExample()
#' small <- reduceHiCDOCDataSet(object, chromosomes = "17", replicates = "1")
reduceHiCDOCDataSet <- function(object,
    chromosomes = NULL,
    conditions = NULL,
    replicates = NULL,
    dropLevels = TRUE) {
    if (!is.null(object@differences)) {
          warning(paste("You should not reduce an HiCDOCDataSet object after",
              "a call to detectCompartments(). All the chromosomes,",
              "conditons and replicates have been used in the",
              "computations",
              sep = " "
          ))
      }
    chromosomes <- testValidId(object, chromosomes, "chromosomes")
    conditions <- testValidId(object, conditions, "conditions")
    replicates <- testValidId(object, replicates, "replicates")

    if (!is.null(chromosomes)) {
        object <- reduceHiCDOCChromosomes(object, chromosomes, dropLevels)
    }

    if (!is.null(conditions)) {
        object <- reduceHiCDOCConditions(object, conditions, dropLevels)
    }

    if (!is.null(replicates)) {
        object <- reduceHiCDOCReplicates(object, replicates, dropLevels)
    }

    object@totalReplicates <- length(object@replicates)
    object@totalReplicatesPerCondition <-
        vapply(c(1, 2), function(x) {
            length(which(object@conditions == x))
        }, FUN.VALUE = 0)

    return(object)
}
