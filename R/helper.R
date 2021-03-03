## - sparseInteractionsToMatrix ---------------------------------------------#
## --------------------------------------------------------------------------#
#' Build the interaction matrix for a chromosome in a condition and replicate.
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param chromosomeName A chromosome.
#' @param conditionName A condition.
#' @param replicateName A replicate.
#' @param filter  Shrink the matrix by removing weak rows/columns.
#'
#' @return A matrix.
#' @keywords internal
#' @noRd
sparseInteractionsToMatrix <- function(
    object,
    chromosomeName,
    conditionName,
    replicateName,
    filter = FALSE
) {
    totalBins <- object@totalBins[[chromosomeName]]

    interactions <-
        object@interactions %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::filter(condition == conditionName) %>%
        dplyr::filter(replicate == replicateName) %>%
        dplyr::select(bin.1, bin.2, interaction) %>%
        dplyr::filter(interaction != 0) %>%
        as.matrix()

    if (nrow(interactions) == 0) return(matrix(0, nrow = 0, ncol = 0))
    if (!is.numeric(interactions)) {
        stop("Non-numeric matrix of interactions.", call. = FALSE)
    }
    result <- matrix(0, nrow = totalBins, ncol = totalBins)
    result[interactions[, c(1, 2)]] <- interactions[, 3]
    result <- result + t(result) - diag(diag(result))
    if (!isSymmetric(result)) {
        stop("Matrix is not symmetric.", call. = FALSE)
    }

    if (filter && length(object@weakBins[[chromosomeName]]) > 0) {
        result <- result[
            -object@weakBins[[chromosomeName]],
            -object@weakBins[[chromosomeName]]
        ]
    }

    if (nrow(result) == 0) {
        message(
            "No interactions for chromosome ",
            chromosomeName,
            ", condition ",
            conditionName,
            ", replicate ",
            replicateName,
            "."
        )
    }

    return(result)
}

## - matrixToSparseInteractions ---------------------------------------------#
## --------------------------------------------------------------------------#
#' Build the interactions tibble for a chromosome in a condition and replicate.
#'
#' @param m                         A matrix.
#' @param object                A \code{HiCDOCDataSet} object.
#' @param chromosomeName    A chromosome.
#' @param conditionName     A condition.
#' @param replicateName     A replicate.
#'
#' @return An interactions tibble.
#' @keywords internal
#' @noRd
matrixToSparseInteractions <- function(
    m,
    object,
    chromosomeName,
    conditionName,
    replicateName
) {
    totalBins <- object@totalBins[[chromosomeName]]
    if (nrow(m) < totalBins) {
        refilled <- matrix(0, nrow = totalBins, ncol = totalBins)
        refilled[
            -object@weakBins[[chromosomeName]],
            -object@weakBins[[chromosomeName]]
        ] <- m
        m <- refilled
    }

    interactions <-
        dplyr::tibble(
            chromosome = chromosomeName,
            bin.1 = rep(seq(totalBins), each = totalBins),
            bin.2 = rep(seq(totalBins), times = totalBins),
            condition = conditionName,
            replicate = replicateName,
            interaction = as.vector(t(m))
        ) %>%
        dplyr::filter(bin.1 <= bin.2 & interaction > 0)

    return(interactions)
}

reduceHiCDOCChromosomes <- function(object, chromosomeNames, dropLevels) {
    chromosomeIds <- which(object@chromosomes %in% chromosomeNames)
    object@chromosomes <- object@chromosomes[chromosomeIds]
    if (dropLevels) {
        object@chromosomes <-
            gtools::mixedsort(as.character(object@chromosomes))
    }

    object@weakBins <- object@weakBins[chromosomeIds]
    object@totalBins <- object@totalBins[chromosomeIds]
    for (
        slotName in c(
            "interactions",
            "distances",
            "selfInteractionRatios",
            "compartments",
            "concordances",
            "differences",
            "centroids",
            "positions"
        )
    ) {
        if (!is.null(slot(object, slotName))) {
            slot(object, slotName) %<>%
                dplyr::filter(chromosome %in% chromosomeNames)
            if (dropLevels) {
                slot(object, slotName) %<>%
                    dplyr::mutate(chromosome = droplevels(chromosome))
            }
        }
    }
    return(object)
}

reduceHiCDOCConditions <- function(object, conditionNames, dropLevels) {
    conditionIds <- which(object@conditions %in% conditionNames)
    object@conditions <- object@conditions[conditionIds]
    object@replicates <- object@replicates[conditionIds]
    for (
        slotName in c(
            "interactions",
            "distances",
            "selfInteractionRatios",
            "compartments",
            "concordances",
            "centroids"
        )
    ) {
        if (!is.null(slot(object, slotName))) {
            slot(object, slotName) %<>%
                dplyr::filter(condition %in% conditionNames)
            if (dropLevels) {
                slot(object, slotName) %<>%
                    dplyr::mutate(condition = droplevels(condition))
            }
        }
    }
    return(object)
}

reduceHiCDOCReplicates <- function(object, replicateNames, dropLevels) {
    replicateIds <- which(object@replicates %in% replicateNames)
    object@conditions <- object@conditions[replicateIds]
    object@replicates <- object@replicates[replicateIds]
    for (
        slotName in c(
            "interactions",
            "distances",
            "selfInteractionRatios",
            "concordances"
        )
    ) {
        if (!is.null(slot(object, slotName))) {
            slot(object, slotName) %<>%
                dplyr::filter(replicate %in% replicateNames)
            if (dropLevels) {
                slot(object, slotName) %<>%
                    dplyr::mutate(replicate = droplevels(replicate))
            }
        }
    }
    return(object)
}

#' reduce a HiCDOCDataSet
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
reduceHiCDOCDataSet <- function(
    object,
    chromosomes = NULL,
    conditions = NULL,
    replicates = NULL,
    dropLevels = TRUE
) {
    if (!is.null(object@differences)) {
        warning(
            "You should not reduce a HiCDOCDataSet object after calling ",
            "'detectCompartments()'. All chromosomes, conditions, and ",
            "replicates have been used in the computations.",
            call. = FALSE
        )
    }
    chromosomeNames <- validateNameOrId(object, chromosomes, "chromosomes")
    conditionNames <- validateNameOrId(object, conditions, "conditions")
    replicateNames <- validateNameOrId(object, replicates, "replicates")

    if (!is.null(chromosomeNames)) {
        object <- reduceHiCDOCChromosomes(object, chromosomeNames, dropLevels)
    }

    if (!is.null(conditionNames)) {
        object <- reduceHiCDOCConditions(object, conditionNames, dropLevels)
    }

    if (!is.null(replicateNames)) {
        object <- reduceHiCDOCReplicates(object, replicateNames, dropLevels)
    }

    return(object)
}
