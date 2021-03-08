.reduceHiCDOCChromosomes <- function(object, chromosomeNames, dropLevels) {
    chromosomeIds <- which(object@chromosomes %in% chromosomeNames)
    object@chromosomes <- object@chromosomes[chromosomeIds]
    if (dropLevels) {
        object@chromosomes <-
            gtools::mixedsort(unique(as.character(object@chromosomes)))
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

.reduceHiCDOCConditions <- function(object, conditionNames, dropLevels) {
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

.reduceHiCDOCReplicates <- function(object, replicateNames, dropLevels) {
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
#' object <- HiCDOCDataSetExample()
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
            "You should not reduce a HiCDOCDataSet after calling ",
            "'detectCompartments()'. All chromosomes, conditions, and ",
            "replicates have been used in the computations.",
            call. = FALSE
        )
    }
    chromosomeNames <- .validateNameOrId(object, chromosomes, "chromosomes")
    conditionNames <- .validateNameOrId(object, conditions, "conditions")
    replicateNames <- .validateNameOrId(object, replicates, "replicates")

    if (
        !is.null(chromosomeNames) &&
        !is.null(conditionNames) &&
        !is.null(replicateNames)
    ) {
        # TODO:
        # Filter the replicates of those conditions and chromosomes only
        # object <- ...
        return(object)
    }

    if (
        !is.null(chromosomeNames) &&
        !is.null(conditionNames)
    ) {
        # TODO:
        # Filter the conditions of those chromosomes only
        # object <- ...
        return(object)
    }

    if (
        !is.null(conditionNames) &&
        !is.null(replicateNames)
    ) {
        # TODO:
        # Filter the replicates of those conditions only
        # object <- ...
        return(object)
    }

    if (!is.null(chromosomeNames)) {
        object <- .reduceHiCDOCChromosomes(object, chromosomeNames, dropLevels)
        return(object)
    }

    if (!is.null(conditionNames)) {
        object <- .reduceHiCDOCConditions(object, conditionNames, dropLevels)
        return(object)
    }

    if (!is.null(replicateNames)) {
        stop("Provide 'conditions' to filter 'replicates'", call. = FALSE)
    }
}
