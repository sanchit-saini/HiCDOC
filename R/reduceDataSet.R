#' @description
#' Reduces a \code{\link{HiCDOCDataSet}} by keeping only given chromosomes.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeNames
#' The names of chromosomes to keep.
#' @param dropLevels
#' Whether or not to also remove unused factor levels after filtering. Should
#' be set to FALSE if the reduced objects are meant to be re-combined later.
#' Defaults to TRUE.
#'
#' @return
#' A reduced \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.reduceHiCDOCChromosomes <- function(object, chromosomeNames, dropLevels) {
    chromosomeIds <- which(object@chromosomes %in% chromosomeNames)
    object@chromosomes <- object@chromosomes[chromosomeIds]
    if (dropLevels) {
        object@chromosomes <-
            gtools::mixedsort(unique(as.character(object@chromosomes)))
    }

    object@weakBins <- object@weakBins[chromosomeIds]
    object@totalBins <- object@totalBins[chromosomeIds]
    object@validReplicates <- object@validReplicates[chromosomeIds]
    object@validConditions <- object@validConditions[chromosomeIds]
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

#' @description
#' Reduces a \code{\link{HiCDOCDataSet}} by keeping only given conditions.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param conditionNames
#' The names of conditions to keep.
#' @param dropLevels
#' Whether or not to also remove unused factor levels after filtering. Should
#' be set to FALSE if the reduced objects are meant to be re-combined later.
#' Defaults to TRUE.
#'
#' @return
#' A reduced \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.reduceHiCDOCConditions <- function(object, conditionNames, dropLevels) {
    conditionIds <- which(object@conditions %in% conditionNames)
    object@conditions <- object@conditions[conditionIds]
    object@replicates <- object@replicates[conditionIds]

    selection <- lapply(object@validConditions,
                        FUN = function(x)
                            which(x %in% conditionNames))
    object@validConditions <- mapply(FUN = function(x, y) x[y],
                                     object@validConditions,
                                     selection,
                                     SIMPLIFY = FALSE)
    object@validReplicates <- mapply(FUN = function(x, y) x[y],
                                     object@validReplicates,
                                     selection,
                                     SIMPLIFY = FALSE)
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

#' @description
#' Reduces a \code{\link{HiCDOCDataSet}} by keeping only given replicates.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param replicateNames
#' The names of replicates to keep.
#' @param dropLevels
#' Whether or not to also remove unused factor levels after filtering. Should
#' be set to FALSE if the reduced objects are meant to be re-combined later.
#' Defaults to TRUE.
#'
#' @return
#' A reduced \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.reduceHiCDOCReplicates <- function(object, replicateNames, dropLevels) {
    replicateIds <- which(object@replicates %in% replicateNames)
    object@conditions <- object@conditions[replicateIds]
    object@replicates <- object@replicates[replicateIds]

    selection <- lapply(object@validReplicates,
                        FUN = function(x)
                            which(x %in% replicateNames))
    object@validConditions <- mapply(FUN = function(x, y) x[y],
                                     object@validConditions,
                                     selection,
                                     SIMPLIFY = FALSE)
    object@validReplicates <- mapply(FUN = function(x, y) x[y],
                                     object@validReplicates,
                                     selection,
                                     SIMPLIFY = FALSE)
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

#' @title
#' Reduce a \code{\link{HiCDOCDataSet}}.
#'
#' @description
#' Reduces a \code{\link{HiCDOCDataSet}} by keeping only given chromosomes,
#' conditions, or replicates.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomes
#' The chromosome names or indices in \code{chromosomes(object)} to keep.
#' Defaults to NULL.
#' @param conditions
#' The condition names in \code{conditions(object)} to keep. Defaults to NULL.
#' @param replicates
#' The replicate names in \code{replicates(object)} to keep. Defaults to NULL.
#' @param dropLevels
#' Whether or not to also remove unused factor levels after filtering. Should
#' be set to FALSE if the reduced objects are meant to be re-combined later.
#' Defaults to TRUE.
#'
#' @return
#' A reduced \code{\link{HiCDOCDataSet}}.
#'
#' @examples
#' data(HiCDOCDataSetExample)
#' reduced <- reduceHiCDOCDataSet(HiCDOCDataSetExample, chromosomes = c(1, 2))
#'
#' @export
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
    chromosomeNames <- .validateNames(object, chromosomes, "chromosomes")
    conditionNames <- .validateNames(object, conditions, "conditions")
    replicateNames <- .validateNames(object, replicates, "replicates")

    if (!is.null(chromosomeNames)) {
        object <- .reduceHiCDOCChromosomes(object, chromosomeNames, dropLevels)
    }

    if (!is.null(conditionNames)) {
        object <- .reduceHiCDOCConditions(object, conditionNames, dropLevels)
    }

    if (!is.null(replicateNames)) {
        object <- .reduceHiCDOCReplicates(object, replicateNames, dropLevels)
    }

    return(object)
}
