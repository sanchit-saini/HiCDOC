checkParameters <- function(value) {
    # TODO
}

##- sparseInteractionsToFullInteractions -------------------------------------#
##----------------------------------------------------------------------------#
#' Build the full interactions tibble filled with zeros for a chromosome in a
#' condition and replicate.
#'
#' @param object A \code{HiCDOCExp} object.
#' @param chromosomeId A chromosome.
#' @param conditionId A condition.
#' @param replicateId A replicate.
#' @param filter Whether to remove weak positions. Defaults to TRUE.
#'
#' @return An interactions tibble.
sparseInteractionsToFullInteractions <- function(object,
                                                 chromosomeId,
                                                 conditionId,
                                                 replicateId,
                                                 filter = TRUE) {
    totalBins <- object@totalBins[[chromosomeId]]
    positions <- (seq_len(totalBins) - 1) * object@binSize

    interactions = object@interactions %>%
        dplyr::filter(chromosome == chromosomeId) %>%
        dplyr::filter(condition == conditionId) %>%
        dplyr::filter(replicate == replicateId)

    mirrorInteractions = interactions %>%
        dplyr::filter(position.1 != position.2) %>%
        dplyr::rename(position.1 = position.2, position.2 = position.1)

    fullInteractions <- tibble(
        chromosome = chromosomeId,
        condition = conditionId,
        replicate = replicateId,
        position.1 = rep(positions, totalBins),
        position.2 = rep(positions, each = totalBins),
        value = 0
    ) %>%
        dplyr::left_join(
            dplyr::bind_rows(interactions, mirrorInteractions),
            by = c(
                "chromosome",
                "condition",
                "replicate",
                "position.1",
                "position.2"
            )
        ) %>%
        dplyr::mutate(value = dplyr::coalesce(value.y, value.x)) %>%
        dplyr::select(-c(value.x, value.y))

    if (filter) {
        weakBins <- object@weakBins[[chromosomeId]]
        weakPositions <- (weakBins - 1) * object@binSize

        fullInteractions %<>%
            dplyr::filter(!(position.1 %in% weakPositions)) %>%
            dplyr::filter(!(position.2 %in% weakPositions))
    }

    return (fullInteractions)
}

##- sparseInteractionsToMatrix -----------------------------------------------#
##----------------------------------------------------------------------------#
#' Build the interaction matrix for a chromosome in a condition and replicate.
#'
#' @param object A \code{HiCDOCExp} object.
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
#' @param object                A \code{HiCDOCExp} object.
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
