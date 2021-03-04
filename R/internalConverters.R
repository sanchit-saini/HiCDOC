## - .sparseInteractionsToMatrix ---------------------------------------------#
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
.sparseInteractionsToMatrix <- function(
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

## - .matrixToSparseInteractions ---------------------------------------------#
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
.matrixToSparseInteractions <- function(
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
