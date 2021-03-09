#' @description
#' Builds the full interaction matrix of a given chromosome, condition, and
#' replicate.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' The name of a chromosome.
#' @param conditionName
#' The name of a condition.
#' @param replicateName
#' The name of a replicate.
#' @param filter
#' Whether or not to shrink the matrix by removing weak rows/columns. Defaults
#' to FALSE.
#'
#' @return
#' A matrix.
#'
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

#' @description
#' Builds the sparse interactions tibble from a matrix for a given chromosome,
#' condition, and replicate.
#'
#' @param m
#' A matrix.
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' The name of a chromosome.
#' @param conditionName
#' The name of a condition.
#' @param replicateName
#' The name of a replicate.
#'
#' @return
#' A tibble of interactions.
#'
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
            condition = conditionName,
            replicate = replicateName,
            bin.1 = rep(seq(totalBins), each = totalBins),
            bin.2 = rep(seq(totalBins), times = totalBins),
            interaction = as.vector(t(m))
        ) %>%
        dplyr::filter(bin.1 <= bin.2) %>%
        dplyr::filter(interaction > 0) %>%
        .sortInteractions(
            object@chromosomes,
            object@conditions,
            object@replicates
        )

    return(interactions)
}
