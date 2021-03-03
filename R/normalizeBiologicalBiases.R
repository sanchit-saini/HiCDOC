# Based on
# https://github.com/dozmorovlab/HiCcompare/blob/master/R/KRnormalization.R
normalizeKnightRuiz <- function(
    A,
    tol = 1e-6,
    minDelta = 0.1,
    maxDelta = 3
) {
    n <- nrow(A)
    e <- matrix(1, nrow = n, ncol = 1)
    x0 <- e
    g <- 0.9
    etamax <- 0.1
    eta <- etamax
    stop_tol <- tol * .5
    x <- x0
    rt <- tol^2
    v <- x * (A %*% x)
    rk <- 1 - v
    rho_km1 <- drop(t(rk) %*% rk)
    rout <- rho_km1
    rold <- rout
    while (rout > rt) {
        k <- 0
        y <- e
        innertol <- max(c(eta^2 * rout, rt))
        while (rho_km1 > innertol) {
            k <- k + 1
            if (k == 1) {
                Z <- rk / v
                p <- Z
                rho_km1 <- drop(t(rk) %*% Z)
            } else {
                beta <- rho_km1 / rho_km2
                p <- Z + beta * p
            }
            w <- x * (A %*% (x * p)) + v * p
            if (max(w) == Inf) {
                warning("KR algorithm diverges.", call. = FALSE)
                return(t(t(x[, 1] * A) * x[, 1]))
            }
            alpha <- rho_km1 / drop(t(p) %*% w)
            ap <- alpha * p
            ynew <- y + ap
            if (min(ynew) <= minDelta) {
                if (minDelta == 0) break()
                ind <- which(ap < 0)
                gamma <- min((minDelta - y[ind]) / ap[ind])
                y <- y + gamma * ap
                break()
            }
            if (max(ynew) >= maxDelta) {
                ind <- which(ynew > maxDelta)
                gamma <- min((maxDelta - y[ind]) / ap[ind])
                y <- y + gamma * ap
                break()
            }
            y <- ynew
            rk <- rk - alpha * w
            rho_km2 <- rho_km1
            Z <- rk / v
            rho_km1 <- drop(t(rk) %*% Z)
        }
        x <- x * y
        v <- x * (A %*% x)
        rk <- 1 - v
        rho_km1 <- drop(t(rk) %*% rk)
        rout <- rho_km1
        rat <- rout / rold
        rold <- rout
        res_norm <- sqrt(rout)
        eta_o <- eta
        eta <- g * rat
        if (g * eta_o^2 > 0.1) eta <- max(c(eta, g * eta_o^2))
        eta <- max(c(min(c(eta, etamax)), stop_tol / res_norm))
    }
    result <- t(t(x[, 1] * A) * x[, 1])
    return(result)
}


## - normalizeBiologicalBiasesOfChromosome -----------------------------------#
## ---------------------------------------------------------------------------#
#' Remove biological biases by normalizing with Knight-Ruiz matrix balancing.
#' @param object A \code{HiCDOCDataSet} object.
#' @param chromosomeName The name or number of the chromosome to plot.
#' If number, will be taken in \code{object@chromosomes[chromosomeName]}
#'
#' @return A \code{HiCDOCDataSet} object, with the normalized matrices.
normalizeBiologicalBiasesOfChromosome <- function(object, chromosomeName) {
    validateSlots(
        object,
        slots = c(
            "interactions",
            "totalBins",
            "conditions",
            "replicates",
            "totalBins"
        )
    )

    message("Chromosome ", chromosomeName, ": normalizing biological biases.")

    if (object@totalBins[[chromosomeName]] == -Inf) return(NULL)

    matrices <-
        mapply(
            function(conditionName, replicateName) {
                sparseInteractionsToMatrix(
                    object,
                    chromosomeName,
                    conditionName,
                    replicateName,
                    filter = TRUE
                )
            },
            object@conditions,
            object@replicates,
            SIMPLIFY = FALSE
        )

    badMatrices <-
        vapply(
            matrices,
            function(matrix) min(rowSums(matrix)) == 0,
            FUN.VALUE = TRUE
        )

    if (any(badMatrices)) {
        stop(
            "Cannot normalize matri",
            if (sum(badMatrices) != 1) "ces" else "x",
            " with empty positions:\n",
            sapply(
                which(badMatrices),
                function(x) {
                    paste0(
                        "chromosome ",
                        chromosomeName,
                        ", condition ",
                        object@conditions[x],
                        ", replicate ",
                        object@replicates[x],
                        ".\n"
                    )
                }
            ),
            "Call 'filterWeakPositions()' before normalizing."
        )
    }

    normalizedMatrices <- lapply(matrices, normalizeKnightRuiz)
    chromosomeInteractions <-
        purrr::pmap_dfr(
            list(
                normalizedMatrices,
                object@conditions,
                object@replicates
            ),
            .f = function(matrix, conditionName, replicateName) {
                matrixToSparseInteractions(
                    matrix,
                    object,
                    chromosomeName,
                    conditionName,
                    replicateName
                )
            }
        )
    return(chromosomeInteractions)
}

## - normalizeBiologicalBiases ------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Remove biological biases on a HiCDOCDataSet object by normalizing
#' with Knight-Ruiz matrix balancing.
#'
#' @rdname normalizeBiologicalBiases
#'
#' @param object A \code{HiCDOCDataSet} object.
#'
#' @return A \code{HiCDOCDataSet} object, with the normalized matrices.
#'
#' @examples
#' object <- HiCDOCExample()
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' @seealso \code{\link[HiCDOC]{normalizeTechnicalBiases}},
#' \code{\link[HiCDOC]{normalizeDistanceEffect}} and
#' \code{\link[HiCDOC]{HiCDOC}} for the recommended pipeline.
#' @export
normalizeBiologicalBiases <- function(object) {
    normalizedInteractions <-
        purrr::map_dfr(
            object@chromosomes,
            function(chromosomeName) {
                normalizeBiologicalBiasesOfChromosome(object, chromosomeName)
            }
        ) %>%
        dplyr::mutate(
            chromosome = factor(chromosome, levels = object@chromosomes)
        ) %>%
        dplyr::mutate(condition = factor(condition)) %>%
        dplyr::mutate(replicate = factor(replicate))
    object@interactions <- normalizedInteractions
    return(object)
}
