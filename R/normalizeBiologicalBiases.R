#" Based on
#" https://github.com/dozmorovlab/HiCcompare/blob/master/R/KRnormalization.R
#"
#' @description
#' Applies the Knight-Ruiz balancing algorithm on the provided matrix. The
#' matrix is transformed such that the sum of each row equals 1 and the sum of
#' each column equals 1.
#'
#' @param m
#' A matrix.
#'
#' @return
#' The transformed matrix.
#'
#' @keywords internal
#' @noRd
.normalizeKnightRuiz <- function(m) {
    tol <- 1e-6
    minDelta <- 0.1
    maxDelta <- 3
    n <- nrow(m)
    e <- matrix(1, nrow = n, ncol = 1)
    x0 <- e
    g <- 0.9
    etamax <- 0.1
    eta <- etamax
    stop_tol <- tol * .5
    x <- x0
    rt <- tol^2
    v <- x * (m %*% x)
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
            w <- x * (m %*% (x * p)) + v * p
            if (max(w) == Inf) {
                warning("KR algorithm diverges.", call. = FALSE)
                return(t(t(x[, 1] * m) * x[, 1]))
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
        v <- x * (m %*% x)
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
    result <- t(t(x[, 1] * m) * x[, 1])
    return(result)
}

#' @description
#' Normalizes the biological biases in the interactions of a given chromosome.
#' Calls \code{.normalizeKnightRuiz} internally.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' The name of a chromosome to normalize.
#'
#' @return
#' A tibble of normalized interactions.
#'
#' @keywords internal
#' @noRd
.normalizeBiologicalBiasesOfChromosome <- function(object, chromosomeName) {

    message("Chromosome ", chromosomeName, ": normalizing biological biases.")

    if (object@totalBins[[chromosomeName]] == -Inf) return(NULL)

    matrices <-
        mapply(
            function(conditionName, replicateName) {
                .sparseInteractionsToMatrix(
                    object,
                    chromosomeName,
                    conditionName,
                    replicateName,
                    filter = TRUE
                )
            },
            object@validConditions[[chromosomeName]],
            object@validReplicates[[chromosomeName]],
            SIMPLIFY = FALSE
        )

    normalizedMatrices <- lapply(matrices, .normalizeKnightRuiz)

    interactions <-
        purrr::pmap_dfr(
            list(
                normalizedMatrices,
                object@validConditions[[chromosomeName]],
                object@validReplicates[[chromosomeName]]
            ),
            .f = function(matrix, conditionName, replicateName) {
                .matrixToSparseInteractions(
                    matrix,
                    object,
                    chromosomeName,
                    conditionName,
                    replicateName
                )
            }
        )
    return(interactions)
}

#' @title
#' Normalize biological biases.
#'
#' @description
#' Normalizes biological biases such as GC content and repeated regions. Uses
#' the Knight-Ruiz balancing algorithm to transform interaction matrices into
#' doubly stochastic matrices, with sum of rows and sum of columns equal to 1.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}} with the normalized interactions.
#'
#' @examples
#' object <- HiCDOCDataSetExample()
#' object <- filterSparseReplicates(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeBiologicalBiases(object)
#'
#' @seealso
#' \code{\link{filterSparseReplicates}},
#' \code{\link{filterWeakPositions}},
#' \code{\link{HiCDOC}}
#'
#' @export
normalizeBiologicalBiases <- function(object) {

    .validateSlots(
        object,
        slots = c(
            "interactions",
            "totalBins",
            "weakBins",
            "validReplicates",
            "validConditions"
        )
    )

    normalizedInteractions <-
        purrr::map_dfr(
            object@chromosomes,
            function(chromosomeName) {
                .normalizeBiologicalBiasesOfChromosome(object, chromosomeName)
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
