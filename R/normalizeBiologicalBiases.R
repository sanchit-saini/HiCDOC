#" Based on
#" https://github.com/dozmorovlab/HiCcompare/blob/master/R/KRnormalization.R
#"
#' @description
#' Applies the Knight-Ruiz balancing algorithm to transform the provided
#' matrix into a doubly stochastic matrix, with sum of each row and sum of each
#' column equal to 1.
#'
#' @param m
#' A matrix.
#'
#' @return
#' The transformed matrix.
#'
#' @keywords internal
#' @noRd
.normalizeKnightRuiz <- function(cm) {
    m <- cm@matrix
    m[is.na(m)] <- 0
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
    cm@matrix <- result
    return(cm)
}

#' @description
#' Normalizes biological biases in the interactions of a given chromosome. Calls
#' \code{.normalizeKnightRuiz} internally.
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
.normalizeBiologicalBiasesOfChromosome <- function(object) {
    
    chromosomeName <- as.character(SummarizedExperiment::mcols(object)$Chr[1])
    message("Chromosome ", chromosomeName, ": normalizing biological biases.")
    if (object@totalBins[[chromosomeName]] <= 0) return(NULL)
    currentOrder <- InteractionSet::anchorIds(object)
    currentAssay <- SummarizedExperiment::assay(object)
    
    # Pass by InteractionSet so we can use inflate/deflate
    isetChromosome <- InteractionSet::InteractionSet(
        currentAssay,
        InteractionSet::interactions(object)
    )
    validAssay <- object@validAssay[[chromosomeName]]
    
    matrices <- 
        lapply(validAssay,
               FUN = function(x) {
                   InteractionSet::inflate(isetChromosome, 
                                           rows = chromosomeName,
                                           columns = chromosomeName,
                                           sample = x,
                                           sparse = FALSE)
                   })
    matrices <- lapply(matrices, function(m) {
        m@matrix[is.na(m@matrix)] <- 0
        return(m)
        })
    
    matrices <- lapply(matrices, .normalizeKnightRuiz)
    matrices <- lapply(matrices, function(m) {
        m@matrix[m@matrix == 0] <- NA
        return(m)
    })
    matrices <- lapply(matrices, InteractionSet::deflate, use.na=TRUE)

    ids <- InteractionSet::anchorIds(matrices[[1]])
    ids <- paste(ids$first, ids$second)
    correctIds <- paste(currentOrder$first, currentOrder$second)
    orderids <- match(correctIds, ids)
    matrices <- lapply(matrices, function(x) SummarizedExperiment::assay(x))
    matrices <- lapply(matrices, function(x) x[orderids,])
    matrices <- do.call(base::"cbind", matrices)
    currentAssay[,validAssay] <- matrices
    return(currentAssay)
}

#' @title
#' Normalize biological biases.
#'
#' @description
#' Normalizes biological biases such as GC content and repeated regions. Uses
#' the Knight-Ruiz balancing algorithm to transform interaction matrices into
#' doubly stochastic matrices, with sum of each row and sum of each column equal
#' to 1.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}} with normalized interactions.
#'
#' @examples
#' data(exampleHiCDOCDataSet)
#' object <- exampleHiCDOCDataSet
#' object <- filterSparseReplicates(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeBiologicalBiases(object)
#'
#' @seealso
#' \code{\link{filterSparseReplicates}},
#' \code{\link{filterWeakPositions}},
#' \code{\link{normalizeTechnicalBiases}},
#' \code{\link{normalizeDistanceEffect}},
#' \code{\link{HiCDOC}}
#'
#' @export
normalizeBiologicalBiases <- function(object, parallel = FALSE) {

    .validateSlots(
        object,
        slots = c(
            "totalBins",
            "validAssay"
        )
    )
    objectChromosomes <- S4Vectors::split(
        object, 
        SummarizedExperiment::mcols(object)$Chr, drop=FALSE)
    
    normAssay <- .internalLapply(parallel,
                                objectChromosomes,
                                FUN = .normalizeBiologicalBiasesOfChromosome)
    normAssay <- do.call("rbind", normAssay)

    SummarizedExperiment::assay(object) <- normAssay

    return(object)
}
