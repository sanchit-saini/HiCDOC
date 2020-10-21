KR <- function(A,
               tol = 1e-6,
               delta = 0.1,
               Delta = 3) {
    n <- nrow(A)
    e <- matrix(1, nrow = n, ncol = 1)
    x0 <- e
    # Inner stopping criterior
    g <- 0.9
    etamax <- 0.1
    eta <- etamax
    stop_tol <- tol * .5
    x <- x0
    rt <- tol ^ 2
    v <- x * (A %*% x)
    rk <- 1 - v
    rho_km1 <- drop(t(rk) %*% rk)
    rout <- rho_km1
    rold <- rout
    MVP <- 0 # We'll count matrix vector products.
    i <- 0 # Outer iteration count.
    while (rout > rt) {
        # Outer iteration
        i <- i + 1
        k <- 0
        y <- e
        innertol <- max(c(eta ^ 2 * rout, rt))
        while (rho_km1 > innertol) {
            # Inner iteration by CG
            k <- k + 1
            if (k == 1) {
                Z <- rk / v
                p <- Z
                rho_km1 <- drop(t(rk) %*% Z)
            }
            else {
                beta <- rho_km1 / rho_km2
                p <- Z + beta * p
            }

            # Update search direction efficiently
            w <- x * (A %*% (x * p)) + v * p
            if (max(w) == Inf) {
                message("Warning: KR algorithm diverges.")
                return(t(t(x[, 1] * A) * x[, 1]))
            }
            alpha <- rho_km1 / drop(t(p) %*% w)
            ap <- alpha * p
            # Test distance to boundary of cone
            ynew <- y + ap
            if (min(ynew) <= delta) {
                if (delta == 0)
                    break()
                ind <- which(ap < 0)
                gamma <- min((delta - y[ind]) / ap[ind])
                y <- y + gamma * ap
                break()
            }
            if (max(ynew) >= Delta) {
                ind <- which(ynew > Delta)
                gamma <- min((Delta - y[ind]) / ap[ind])
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
        MVP <- MVP + k + 1
        # Update inner iteration stopping criterion.
        rat <- rout / rold
        rold <- rout
        res_norm <- sqrt(rout)
        eta_o <- eta
        eta <- g * rat
        if (g * eta_o ^ 2 > 0.1) {
            eta <- max(c(eta, g * eta_o ^ 2))
        }
        eta <- max(c(min(c(eta, etamax)), stop_tol / res_norm))
    }

    result <- t(t(x[, 1] * A) * x[, 1])

    return(result)
}


##- normalizeBiologicalBiasesChr ---------------------------------------------#
##----------------------------------------------------------------------------#
#' Remove biological biases by normalizing with Knight-Ruiz matrix balancing.
#'
#' @rdname normalizeBiologicalBiasesChr
#'
#' @param object A \code{HiCDOCExp} object.
#' @param chromosomeId The name or number of the chromosome to plot.
#' If number, will be taken in \code{object@chromosomes[chromosomeId]}
#'
#' @return A \code{HiCDOCExp} object, with the normalized matrices.
normalizeBiologicalBiasesChr <- function(object, chromosomeId) {
    testSlotsHiCDOCExp(
        object,
        slots = c(
            "interactions",
            "totalBins",
            "conditions",
            "replicates",
            "totalBins"
        )
    )
    chr <- testChromosome(object, chromosomeId)
    message("Chromosome: ", chr)

    if (object@totalBins[[chr]] == -Inf)
        return(NULL)

    rawMatrices    <-
        mapply(
            function(x, y)
                sparseInteractionsToMatrix(object, chr, x, y, filter = TRUE),
            object@conditions,
            object@replicates,
            SIMPLIFY = FALSE
        )

    testEmptyRow <-
        vapply(rawMatrices, function(x)
            min(rowSums(x)), FUN.VALUE = c(0))
    if (length(testEmptyRow[testEmptyRow == 0]) > 0) {
        id <- which(testEmptyRow == 0)
        stop(
            "A matrix has an empty row",
            paste0(
                "condition ",
                object@conditions[id],
                " replicate ",
                object@replicates[id]
            )
        )
    }

    normalizedMatrices <- lapply(rawMatrices, KR)
    sparseIntMatrices <-
        mapply(
            function(x, y, z)
                matrixToSparseInteractions(x, object, chr, y, z),
            normalizedMatrices,
            object@conditions,
            object@replicates,
            SIMPLIFY = FALSE
        )

    interactionsChr <- dplyr::bind_rows(sparseIntMatrices)
    return(interactionsChr)
}

##- normalizeBiologicalBiases ------------------------------------------------#
##----------------------------------------------------------------------------#
#' Remove biological biases by normalizing with Knight-Ruiz matrix balancing.
#'
#' @rdname normalizeBiologicalBiases
#'
#' @param object A \code{HiCDOCExp} object.
#'
#' @return A \code{HiCDOCExp} object, with the normalized matrices.
#'
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' @export
normalizeBiologicalBiases <- function(object) {
    interactionsNorm <-
        purrr::map_dfr(object@chromosomes,
                       function(x) normalizeBiologicalBiasesChr(object, x)) %>%
        dplyr::mutate(chromosome = factor(chromosome, levels = object@chromosomes)) %>%
        dplyr::mutate(condition = factor(condition)) %>%
        dplyr::mutate(replicate = factor(replicate))
    object@interactions <- interactionsNorm
    return(object)
}
