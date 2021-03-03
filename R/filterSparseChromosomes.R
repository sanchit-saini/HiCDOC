#' sparsityChromosome
#' Compute the maximum of sparsity of a given chromosome
#'
#' @param object a HiCDOCDataSet object
#' @param chromosomeName the name of the chromosome
#'
#' @return a one line tibble
#' @keywords internal
sparsityChromosome <- function(
    object,
    chromosomeName
) {
    interactionsChr <-
        object@interactions %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::filter(interaction > 0)
    totalCells <- object@totalBins[chromosomeName]^2

    pctFill <-
        interactionsChr %>%
        dplyr::mutate(interaction = ifelse(bin.1 == bin.2, 1, 2)) %>%
        dplyr::group_by(replicate, condition) %>%
        dplyr::summarise(
            pctSparse = 1 - (0 + sum(interaction)) / totalCells,
            .groups = "keep"
        ) %>%
        dplyr::ungroup()

    pctFill %<>%
        dplyr::summarise(maxSparsity = max(pctSparse, na.rm = TRUE)) %>%
        dplyr::mutate(
            totalBins = object@totalBins[chromosomeName],
            totalCells = totalCells,
            chromosome = factor(chromosomeName, levels = object@chromosomes)
        ) %>%
        dplyr::mutate(quality = dplyr::case_when(
            maxSparsity < 0.001 ~ "***",
            maxSparsity < 0.01 ~ "**",
            maxSparsity < 0.05 ~ "*",
            maxSparsity < 0.1 ~ ".",
            TRUE ~ ""
        ))

    # reordering columns
    pctFill %<>% dplyr::select(chromosome, totalCells, maxSparsity, quality)
    return(pctFill)
}

#' #filterSparseChromosomes
#' The function compute the sparsity of the matrix, by condition and
#' replicate, for each chromosome. Then it remove (if \code{removeChromosomes}
#' is TRUE) the chromosomes with a sparcity >= \code{threshold}
#' on at least one condition and one replicate.
#' @param object A HiCDOCDataSet object
#' @param threshold A numerical value. The sparsity threshold
#' from which a chromosome will be removed, if \code{removeChromosomes} is TRUE.
#' If NULL, default to the first not NULL of \code{object$sparseReplicateThreshold.} and
#' \code{HiCDOCDefaultParameters$sparseReplicateThreshold.}.
#' @param removeChromosomes Logical, default to TRUE. Should the function
#' remove sparse chromosomes ?
#'
#' @return A HiCDOCDataSet object
#' @export
#' @seealso \code{\link[HiCDOC]{filterSmallChromosomes}},
#' \code{\link[HiCDOC]{filterWeakPositions}} and
#' \code{\link[HiCDOC]{HiCDOC}} for the recommended pipeline.
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSparseChromosomes(object)
filterSparseChromosomes <- function(
    object,
    threshold = NULL,
    removeChromosomes = TRUE
) {
    # Parameters
    if (!is.null(threshold)) {
        object@parameters$sparseReplicateThreshold <- threshold
    }
    object@parameters <- validateParameters(object@parameters)
    thresh <- object@parameters$sparseReplicateThreshold
    # Sparsity
    chrQuality <-
        purrr::map_dfr(
            object@chromosomes,
            function(x) sparsityChromosome(object, x)
        )

    if (removeChromosomes) {
        sparseChromosomes <-
            chrQuality %>%
            dplyr::filter(maxSparsity >= thresh | is.infinite(maxSparsity)) %>%
            dplyr::pull(chromosome) %>%
            as.character()
        if (length(sparseChromosomes) == 0) {
            message(
                "No chromosome removed (threshold: ",
                round(100 * thresh, 2),
                "%)"
            )
        } else {
            message(
                length(sparseChromosomes),
                " chromosome(s) removed: ",
                sparseChromosomes,
                " (threshold : ",
                round(100 * thresh, 2),
                "%)"
            )
            nonSparseChr <-
                object@chromosomes[
                    !(object@chromosomes %in% sparseChromosomes)
                ]
            object <-
                reduceHiCDOCDataSet(object, chromosomes = nonSparseChr)
        }
    }
    chrQuality %<>%
        dplyr::mutate(
            maxSparsity = paste0(round(100 * maxSparsity, 2), "%")
        )
    print(data.frame(chrQuality)) # In data.frame for correct alignement
    cat("\n")
    return(object)
}
