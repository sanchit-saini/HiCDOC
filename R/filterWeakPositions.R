#' Look for the weak positions for a given chromosome
#'
#' The function identifies and remove the weak positions
#' from the interaction matrix for a given chromosome, in a recursive way.
#' To be kept, the bins (positions) must have a mean value greater than
#' \code{threshold}, on all replicates and all conditions.
#' The mean is computed on the row of the reconstructed full interaction
#' matrix for the chromosome, 1 condition and 1 replicate.
#'
#' @param object Interaction matrix for the chromosome
#' @param chromosomeName A character or numeric value.
#' Name or number of the chromosome
#' @param threshold Numeric default to 0.
#'
#' @return list of length 2 : \code{"pos"} = the weak positions,
#' \code{"interactions"} the interactions matrix for the chromosome,
#' whithout the weak bins.
.filterWeakPositionsOfChromosome <- function(
    object,
    chromosomeName,
    threshold = 0
) {
    # Initialization
    interactions <-
        object@interactions %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::filter(interaction > 0)

    totalBins <- object@totalBins[[chromosomeName]]
    allBins <- seq_len(totalBins)
    removedBins <-
        allBins[
            !(allBins %in% unique(c(interactions$bin.1, interactions$bin.2)))
        ]

    totalNewWeakBins <- 1
    totalRemovedBins <- 0

    # Recursive removal of bins - deleting a bin can create a new weak bin.
    while (totalNewWeakBins > 0 & totalRemovedBins <= totalBins) {
        diagonalInteractions <-
            interactions %>%
            dplyr::filter(bin.1 == bin.2) %>%
            dplyr::rename(bin = bin.1) %>%
            dplyr::select(-chromosome, -bin.2)

        rowInteractions <-
            interactions %>%
            dplyr::filter(bin.1 != bin.2) %>%
            tidyr::pivot_longer(
                cols = c(`bin.1`, `bin.2`),
                values_to = "bin"
            ) %>%
            dplyr::select(-name, -chromosome) %>%
            dplyr::bind_rows(diagonalInteractions) %>%
            tidyr::complete(
                bin,
                tidyr::nesting(condition, replicate),
                fill = list(interaction = 0)
            )

        weakBins <-
            rowInteractions %>%
            dplyr::group_by(replicate, condition, bin) %>%
            dplyr::mutate(mean = sum(interaction) / totalBins) %>%
            dplyr::filter(mean <= threshold) %>%
            dplyr::pull(bin) %>%
            unique() %>%
            sort()

        totalNewWeakBins <- length(weakBins) - totalRemovedBins
        removedBins <- c(removedBins, weakBins)

        # Remove interactions of weak bins
        if (totalNewWeakBins > 0) {
            interactions <- interactions[!(interactions$bin.1 %in% weakBins), ]
            interactions <- interactions[!(interactions$bin.2 %in% weakBins), ]
            totalRemovedBins <- totalRemovedBins + totalNewWeakBins
        }
    }

    message(
        "Chromosome ",
        chromosomeName,
        ": ",
        length(removedBins),
        " position",
        if (length(removedBins) != 1) "s",
        " removed, ",
        length(allBins) - length(removedBins),
        " position",
        if (length(allBins) - length(removedBins) != 1) "s",
        " remaining."
    )
    return(list("removedBins" = removedBins, "interactions" = interactions))
}


#' Remove the weak bins of an HiCDOCDataSet object
#'
#' The function indentifies and remove the weak bins of an HiCDOCDataSet object,
#' chromosome by chromosome.
#' To be kept, the bins (positions)of a chromosome must have a mean value
#' greater than \code{threshold}, on all replicates and all conditions.
#' The mean is computed on the row of the reconstructed full interaction
#' matrix for 1 chromosome, 1 condition and 1 replicate. The weak bins will
#' be discarded, and added in object@weakBins.
#'
#' @param object A HiCDOCDataSet object
#' @param threshold Numeric. Threshold to consider bins as 'weaks'. If NULL,
#' default to the first not NULL of \code{object$weakPositionThreshold} and
#' \code{HiCDOCDefaultParameters$weakPositionThreshold}.
#'
#' @return A \code{HiCDOCDataSet} object with a reduced \code{interactions}
#' slot, and the weak bins identified by chromosomes in \code{object@weakBins}.
#' @export
#'
#' @seealso \code{\link[HiCDOC]{filterSmallChromosomes}},
#' \code{\link[HiCDOC]{filterSparseChromosomes}} and
#' \code{\link[HiCDOC]{HiCDOC}} for the recommended pipeline.
#' @examples
#' exp <- HiCDOCExample()
#' exp <- filterWeakPositions(exp)
filterWeakPositions <- function(object, threshold = NULL) {
    .validateSlots(
        object,
        slots = c(
            "interactions",
            "binSize",
            "chromosomes",
            "replicates",
            "conditions",
            "totalBins",
            "parameters"
        )
    )
    if (!is.null(threshold)) {
        object@parameters$weakPositionThreshold <- threshold
    }
    object@parameters <- .validateParameters(object@parameters)

    message(
        "Keeping positions with interactions average greater than ",
        object@parameters$weakPositionThreshold,
        "."
    )

    results <-
        lapply(
            object@chromosomes,
            function(chromosomeName) {
                .filterWeakPositionsOfChromosome(
                    object,
                    chromosomeName,
                    object@parameters$weakPositionThreshold
                )
            }
        )
    names(results) <- object@chromosomes
    weakBins <- results %>% purrr::map("removedBins")
    intactChromosomes <-
        vapply(
            weakBins,
            function(x) length(x) == 0,
            FUN.VALUE = TRUE
        )
    weakBins[intactChromosomes] <- list(NULL)
    interactions <- results %>% purrr::map_dfr("interactions")

    object@weakBins <- weakBins
    object@interactions <- interactions
    totalWeakBins <- sum(vapply(weakBins, length, FUN.VALUE = 0))
    message(
        "Removed ",
        totalWeakBins,
        " position",
        if (totalWeakBins != 1) "s",
        " in total."
    )

    if (totalWeakBins >= sum(unlist(object@totalBins))) message("No data left!")

    return(object)
}
