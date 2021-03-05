#' @description
#' Removes weak positions of a given chromosome.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' The name of a chromosome.
#' @param threshold
#' The minimum average interaction for a position to be kept.
#'
#' @return
#' A list of:
#' - The weak positions.
#' - The filtered interactions.
#'
#' @md
#' @keywords internal
#' @noRd
.filterWeakPositionsOfChromosome <- function(
    object,
    chromosomeName,
    threshold
) {

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
            dplyr::filter(mean < threshold) %>%
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

#' @title
#' Filter weak positions.
#'
#' @description
#' Removes weak genomic positions whose interactions average is lower than the
#' threshold.
#'
#' @details
#' Detects weak genomic positions in each replicate, and removes them from all
#' replicates to guarantee comparability across conditions when calling
#' \code{\link{detectCompartments}}.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param threshold
#' The minimum average interaction for a position to be kept. If a position's
#' average interaction with the entire chromosome is lower than this value in
#' any of the replicates, it is removed from all replicates and conditions.
#' Defaults to \code{object$smallChromosomeThreshold} which is
#' originally set to \code{defaultHiCDOCParameters$smallChromosomeThreshold}.
#'
#' @return
#' A filtered \code{\link{HiCDOCDataSet}}.
#'
#' @seealso
#' \code{\link{filterSmallChromosomes}},
#' \code{\link{filterSparseReplicates}},
#' \code{\link{HiCDOC}}
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- filterWeakPositions(exp)
#'
#' @usage
#' filterWeakPositions(object, threshold = 1)
#'
#' @export
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
    threshold <- object@parameters$weakPositionThreshold

    message(
        "Keeping positions with interactions average greater or equal to ",
        threshold,
        "."
    )

    results <-
        lapply(
            object@chromosomes,
            function(chromosomeName) {
                .filterWeakPositionsOfChromosome(
                    object,
                    chromosomeName,
                    threshold
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

    if (nrow(object@interactions) == 0) message("No data left!")

    return(object)
}
