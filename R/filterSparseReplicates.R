#' @description
#' Removes sparse replicates of a given chromosome.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' The name of a chromosome.
#' @param threshold
#' The minimum percentage of non-zero interactions for a replicate to be kept.
#'
#' @return
#' A list of:
#' - The sparse condition names repeated along the sparse replicates.
#' - The sparse replicate names repeated along the sparse conditions.
#' - The filtered interactions.
#'
#' @md
#' @keywords internal
#' @noRd
.filterSparseReplicatesOfChromosome <- function(
    object,
    chromosomeName,
    threshold
) {

    chromosomeInteractions <-
        object@interactions %>%
        dplyr::filter(chromosome == chromosomeName)

    presentInteractions <-
        chromosomeInteractions %>%
        dplyr::filter(interaction > 0) %>%
        dplyr::mutate(total.interactions = ifelse(bin.1 == bin.2, 1, 2)) %>%
        dplyr::select(-interaction)

    totalBins <- object@totalBins[chromosomeName]
    totalCells <- totalBins^2

    allInteractions <-
        dplyr::tibble(
            chromosome = factor(chromosomeName, levels = object@chromosomes),
            condition = rep(object@conditions, each = totalCells),
            replicate = rep(object@replicates, each = totalCells),
            bin.1 = rep(
                seq(totalBins),
                each = totalBins,
                times = length(object@conditions)
            ),
            bin.2 = rep(
                seq(totalBins),
                times = totalBins * length(object@conditions)
            ),
            total.interactions = 0
        ) %>%
        dplyr::filter(bin.2 >= bin.1) %>%
        dplyr::left_join(
            presentInteractions,
            by = c("chromosome", "condition", "replicate", "bin.1", "bin.2")
        ) %>%
        dplyr::mutate(
            total.interactions = dplyr::coalesce(
                total.interactions.y,
                total.interactions.x
            )
        ) %>%
        dplyr::select(-total.interactions.y, -total.interactions.x)

    fillPercentages <-
        allInteractions %>%
        dplyr::group_by(condition, replicate) %>%
        dplyr::summarise(
            fill = sum(total.interactions) / totalCells,
            .groups = "keep"
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            chromosome = factor(chromosomeName, levels = object@chromosomes)
        )

    chromosomeInteractions %<>%
        dplyr::left_join(
            fillPercentages,
            by = c("chromosome", "condition", "replicate")
        ) %>%
        dplyr::filter(fill >= threshold) %>%
        dplyr::filter(interaction > 0) %>%
        dplyr::select(-fill)

    removedReplicates <- fillPercentages %>% dplyr::filter(fill < threshold)

    if (nrow(removedReplicates) > 0) {
        message(
            paste0(
                "Removed interactions matrix of chromosome ",
                removedReplicates$chromosome,
                ", condition ",
                removedReplicates$condition,
                ", replicate ",
                removedReplicates$replicate,
                " filled at ",
                round(removedReplicates$fill, digits = 5) * 100,
                "%.",
                collapse = "\n"
            )
        )
    }

    sparseConditions <- removedReplicates %>% dplyr::pull(condition)
    sparseReplicates <- removedReplicates %>% dplyr::pull(replicate)

    return(list(
        "sparseConditions" = sparseConditions,
        "sparseReplicates" = sparseReplicates,
        "interactions" = chromosomeInteractions
    ))
}

#' @title
#' Filter sparse replicates.
#'
#' @description
#' Removes chromosome replicates whose percentage of non-zero interactions is
#' smaller than the threshold.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param threshold
#' The minimum percentage of non-zero interactions for a chromosome replicate to
#' be kept. If a chromosome replicate's percentage of non-zero interactions is
#' lower than this value, it is removed. Defaults to
#' \code{object$smallChromosomeThreshold} which is originally set to
#' \code{defaultHiCDOCParameters$smallChromosomeThreshold} = 0.05.
#'
#' @return
#' A filtered \code{\link{HiCDOCDataSet}}.
#'
#' @seealso
#' \code{\link{filterSmallChromosomes}},
#' \code{\link{filterWeakPositions}},
#' \code{\link{HiCDOC}}
#'
#' @examples
#' object <- HiCDOCDataSetExample()
#' object <- filterSparseReplicates(object)
#'
#' @usage
#' filterSparseReplicates(object, threshold = 0.05)
#'
#' @export
filterSparseReplicates <- function(
    object,
    threshold = NULL
) {

    if (!is.null(threshold)) {
        object@parameters$sparseReplicateThreshold <- threshold
    }
    object@parameters <- .validateParameters(object@parameters)
    threshold <- object@parameters$sparseReplicateThreshold

    message(
        "Keeping replicates filled with at least ",
        threshold * 100,
        "% non-zero interactions."
    )

    results <-
        lapply(
            object@chromosomes,
            function(chromosomeName) {
                .filterSparseReplicatesOfChromosome(
                    object,
                    chromosomeName,
                    threshold
                )
            }
        )

    names(results) <- object@chromosomes

    sparseConditions <- results %>% purrr::map("sparseConditions")
    sparseReplicates <- results %>% purrr::map("sparseReplicates")
    intactChromosomes <-
        vapply(
            sparseReplicates,
            function(x) length(x) == 0,
            FUN.VALUE = TRUE
        )
    sparseConditions[intactChromosomes] <- list(NULL)
    sparseReplicates[intactChromosomes] <- list(NULL)

    interactions <- results %>% purrr::map_dfr("interactions")

    object@sparseConditions <- sparseConditions
    object@sparseReplicates <- sparseReplicates
    object@interactions <- interactions

    totalSparseReplicates <-
        sum(vapply(sparseReplicates, length, FUN.VALUE = 0))

    message(
        "Removed ",
        totalSparseReplicates,
        " replicate",
        if (totalSparseReplicates != 1) "s",
        " in total."
    )

    if (nrow(object@interactions) == 0) {
        warning("No data left!", call. = FALSE)
    }

    return(object)
}
