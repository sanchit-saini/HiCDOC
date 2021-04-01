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
            fillPct = sum(total.interactions) / totalCells,
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
        dplyr::filter(fillPct >= threshold) %>%
        dplyr::filter(interaction > 0) %>%
        dplyr::select(-fillPct)

    removed <- fillPercentages %>% dplyr::filter(fillPct < threshold)

    if (nrow(removed) > 0) {
        message(
            paste0(
                "Removed interactions matrix of chromosome ",
                removed$chromosome,
                ", condition ",
                removed$condition,
                ", replicate ",
                removed$replicate,
                " filled at ",
                round(removed$fillPct, digits = 5) * 100,
                "%.",
                collapse = "\n"
            )
        )
    }

    valid <- fillPercentages %>% dplyr::filter(fillPct >= threshold)

    validConditions <- valid %>% dplyr::pull(condition)
    validReplicates <- valid %>% dplyr::pull(replicate)

    return(list(
        "validConditions" = validConditions,
        "validReplicates" = validReplicates,
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
#' data(HiCDOCDataSetExample)
#' object <- filterSparseReplicates(HiCDOCDataSetExample)
#'
#' @usage
#' filterSparseReplicates(object, threshold = 0.05)
#'
#' @export
filterSparseReplicates <- function(object, threshold = NULL) {

    .validateSlots(
        object,
        slots = c(
            "interactions",
            "chromosomes",
            "conditions",
            "replicates",
            "parameters"
        )
    )

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

    validConditions <- results %>% purrr::map("validConditions")
    validReplicates <- results %>% purrr::map("validReplicates")
    badChromosomes <-
        vapply(
            validReplicates,
            function(x) length(x) == 0,
            FUN.VALUE = TRUE
        )
    validConditions[badChromosomes] <- list(NULL)
    validReplicates[badChromosomes] <- list(NULL)

    interactions <-
        results %>%
        purrr::map_dfr("interactions") %>%
        .sortInteractions(
            object@chromosomes,
            object@conditions,
            object@replicates
        )

    object@validConditions <- validConditions
    object@validReplicates <- validReplicates
    object@interactions <- interactions

    totalSparseReplicates <-
        sum(
            vapply(
                validReplicates,
                function(replicates) {
                    length(object@replicates) - length(replicates)
                },
                FUN.VALUE = 0
            )
        )

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
