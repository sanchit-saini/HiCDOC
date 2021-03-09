#' @description
#' Determines the number of bases per bin by finding the minimum distance
#' between \code{position.1} and \code{position.2} that is not zero.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' The number of bases per bin.
#'
#' @keywords internal
#' @noRd
.determineBinSize <- function(object) {

    binSize <-
        min(abs(
            object@interactions$position.2[
                object@interactions$position.1 != object@interactions$position.2
            ] -
            object@interactions$position.1[
                object@interactions$position.1 != object@interactions$position.2
            ]
        ))

    return(binSize)
}

#' @description
#' Determines the position corresponding to each bin of each chromosome.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A tibble of 4 columns: chromosome, bin, start, end.
#'
#' @keywords internal
#' @noRd
.determinePositions <- function(object) {

    maxPositions <-
        object@interactions %>%
        dplyr::group_by(chromosome) %>%
        dplyr::summarize(
            maxPosition.1 = max(position.1),
            maxPosition.2 = max(position.2)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            maxPosition = max(maxPosition.1, maxPosition.2),
            .keep = "unused"
        )

    positions <-
        lapply(
            object@chromosomes,
            function(chromosomeName) {
                return(data.frame(
                    "chromosome" = chromosomeName,
                    "start" = seq(
                        0,
                        maxPositions[
                            maxPositions$chromosome == chromosomeName,
                        ]$maxPosition,
                        by = object@binSize
                    )
                ))
            }
        ) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(
            chromosome = factor(chromosome, levels = object@chromosomes)
        ) %>%
        dplyr::group_by(chromosome) %>%
        dplyr::arrange(chromosome, start, .by_group = TRUE) %>%
        dplyr::mutate(end = dplyr::lead(start) - 1) %>%
        dplyr::mutate(bin = dplyr::row_number()) %>%
        dplyr::ungroup()

    positions %<>%
        dplyr::mutate(
            end = ifelse(is.na(end), start + object@binSize - 1, end)
        ) %>%
        dplyr::select(chromosome, bin, start, end)

    return(positions)
}

#' @description
#' Replaces positions by their corresponding bins.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A tibble of interactions with bin columns replacing position columns.
#'
#' @keywords internal
#' @noRd
.replacePositionsByBins <- function(object) {

    interactions <-
        object@interactions %>%
        dplyr::left_join(
            object@positions %>% dplyr::select(
                chromosome,
                position.1 = start,
                bin.1 = bin
            ),
            by = c("chromosome", "position.1")
        ) %>%
        dplyr::left_join(
            object@positions %>% dplyr::select(
                chromosome,
                position.2 = start,
                bin.2 = bin
            ),
            by = c("chromosome", "position.2")
        ) %>%
        dplyr::select(
            chromosome,
            condition,
            replicate,
            bin.1,
            bin.2,
            interaction
        )

    return(interactions)
}

#' @description
#' Determines the number of bins per chromosome.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A list of the number of bins per chromosome.
#'
#' @keywords internal
#' @noRd
.determineChromosomeSizes <- function(object) {

    totalBins <- vapply(
        object@chromosomes,
        function(chromosomeName) {
            nrow(
                object@positions[
                    object@positions$chromosome == chromosomeName,
                ]
            )
        },
        FUN.VALUE = 0
    )

    return(totalBins)
}

#' @description
#' Reformats the interactions to fit them into a numeric upper triangular
#' sparse matrix. If duplicate interactions are found, they are averaged.
#'
#' @param interactions
#' A tibble of interactions.
#'
#' @return
#' A reformated tibble of interactions.
#'
#' @keywords internal
#' @noRd
.reformatInteractions <- function(interactions) {

    lowerTriangle <- which(interactions$bin.1 > interactions$bin.2)
    if (length(lowerTriangle) > 0) {
        lowerTriangleBins1 <-
            interactions[lowerTriangle, "bin.1"] %>%
            dplyr::pull()
        lowerTriangleBins2 <-
            interactions[lowerTriangle, "bin.2"] %>%
            dplyr::pull()
        interactions[lowerTriangle, ]$bin.1 <- lowerTriangleBins2
        interactions[lowerTriangle, ]$bin.2 <- lowerTriangleBins1
    }

    if (!is.numeric(interactions$interaction)) {
        interactions$interaction <- as.numeric(interactions$interaction)
    }

    if (
        nrow(interactions) != nrow(unique(
            interactions[,
                c("chromosome", "condition", "replicate", "bin.1", "bin.2")
            ]
        ))
    ) {
        message("Averaging duplicate interactions.")
        interactions %<>%
            dplyr::group_by(chromosome, condition, replicate, bin.1, bin.2) %>%
            dplyr::mutate(interaction = mean(interaction)) %>%
            dplyr::ungroup()
    }

    interactions %<>% dplyr::filter(interaction > 0)

    return(interactions)
}

#' @description
#' Sorts interactions and factorizes chromosomes, conditions, and replicates.
#'
#' @param interactions
#' A tibble of interactions.
#' @param chromosomeNames
#' The names of chromosomes to sort with.
#' @param conditionNames
#' The names of conditions to sort with.
#' @param replicateNames
#' The names of replicates to sort with.
#'
#' @return
#' A sorted tibble of interactions.
#'
#' @keywords internal
#' @noRd
.sortInteractions <- function(
    interactions,
    chromosomeNames,
    conditionNames,
    replicateNames
) {

    interactions %<>%
        dplyr::mutate(
            chromosome = factor(
                chromosome,
                levels = gtools::mixedsort(unique(chromosomeNames))
            ),
            condition = factor(
                condition,
                levels = gtools::mixedsort(unique(conditionNames))
            ),
            replicate = factor(
                replicate,
                levels = gtools::mixedsort(unique(replicateNames))
            ),
            interaction = as.numeric(interaction)
        ) %>%
        dplyr::filter(interaction > 0)

    if ("bin.1" %in% colnames(interactions)) {
        interactions %<>%
            dplyr::arrange(
                chromosome,
                condition,
                replicate,
                bin.1,
                bin.2
            )
    } else if ("position.1" %in% colnames(interactions)) {
        interactions %<>%
            dplyr::arrange(
                chromosome,
                condition,
                replicate,
                position.1,
                position.2
            )
    } else {
        interactions %<>%
            dplyr::arrange(
                chromosome,
                condition,
                replicate
            )
    }

    return(interactions)
}

#' @description
#' Fills parameters and slots describing the data. Called by a
#' \code{\link{HiCDOCDataSet}} constructor.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}} ready for analysis.
#'
#' @keywords internal
#' @noRd
.fillHiCDOCDataSet <- function(object) {

    .validateSlots(object, "interactions")

    object@chromosomes <-
        gtools::mixedsort(unique(as.character(object@interactions$chromosome)))

    object@interactions %<>%
        .sortInteractions(
            object@chromosomes,
            object@conditions,
            object@replicates
        )

    if ('position.1' %in% colnames(object@interactions)) {
        object@binSize <- .determineBinSize(object)
        object@positions <- .determinePositions(object)
        object@interactions <- .replacePositionsByBins(object)
    }

    if (!is.null(object@positions)) {
        object@totalBins <- .determineChromosomeSizes(object)
    }

    object@interactions <- .reformatInteractions(object@interactions)

    object@parameters <- defaultHiCDOCParameters

    return(object)
}
