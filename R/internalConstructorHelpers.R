# Helpers for HiCDOCDataSet constructor
# For internal use only
# These functions are called by the .fillHiCDOCDataSet() function

#### .computeBinSize ##########################################################
#' Compute the binSize
#'
#' The function takes a HiCDOCDataset object and computes the binSize from
#' the positions in the interactions matrix.
#' It takes the minimum of difference between position.1 and position.2
#' @param object a HiCDOCDataSet object.
#'
#' @return binSize, an integer.
#' @keywords internal
#' @noRd
.computeBinSize <- function(object) {

    if (is.null(object@interactions)) return(NULL)

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

#### .determinePositions ######################################################
#' Determine the positions corresponding to the bins
#'
#' The function determines, for each chromosome, the positions corresponding
#' to the bins from the interaction matrix and the binSize. The positions will
#' be stored in positions(object).
#'
#' @param object a HiCDOCDataSet
#'
#' @return a tibble, with 4 columns : chromosome, bin, start, end.
#' @keywords internal
#' @noRd
.determinePositions <- function(object) {
    .validateSlots(object, c("chromosomes", "interactions"))
    chromosomes <- object@chromosomes

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
            chromosomes,
            function(x) {
                return(data.frame(
                    "chromosome" = x,
                    "start" = seq(
                        0,
                        maxPositions[
                            maxPositions$chromosome == x,
                        ]$maxPosition,
                        by = object@binSize
                    )
                ))
            }
        ) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(chromosome = factor(chromosome, levels = chromosomes)) %>%
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

#### .replacePositionsByBins ##################################################
#' Replace the positions by the corresponding bins in the interaction matrix
#'
#' @param object a HiCDOCDataSet object
#'
#' @return a tibble, of the modified interactions matrix.
#' @keywords internal
#' @noRd
.replacePositionsByBins <- function(object) {
    interactions <-
        object@interactions %>%
        dplyr::mutate(
            chromosome = factor(
                chromosome,
                levels = object@chromosomes
            )
        ) %>%
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
            bin.1,
            bin.2,
            condition,
            replicate,
            interaction
        )
    return(interactions)
}

#### .reformatInteractions ##################################################
#' Reformat the interaction matrix
#'
#'
#' The function takes an interaction matrix and performs some correction
#' - the matrix is putted in upper format only (the values in the lower side
#' will be reaffected in the upper side)
#' - the interaction column is transformed in numeric
#' - the duplicated combinations of positions are aggregated by mean
#'
#' @param interactions the slot interactions of a HiCDOCDataSet object.
#'
#' @return a tibble, of the interactions corrected
#' @keywords internal
#' @noRd
.reformatInteractions <- function(interactions) {
    # Move lower triangle interactions to upper triangle
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
    # Cast values to numeric
    if (!is.numeric(interactions$interaction)) {
        interactions$interaction <- as.numeric(interactions$interaction)
    }
    # Average duplicate interactions
    if (
        nrow(interactions) != nrow(unique(
            interactions[,
                c("chromosome", "bin.1", "bin.2", "condition", "replicate")
            ]
        ))
    ) {
        message("Averaging duplicate interactions.")
        interactions %<>%
            dplyr::group_by(chromosome, bin.1, bin.2, condition, replicate) %>%
            dplyr::mutate(interaction = mean(interaction)) %>%
            dplyr::ungroup()
    }
    return(interactions)
}

## - HiCDOCDataSet S4 class constructor -------------------------------------#
## --------------------------------------------------------------------------#
#' Constructor for the HiCDOCDataSet class.
#'
#' This function should not be called directly. It is called by the functions
#' \code{\link{HiCDOCDataSetFromTabular}},
#' \code{\link{HiCDOCDataSetFromCool}},
#' \code{\link{HiCDOCDataSetFromHiC}} and
#' \code{\link{HiCDOCDataSetFromHiCPro}}.
#' @name HiCDOCDataSet-constructor
#' @rdname HiCDOCDataSet-constructor
#' @aliases HiCDOCDataSet-constructor
#' @param object A prefilled \code{HiCDOCDataSet} object.
#'
#' @return \code{HiCDOCDataSet} constructor returns an \code{HiCDOCDataSet}
#'         object of class S4.
#' @keywords internal
#' @noRd
.fillHiCDOCDataSet <- function(object) {

    # Fill binSize slot
    if (is.null(object@binSize)) object@binSize <- .computeBinSize(object)

    if (!is.null(object@interactions)) {
        chromosomes <-
            gtools::mixedsort(
                unique(as.character(object@interactions$chromosome))
            )
        object@chromosomes <- chromosomes

        # Determine positions
        if (is.null(object@positions)) {
            object@positions <- .determinePositions(object)
        }

        # Replace positions by bins
        object@interactions <- .replacePositionsByBins(object)

        object@weakBins <- vector("list", length(chromosomes))
        names(object@weakBins) <- chromosomes

        object@sparseConditions <- vector("list", length(chromosomes))
        names(object@sparseConditions) <- chromosomes

        object@sparseReplicates <- vector("list", length(chromosomes))
        names(object@sparseReplicates) <- chromosomes

        # Make upper triangular matrix, values as numeric, average duplicates
        object@interactions <- .reformatInteractions(object@interactions)
    }

    # Fill totalBins slot
    if (!is.null(object@binSize) & !is.null(object@interactions)) {
        object@totalBins <-
            vapply(
                as.character(object@chromosomes),
                function(x) {
                    nrow(object@positions[object@positions$chromosome == x, ])
                },
                FUN.VALUE = 0
            )
    }

    object@parameters <- defaultHiCDOCParameters
    return(object)
}
