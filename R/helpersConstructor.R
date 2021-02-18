# Helpers for HiCDOCDataSet constructor
# For internal use only
# These functions are called by HiCDOCDataSet() function

#### computeBinSize ##########################################################
#' Compute the binSize 
#'
#' The function takes a HiCDOCDataset object and compute the binSize from
#' the positions in the interactions matrix.
#' It takes the minimum of difference between position.1 and position.2
#' @param object a HiCDOCDataSet object.
#'
#' @return binSize, an integer.
#' @keywords internal
#' @noRd
computeBinSize <- function(object) {
    if (!is.null(object@interactions)) {
        binSize <- min(abs(
            object@interactions$position.2[object@interactions$position.1
            != object@interactions$position.2]
            - object@interactions$position.1[object@interactions$position.1
                != object@interactions$position.2]
        ))
        return(binSize)
    } else {
        return(NULL)
    }
}

#### determinePositions ######################################################
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
determinePositions <- function(object) {
    testSlotsHiCDOC(object, c("chromosomes", "interactions"))
    chromosomes <- object@chromosomes
    # Determine positions
    minMaxPos <- object@interactions %>%
        dplyr::group_by(chromosome) %>%
        dplyr::summarize(
            maxpos.1 = max(position.1),
            maxpos.2 = max(position.2)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            maxpos = max(maxpos.1, maxpos.2),
            .keep = "unused"
        )

    positions <- lapply(chromosomes, function(x) {
          return(
              data.frame(
                  "chromosome" = x,
                  "start" = seq(0,
                      minMaxPos[minMaxPos$chromosome == x, ]$maxpos,
                      by = object@binSize
                  )
              )
          )
      }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(
            chromosome =
                factor(chromosome, levels = chromosomes)
        ) %>%
        dplyr::group_by(chromosome) %>%
        dplyr::arrange(chromosome, start, .by_group = TRUE) %>%
        dplyr::mutate(end = dplyr::lead(start) - 1) %>%
        dplyr::mutate(bin = dplyr::row_number()) %>%
        dplyr::ungroup()

    positions %<>%
        dplyr::mutate(end = ifelse(is.na(end),
            start + object@binSize - 1, end
        )) %>%
        dplyr::select(chromosome, bin, start, end)
    return(positions)
}

#### replacePositionsByBins ##################################################
#' Replace the positions by the corresponding bins in the interaction matrix
#'
#' @param object a HiCDOCDataSet object
#'
#' @return a tibble, of the modified interactions matrix.
#' @keywords internal
#' @noRd
replacePositionsByBins <- function(object) {
    chromosomes <- object@chromosomes
    interactions <- object@interactions %>%
        dplyr::mutate(
            chromosome =
                factor(chromosome, levels = chromosomes)
        ) %>%
        dplyr::left_join(object@positions %>%
            dplyr::select(chromosome,
                position.1 = start,
                bin.1 = bin
            ),
        by = c("chromosome", "position.1")
        ) %>%
        dplyr::left_join(object@positions %>%
            dplyr::select(chromosome,
                position.2 = start,
                bin.2 = bin
            ),
        by = c("chromosome", "position.2")
        ) %>%
        dplyr::select(chromosome, bin.1, bin.2, condition, replicate, value)
    return(interactions)
}

#### reformatInteractions ##################################################
#' Reformat the interaction matrix
#' 
#' 
#' The function takes an interaction matrix and performs some correction
#' - the matrix is putted in upper format only (the values in the lower side
#' will be reaffected in the upper side)
#' - the value column is transformed in numeric
#' - the duplicated combinations of positions are aggregated by mean
#'
#' @param interactions the slot interactions of a HiCDOCDataSet object.
#'
#' @return a tibble, of the interactions corrected
#' @keywords internal
#' @noRd
reformatInteractions <- function(interactions) {
    # Put matrix in upper only
    low <- which(interactions$bin.1 > interactions$bin.2)
    if (length(low) > 0) {
        lowpos1 <- interactions[low, "bin.1"] %>% pull()
        lowpos2 <- interactions[low, "bin.2"] %>% pull()
        interactions[low, ]$bin.1 <- lowpos2
        interactions[low, ]$bin.2 <- lowpos1
    }
    # Transform value in numeric (avoid arrays, integer -> uniform)
    if (!is.numeric(interactions$value)) {
        interactions$value <- as.numeric(interactions$value)
    }
    # Check if unique combination of positions, if not mean(value)
    if (nrow(unique(interactions[, c(
        "chromosome",
        "bin.1",
        "bin.2",
        "condition",
        "replicate"
    )])) !=
        nrow(interactions)) {
        message("Not unique positions, replacing value by mean(value)")
        interactions %<>%
            dplyr::group_by(chromosome, bin.1, bin.2, condition, replicate) %>%
            dplyr::mutate(value = mean(value)) %>%
            dplyr::ungroup()
    }
    return(interactions)
}
