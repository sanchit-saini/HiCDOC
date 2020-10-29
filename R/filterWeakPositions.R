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
#' @param chromosomeId A character or numeric value.
#' Name or number of the chromosome
#' @param threshold Numeric default to 0.
#'
#' @return list of length 2 : \code{"pos"} = the weak positions,
#' \code{"interactions"} the interactions matrix for the chromosome,
#' whithout the weak bins.
filterWeakChr <- function(object, chromosomeId, threshold = 0) {
    testSlotsHiCDOC(
        object,
        slots = c(
            "interactions",
            "binSize",
            "replicates",
            "conditions",
            "totalBins"
        )
    )
    chr <- testChromosome(object, chromosomeId)

    # Initialization
    interChr <- object@interactions %>%
        dplyr::filter(chromosome == chr) %>%
        dplyr::filter(value > 0)

    totalBinsChr <- object@totalBins[[chr]]
    fullPos <- (seq_len(totalBinsChr) - 1) * object@binSize
    removedBins <-
        fullPos[!(fullPos %in%
                      unique(c(interChr$position.1, interChr$position.2))
               )]

    nbNewEmptyBins <- 1
    nbRemovedBins <- 0

    # Recursive removal of bins - deleting a bin can create a new weak bin.
    while (nbNewEmptyBins > 0 & nbRemovedBins <= totalBinsChr) {
        interChrDiag <- interChr %>%
          dplyr::filter(position.1 == position.2) %>%
          dplyr::rename(position = position.1) %>%
          dplyr::select(-chromosome, -position.2)

        existing <- interChr %>%
            dplyr::filter(position.1 != position.2) %>%
            tidyr::pivot_longer(
                cols = c(`position.1`, `position.2`),
                names_to = "pos_1or2",
                names_prefix = "position.",
                values_to = "position"
            ) %>%
            dplyr::select(-pos_1or2, -chromosome) %>%
          bind_rows(interChrDiag) %>%
            tidyr::complete(position,
                     nesting(condition, replicate),
                     fill = list(value = 0)
            )
        weakdataset <- existing %>%
            dplyr::group_by(replicate, condition, position) %>%
            dplyr::mutate(mean = sum(value) / totalBinsChr) %>%
            dplyr::filter(mean <= threshold)

        weakposChr <- weakdataset %>%
            dplyr::pull(position) %>%
            unique() %>%
            sort()

        nbNewEmptyBins <- length(weakposChr) - nbRemovedBins
        removedBins <- c(removedBins, weakposChr)

        # Remove the positions in rows and columns in empty bins
        if (nbNewEmptyBins > 0) {
          interChr <- interChr[!(interChr$position.1 %in% weakposChr),]
          interChr <- interChr[!(interChr$position.2 %in% weakposChr),]
          nbRemovedBins <- nbRemovedBins + nbNewEmptyBins
        }
    }

    message(paste0(
        "Chromosome ",
        chr,
        ": ",
        length(removedBins),
        " position(s) removed"
    ))
    return(list("pos" = removedBins, "interactions" = interChr))
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
#' @param threshold Numeric, default to HiCDOCDefaultParameters$weakPosThreshold
#'
#' @return A \code{HiCDOCDataSet} object with a reduced \code{interactions} slot,
#' and the weak bins identified by chromosomes in \code{object@weakBins}.
#' @export
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- filterWeakPositions(exp)
filterWeakPositions <- function(object, threshold = NULL) {
    if(! is.null(threshold) ) object@parameters$weakPosThreshold <- threshold
    object@parameters <- checkParameters(object@parameters,
                                       c("weakPosThreshold"))
    weakPositions <- lapply(object@chromosomes,
        function(chr) filterWeakChr(object,
                                    chr,
                                    object@parameters$weakPosThreshold)
        )
    names(weakPositions) <- object@chromosomes

    weakBins <- weakPositions %>% purrr::map("pos")
    weakBins <- lapply(weakBins, function(x)
        x / object@binSize + 1)
    nullweak <-
        vapply(weakBins, function(x)
            length(x) == 0, FUN.VALUE = TRUE)
    weakBins[nullweak] <- list(NULL)

    interactions <- weakPositions %>% purrr::map_dfr("interactions")

    # Save new values
    object@weakBins <- weakBins
    object@interactions <- interactions

    nbweak <- sum(vapply(weakBins, length, c(0)))

    message("Removed ",
            nbweak,
            " position",
            if (nbweak != 1) "s"
    )

    if (nbweak >= sum(unlist(object@totalBins))) {
        message('No data left!')
    }

  return (object)
}
