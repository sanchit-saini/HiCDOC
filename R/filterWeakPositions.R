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
    # Initialization
    interChr <- object@interactions %>%
        dplyr::filter(chromosome == chromosomeId) %>%
        dplyr::filter(value > 0)

    totalBinsChr <- object@totalBins[[chromosomeId]]
    fullPos <- seq_len(totalBinsChr)
    removedBins <-
        fullPos[!(fullPos %in%
            unique(c(interChr$bin.1, interChr$bin.2))
        )]

    nbNewEmptyBins <- 1
    nbRemovedBins <- 0

    # Recursive removal of bins - deleting a bin can create a new weak bin.
    while (nbNewEmptyBins > 0 & nbRemovedBins <= totalBinsChr) {
        interChrDiag <- interChr %>%
            dplyr::filter(bin.1 == bin.2) %>%
            dplyr::rename(bin = bin.1) %>%
            dplyr::select(-chromosome, -bin.2)

        existing <- interChr %>%
            dplyr::filter(bin.1 != bin.2) %>%
            tidyr::pivot_longer(
                cols = c(`bin.1`, `bin.2`),
                names_to = "pos_1or2",
                names_prefix = "bin.",
                values_to = "bin"
            ) %>%
            dplyr::select(-pos_1or2, -chromosome) %>%
            dplyr::bind_rows(interChrDiag) %>%
            tidyr::complete(bin,
                tidyr::nesting(condition, replicate),
                fill = list(value = 0)
            )
        weakdataset <- existing %>%
            dplyr::group_by(replicate, condition, bin) %>%
            dplyr::mutate(mean = sum(value) / totalBinsChr) %>%
            dplyr::filter(mean <= threshold)

        weakposChr <- weakdataset %>%
            dplyr::pull(bin) %>%
            unique() %>%
            sort()

        nbNewEmptyBins <- length(weakposChr) - nbRemovedBins
        removedBins <- c(removedBins, weakposChr)

        # Remove the positions in rows and columns in empty bins
        if (nbNewEmptyBins > 0) {
            interChr <- interChr[!(interChr$bin.1 %in% weakposChr), ]
            interChr <- interChr[!(interChr$bin.2 %in% weakposChr), ]
            nbRemovedBins <- nbRemovedBins + nbNewEmptyBins
        }
    }

    message(paste0(
        "Chromosome ",
        chromosomeId,
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
#' @param threshold Numeric. Threshold to consider bins as 'weaks'. If NULL,
#' default to the first not NULL of \code{object$weakPosThreshold} and
#' \code{HiCDOCDefaultParameters$weakPosThreshold}.
#'
#' @return A \code{HiCDOCDataSet} object with a reduced \code{interactions}
#' slot, and the weak bins identified by chromosomes in \code{object@weakBins}.
#' @export
#'
#' @seealso \code{\link[HiCDOC]{filterSmallChromosomes}}, 
#' \code{\link[HiCDOC]{filterSparseChromosomes}} and 
#' \code{\link[HiCDOC]{runHiCDOC}} for the recommended pipeline.
#' @examples
#' exp <- HiCDOCExample()
#' exp <- filterWeakPositions(exp)
filterWeakPositions <- function(object, threshold = NULL) {
    testSlotsHiCDOC(
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
        object@parameters$weakPosThreshold <- threshold
    }
    object@parameters <- checkParameters(object@parameters)
    weakPositions <- lapply(
        object@chromosomes,
        function(chr) {
            filterWeakChr(
                object,
                chr,
                object@parameters$weakPosThreshold
            )
        }
    )
    names(weakPositions) <- object@chromosomes

    weakBins <- weakPositions %>% purrr::map("pos")
    nullweak <-
        vapply(weakBins, function(x) {
              length(x) == 0
          }, FUN.VALUE = TRUE)
    weakBins[nullweak] <- list(NULL)
    interactions <- weakPositions %>% purrr::map_dfr("interactions")

    # Save new values
    object@weakBins <- weakBins
    object@interactions <- interactions
    nbweak <- sum(vapply(weakBins, length, c(0)))
    message(
        "Removed ",
        nbweak,
        " position",
        if (nbweak != 1) "s"
    )
    if (nbweak >= sum(unlist(object@totalBins))) {
        message("No data left!")
    }

    return(object)
}
