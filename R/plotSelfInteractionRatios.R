#' @title
#' Plot boxplots of self interaction ratios.
#'
#' @description
#' Plots the boxplots of self interaction ratios, which are the differences
#' between self interaction and median of other interactions for each genomic
#' position. Since the A compartment is open with more interactions overall, it
#' is assumed that self interaction ratios in compartment A are smaller than in
#' compartment B.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosome
#' A chromosome name or index in \code{chromosomes(object)}.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' data(HiCDOCDataSetExample)
#' object <- HiCDOC(HiCDOCDataSetExample)
#' plotSelfInteractionRatios(object, chromosome = 1)
#'
#' @export
plotSelfInteractionRatios <- function(object, chromosome) {
    .validateSlots(object, slots = c("selfInteractionRatios", "compartments"))
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")

    data <-
        object@selfInteractionRatios %>%
        dplyr::left_join(
            object@compartments,
            by = c("chromosome", "condition", "bin")
        ) %>%
        dplyr::filter(chromosome == chromosomeName)

    plot <-
        ggplot(data, aes(x = compartment, y = ratio)) +
        geom_jitter(aes(color = compartment)) +
        geom_boxplot(
            outlier.colour = NA,
            fill = NA,
            colour = "grey20"
        ) +
        labs(
            color = "Compartment",
            x = "Compartment",
            y = "Interaction difference",
            title = paste0(
                "Differences between self-interactions ",
                "and other interactions"
            ),
            subtitle = paste0("Chromosome ", chromosomeName)
        )
    return(plot)
}
