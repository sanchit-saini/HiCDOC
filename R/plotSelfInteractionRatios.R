#' Plot boxplots of the diagonal ratios in each compartment,
#' for a given chromosome.
#'
#' @param object an \code{HiCDOCDataSet} object
#' @param chromosome character or numeric value, name or number of chromosome
#' @return A list of \code{ggplot}, one for each chromosome.
#' @examples
#' object <- HiCDOCDataSetExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotSelfInteractionRatios(object, 1)
#' @export
plotSelfInteractionRatios <- function(object, chromosome) {
    .validateSlots(object, slots = c("selfInteractionRatios", "compartments"))
    chromosomeName <- .validateNameOrId(object, chromosome, "chromosomes")

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
                "and other interactions in chromosome ",
                chromosomeName
            )
        )
    return(plot)
}
