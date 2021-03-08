#' Plot the A & B compartments
#'
#' Plot the A and B compartments after \code{detectCompartments()} on
#' a HiCDOCDataSet object, for a chromosome.
#'
#' @param object A \code{HiCDOCDataSet} object on which
#' \code{detectCompartments()} has run.
#' @param chromosomeId The name or number of the chromosome to plot.
#' If number, will be taken in \code{object@chromosomes[chromosomeId]}
#' @param xlim A numeric-value pair, indicating the interval of positions
#' to represent.
#' Default to NULL = all positions.
#'
#' @return A ggplot.
#' @export
#'
#' @examples
#' object <- HiCDOCDataSetExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotCompartments(object, 1)
#' @export
plotCompartments <- function(
    object,
    chromosome,
    xlim = NULL
) {

    .validateSlots(object, slots = c("compartments", "positions"))
    chromosomeName <- .validateNameOrId(object, chromosome, "chromosomes")
    xlim <- .validateXlim(xlim, object, chromosomeName)

    compartments <-
        object@compartments %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::left_join(
            object@positions %>%
            dplyr::filter(chromosome == chromosomeName),
            by = c("chromosome", "bin")
        ) %>%
        dplyr::filter(start >= xlim[1] & end <= xlim[2]) %>%
        dplyr::mutate(compartment = factor(compartment)) %>%
        dplyr::mutate(position = start + 0.5 * object@binSize)

    plot <-
        ggplot(
            data = compartments,
            aes(x = position, fill = compartment)
        ) +
        geom_histogram(
            binwidth = object@binSize,
            colour = "gray90",
            size = 0.05
        ) +
        xlim(xlim[1] - 0.5 * object@binSize, xlim[2] + 0.5 * object@binSize) +
        facet_grid(rows = vars(condition), margins = FALSE, switch = "y") +
        theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            legend.position = "bottom",
            strip.placement = "outside"
        )
    return(plot)
}
