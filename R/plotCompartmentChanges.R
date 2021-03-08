#' Run plotCompartments() and plotConcordances() and assemble them on the
#' same plot
#'
#' @param object An HiCDOCDataSet object, after a detectCompartments() run
#' @param chromosomeId Name or number of the chromosome, like in
#' object@chromosome
#' @param threshold Significance threshold for the changes. Default to 0.05.
#' @param xlim A numeric-value pair, indicating the interval of positions
#' to represent.
#' Default to NULL = all positions.
#' @param points Logical (default to FALSE). If TRUE, points will be added
#' on the concordance lines.
#'
#' @return A ggplot object.
#'
#' @examples
#' object <- HiCDOCDataSetExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotCompartmentChanges(object, 1)
#' @export
plotCompartmentChanges <- function(
    object,
    chromosome,
    threshold = 0.05,
    xlim = NULL,
    points = FALSE
) {

    .validateSlots(
        object,
        slots = c(
            "concordances",
            "compartments",
            "differences"
        )
    )
    chromosomeName <- .validateNameOrId(object, chromosome, "chromosomes")

    concordancesPlot <- plotConcordances(
        object,
        chromosomeName,
        xlim,
        threshold,
        points
    )

    compartmentsPlot <- plotCompartments(
        object,
        chromosomeName,
        xlim
    )

    # Horizontal alignment of the sub-graphs (change width of the plots)
    plotsGrobs <-
        lapply(
            list(
                compartmentsPlot + theme(legend.position = "none"),
                concordancesPlot + theme(legend.position = "none")
            ),
            ggplot2::ggplotGrob
        )

    commonWidths <- plotsGrobs[[length(plotsGrobs)]]$widths
    plotsGrobs <-
        lapply(
            plotsGrobs,
            function(x) {
                x$widths <- commonWidths
                return(x)
            }
        )

    plot <-
        ggpubr::as_ggplot(
            gridExtra::arrangeGrob(
                gridExtra::arrangeGrob(
                    plotsGrobs[[1]],
                    plotsGrobs[[2]],
                    heights = c(1, 5),
                    padding = unit(0, "cm")
                ),
                gridExtra::arrangeGrob(
                    .extractLegends(ggplotGrob(compartmentsPlot)),
                    .extractLegends(ggplotGrob(concordancesPlot)),
                    ncol = 2
                ),
                heights = c(10, 1),
                top = paste0("Chromosome ", chromosomeName)
            )
        )
    return(plot)
}
