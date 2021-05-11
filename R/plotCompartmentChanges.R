#' @title
#' Plot compartment changes.
#'
#' @description
#' Plots the predicted compartments, along with their concordance in each
#' replicate, and significant changes between experiment conditions.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosome
#' A chromosome name or index in \code{chromosomes(object)}.
#' @param threshold
#' Significance threshold for the compartment changes. Defaults to 0.05.
#' @param xlim
#' A vector of the minimum and maximum positions to display. If NULL, displays
#' all positions. Defaults to NULL.
#' @param points
#' Whether or not to add points to the concordances. Defaults to FALSE.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' data(exampleHiCDOCDataSet)
#' object <- HiCDOC(exampleHiCDOCDataSet)
#' plotCompartmentChanges(object, chromosome = 1)
#'
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
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")

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

    if (is.null(compartmentsPlot) || is.null(concordancesPlot)) {
        return(NULL)
    }

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
                top = paste0(
                    "Compartments and concordances of chromosome ",
                    chromosomeName
                )
            )
        )
    return(plot)
}
