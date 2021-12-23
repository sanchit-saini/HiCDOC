#' @title
#' Plot A/B compartments.
#'
#' @description
#' Plots the predicted compartments in each experiment condition.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosome
#' A chromosome name or index in \code{chromosomes(object)}.
#' @param xlim
#' A vector of the minimum and maximum positions to display. If NULL, displays
#' all positions. Defaults to NULL.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' data(exampleHiCDOCDataSetProcessed)
#' plotCompartments(exampleHiCDOCDataSetProcessed, chromosome = 1)
#'
#' @export
plotCompartments <- function(
    object,
    chromosome,
    xlim = NULL
) {

    .validateSlots(object, slots = c("compartments"))
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")
    # xlim <- .validateXlim(xlim, object, chromosomeName)
    
    compartments <- object@compartments[chromosome == chromosomeName]
    # TODO : resolve xlim in case of compartemnts (index)
    # Add index param in validxlim ?
    
    if (nrow(compartments) == 0) {
        message("No compartments for chromosome ", chromosomeName, ".")
        return(NULL)
    }
    
    plot <-
        ggplot(
            data = compartments,
            aes(x = index, fill = compartment)
        ) +
        geom_histogram(
            binwidth = 1,
            colour = "gray90",
            size = 0.05
        ) +
        # xlim(xlim[1] - 0.5 * object@binSize, xlim[2] + 0.5 * object@binSize) +
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
