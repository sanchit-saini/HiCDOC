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
#' object <- HiCDOCDataSetExample()
#' object <- HiCDOC(object)
#' plotCompartments(object, chromosome = 1)
#'
#' @export
plotCompartments <- function(
    object,
    chromosome,
    xlim = NULL
) {

    .validateSlots(object, slots = c("compartments", "positions"))
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")
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
        dplyr::mutate(position = start + 0.5 * object@resolution)

    plot <-
        ggplot(
            data = compartments,
            aes(x = position, fill = compartment)
        ) +
        geom_histogram(
            binwidth = object@resolution,
            colour = "gray90",
            size = 0.05
        ) +
        xlim(xlim[1] - 0.5 * object@resolution,
             xlim[2] + 0.5 * object@resolution) +
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
