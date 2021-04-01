#' @title
#' Plot concordances.
#'
#' @description
#' Plots the concordances of each replicate in each experiment condition. A
#' concordance can be understood as a confidence in a genomic position's
#' assigned compartment. Mathematically, it is the log ratio of a genomic
#' position's distance to each compartment's centroid, normalized by the
#' distance between both centroids, and min-maxed to a [-1,1] interval.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosome
#' A chromosome name or index in \code{chromosomes(object)}.
#' @param xlim
#' A vector of the minimum and maximum positions to display. If NULL, displays
#' all positions. Defaults to NULL.
#' @param threshold
#' Significance threshold for the compartment changes. Defaults to 0.05.
#' @param points
#' Whether or not to add points to the concordances. Defaults to FALSE.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' data(HiCDOCDataSetExample)
#' object <- HiCDOC(HiCDOCDataSetExample)
#' plotConcordances(object, chromosome = 1)
#'
#' @export
plotConcordances <- function(
    object,
    chromosome,
    xlim = NULL,
    threshold = 0.05,
    points = FALSE
) {

    .validateSlots(object, slots = c("concordances", "differences"))
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")
    xlim <- .validateXlim(xlim, object, chromosomeName)

    concordance <-
        object@concordances %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::left_join(
            object@positions %>% dplyr::filter(chromosome == chromosomeName),
            by = c("chromosome", "bin")
        ) %>%
        dplyr::filter(start >= xlim[1] & end <= xlim[2]) %>%
        dplyr::mutate(condition = paste0("Concordances\n", condition)) %>%
        dplyr::mutate(position = start + 0.5 * object@resolution)

    # Significant differences
    differences <-
        object@differences %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::left_join(
            object@positions %>% dplyr::filter(chromosome == chromosomeName),
            by = c("chromosome", "bin")
        ) %>%
        dplyr::filter(start >= xlim[1] & end <= xlim[2]) %>%
        dplyr::filter(pvalue.adjusted <= threshold) %>%
        tidyr::pivot_longer(
            cols = tidyr::starts_with("condition"),
            values_to = "condition"
        ) %>%
        dplyr::mutate(condition = paste0("Concordances\n", condition))

    caption <- "The grey areas are significant changes"
    if (nrow(differences) == 0) caption <- "No change is significant"
    caption <- paste0(
        caption,
        " (adjusted p-value <= ",
        round(100 * threshold, 2),
        "%)"
    )

    ylim <-
        c(
            min(concordance$concordance, na.rm = TRUE),
            max(concordance$concordance, na.rm = TRUE)
        )

    plot <- ggplot()
    if (nrow(differences) > 0) {
        plot <-
            plot +
            geom_rect(
                data = differences,
                aes(
                    xmin = start,
                    xmax = end,
                    ymin = ylim[1],
                    ymax = ylim[2]
                ),
                color = NA,
                fill = "gray80"
            )
    }
    plot <-
        plot +
        geom_line(
            data = concordance,
            aes(x = position, y = concordance, color = replicate)
        )
    if (points) {
        plot <-
            plot +
            geom_point(
                data = concordance,
                aes(
                    x = position,
                    y = concordance,
                    color = replicate
                )
            )
    }
    plot <-
        plot +
        labs(caption = caption) +
        xlim(xlim[1], xlim[2] + object@resolution) +
        ylim(ylim) +
        geom_hline(yintercept = 0.0, size = 0.1) +
        facet_grid(rows = vars(condition), margins = FALSE, switch = "y") +
        theme_minimal() +
        theme(
            axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "bottom",
            strip.placement = "outside"
        )
    return(plot)
}
