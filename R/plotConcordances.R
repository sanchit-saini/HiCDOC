#' Plot the concordance
#'
#' Plot the concordance of all the replicates for one chromosome.
#'
#' @param object A \code{HiCDOCDataSet} object on which
#' \code{detectCompartments()} has run.
#' @param chromosomeId The name or number of the chromosome to plot.
#' If number, will be taken in \code{object@chromosomes[chromosomeId]}
#' @param xlim A numeric-value pair, indicating the interval of positions
#' to represent.
#' Default to NULL = all positions.
#' @param threshold Significance threshold for the changes. Default to 0.05.
#' @param points Logical (default to FALSE). If TRUE, points will be added
#' on the concordance lines.
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
#' plotConcordances(object, 1)
#' @export
plotConcordances <- function(
    object,
    chromosome,
    xlim = NULL,
    threshold = 0.05,
    points = FALSE
) {

    .validateSlots(object, slots = c("concordances", "differences"))
    chromosomeName <- .validateNameOrId(object, chromosome, "chromosomes")
    xlim <- .validateXlim(xlim, object, chromosomeName)

    concordance <-
        object@concordances %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::left_join(
            object@positions %>% dplyr::filter(chromosome == chromosomeName),
            by = c("chromosome", "bin")
        ) %>%
        dplyr::filter(start >= xlim[1] & end <= xlim[2]) %>%
        dplyr::mutate(condition = paste0("confidence, cond. ", condition)) %>%
        dplyr::mutate(position = start + 0.5 * object@binSize)

    # Significant differences
    differences <-
        object@differences %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::left_join(
            object@positions %>% dplyr::filter(chromosome == chromosomeName),
            by = c("chromosome", "bin")
        ) %>%
        dplyr::filter(start >= xlim[1] & end <= xlim[2]) %>%
        dplyr::filter(pvalue.adjusted < threshold) %>%
        tidyr::pivot_longer(
            cols = tidyr::starts_with("condition"),
            values_to = "condition"
        ) %>%
        dplyr::mutate(condition = paste0("confidence, cond. ", condition))

    caption <- "The grey areas are significant changes"
    if (nrow(differences) == 0) caption <- "No change is significant"
    caption <- paste0(
        caption,
        " (adjusted p-value < ",
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
        xlim(xlim[1], xlim[2] + object@binSize) +
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
