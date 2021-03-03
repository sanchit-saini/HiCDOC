#' Return valid xlim
#'
#' test if the xlim is valid (length(2)) or \code{NULL} and return
#' the xlim for the plot.
#' If \code{NULL} or invalid return the min and max of \code{positions}.
#'
#' @param xlim The entred xlim
#' @param positions The positions from the data
#'
#' @return A length 2 numerical vector
#' @keywords internal
validateXlim <- function(xlim, object, chromosomeName) {
    if (!is.null(xlim)) {
        if (length(xlim) != 2) {
            message(
                "Incorrect values for xlim (numerical of length 2 expected). ",
                "Setting xlim to NULL."
            )
            xlim <- NULL
        } else {
            xlim <- sort(xlim)
        }
    }
    if (is.null(xlim)) {
        positions <-
            object@positions %>%
            dplyr::filter(chromosome == chromosomeName) %>%
            dplyr::select(start) %>%
            dplyr::pull()
        xlim <- c(min(positions, na.rm = TRUE), max(positions, na.rm = TRUE))
    }
    return(xlim)
}


#' Function to extract legends from ggplot2 objects
#'
#' Inspired from gtable::gtable_filter (not used to avoid dependency)
#'
#' @param x a grob
#'
#' @return the legend as a grob
#' @keywords internal
extractLegends <- function(x) {
    matches <- grepl("guide-box", .subset2(x$layout, "name"), fixed = FALSE)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    return(x)
}

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
#' object <- HiCDOCExample()
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

    validateSlots(object, slots = c("concordances", "differences"))
    chromosomeName <- validateNameOrId(object, chromosome, "chromosomes")
    xlim <- validateXlim(xlim, object, chromosomeName)

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
#' object <- HiCDOCExample()
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

    validateSlots(object, slots = c("compartments", "positions"))
    chromosomeName <- validateNameOrId(object, chromosome, "chromosomes")
    xlim <- validateXlim(xlim, object, chromosomeName)

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
#' object <- HiCDOCExample()
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

    validateSlots(
        object,
        slots = c(
            "concordances",
            "compartments",
            "differences"
        )
    )
    chromosomeName <- validateNameOrId(object, chromosome, "chromosomes")

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
                    extractLegends(ggplotGrob(compartmentsPlot)),
                    extractLegends(ggplotGrob(concordancesPlot)),
                    ncol = 2
                ),
                heights = c(10, 1),
                top = paste0("Chromosome ", chromosomeName)
            )
        )
    return(plot)
}
