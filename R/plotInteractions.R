#' @title
#' Plot interaction matrices.
#'
#' @description
#' Plots the interaction matrices as heatmaps.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosome
#' A chromosome name or index in \code{chromosomes(object)}.
#' @param transform
#' Transformation of the color scale. Set to NULL for no transformation. See
#' \code{\link[ggplot2]{scale_fill_gradientn}} for other accepted values.
#' Defaults to "log2".
#' @param colours
#' A character vector of colours to use for the gradient. See
#' \code{\link[ggplot2]{scale_fill_gradientn}} for more info. Defaults to
#' \code{c("#000066", "#ffffbf", "#990000")}.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' object <- HiCDOCDataSetExample()
#' p <- plotInteractions(object, chromosome = 1)
#'
#' @export
plotInteractions <- function(
    object,
    chromosome,
    transform = "log2",
    colours = c("#000066", "#ffffbf", "#990000")
) {

    .validateSlots(
        object,
        slots = c(
            "interactions",
            "conditions",
            "totalBins",
            "binSize",
            "positions"
        )
    )
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")
    if (is.null(transform)) transform <- "Identity"

    positions <-
        object@positions %>%
        dplyr::filter(chromosome == chromosomeName)
    interactions <-
        object@interactions %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::filter(interaction > 0) %>%
        dplyr::left_join(
            positions %>% dplyr::select(
                bin.1 = bin,
                position.1 = start
            ),
            by = "bin.1"
        ) %>%
        dplyr::left_join(
            positions %>% dplyr::select(
                bin.2 = bin,
                position.2 = start
            ),
            by = "bin.2"
        )

    if (nrow(interactions) == 0) {
        message("No interaction data to plot.")
        return(NULL)
    }

    totalLevels <- table(object@conditions)
    totalRows <- 1
    if (max(totalLevels) == min(totalLevels)) {
        totalRows <- length(unique(object@conditions))
    }
    xylim <- c(min(positions$start), max(positions$start))

    plot <-
        ggplot(
            data = interactions,
            aes(x = position.1, y = position.2, z = interaction)
        ) +
        geom_raster(aes(fill = interaction), na.rm = TRUE) +
        geom_raster(
            data = interactions[interactions$bin.1 != interactions$bin.2, ],
            aes(x = position.2, y = position.1, fill = interaction),
            na.rm = TRUE
        ) +
        coord_fixed(ratio = 1) +
        theme_bw() +
        xlim(xylim) +
        scale_y_reverse(limits = rev(xylim)) +
        facet_wrap(
            condition ~ replicate,
            nrow = totalRows,
            labeller = label_wrap_gen(multi_line = FALSE)
        ) +
        labs(title = paste("Chromosome ", chromosomeName), x = "", y = "") +
        scale_fill_gradientn(
            colours = colours,
            trans = transform,
            name = "Intensity",
            na.value = "transparent"
        )
    return(plot)
}
