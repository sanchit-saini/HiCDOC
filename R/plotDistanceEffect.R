#' @title
#' Plot the distance effect.
#'
#' @description
#' Plots the distance effect on proportion of interactions.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' data(exampleHiCDOCDataSet)
#' plotDistanceEffect(exampleHiCDOCDataSet)
#'
#' @export
plotDistanceEffect <- function(object) {
    .validateSlots(object, slots = c("interactions"))
    distances <- InteractionSet::pairdist(object, type="mid")
    matAssay <- SummarizedExperiment::assay(object)
    dfDistance <- data.table("distance" = rep(distances, ncol(matAssay)),
                             "interaction" = as.vector(matAssay))
    dfDistance <- dfDistance[!is.na(interaction)]
    
    plot <-
        ggplot(dfDistance, aes(x = distance, y = interaction)) +
        geom_bin2d() +
        scale_fill_gradient(
            low = "white",
            high = "blue",
            trans = "log2"
        ) +
        geom_point(col = "transparent") + # necessary for geom_smooth
        geom_smooth(col = "red") +
        labs(title = "Distance effect")
    plot <-
        ggExtra::ggMarginal(
            plot,
            margins = "x",
            type = "histogram",
            fill = "transparent",
            lwd = 0.5
        )
    return(plot)
}
