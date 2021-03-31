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
#' object <- HiCDOCDataSetExample()
#' plotDistanceEffect(object)
#'
#' @export
plotDistanceEffect <- function(object) {
    .validateSlots(object, slots = c("interactions", "resolution"))
    data <-
        object@interactions %>%
        dplyr::mutate(distance = (bin.2 - bin.1) * object@resolution)
    plot <-
        ggplot(data, aes(x = distance, y = interaction)) +
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
