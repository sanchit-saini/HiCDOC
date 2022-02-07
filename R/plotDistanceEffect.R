#' @title
#' Plot the distance effect.
#'
#' @description
#' Plots the distance effect on proportion of interactions.
#' Each point is a cell in the interaction matrix, such that
#' the x-axis is the distance with respect to the diagonal,
#' the y-axis is number of counts.
#' Dots are binned.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosome
#' Name (character) or index of the chromosome, if the plot should be
#' restricted to only one chromosome. Default to NULL.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' data(exampleHiCDOCDataSet)
#' plotDistanceEffect(exampleHiCDOCDataSet)
#'
#' @export
plotDistanceEffect <- function(object, chromosome = NULL) {
    .validateSlots(object, slots = c("interactions"))
    if (!is.null(chromosome)) {
        if (length(chromosome) > 1) {
            warning("`chromosome` should be of length 1, ",
                    "taking the first one")
            chromosome < chromosome[1]
        }
        chromosomeName <-
            .validateNames(object, chromosome, "chromosomes")
        rowsId <-
            as.logical(S4Vectors::mcols(object)$Chr == chromosomeName)
        addTitle <- paste(", chromosome", chromosomeName)
    } else {
        rowsId <- rep(TRUE, length(object))
        addTitle <- ""
    }
    
    distances <-
        InteractionSet::pairdist(object, type = "mid")[rowsId]
    matAssay <- SummarizedExperiment::assay(object)[rowsId, ]
    dfDistance <-
        data.table::data.table("distance" = rep(distances, ncol(matAssay)),
                               "interaction" = as.vector(matAssay))
    dfDistance <- dfDistance[!is.na(interaction)]
    
    plot <-
        ggplot(dfDistance, aes(x = distance, y = interaction)) +
        geom_bin2d() +
        scale_fill_gradient(low = "white",
                            high = "blue",
                            trans = "log2") +
        geom_point(col = "transparent") + # necessary for geom_smooth
        geom_smooth(col = "red") +
        labs(title = paste0("Distance effect", addTitle))
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
