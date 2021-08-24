#' @title
#' Plot centroids.
#'
#' @description
#' Plots the result of the PCA on the compartments' centroids.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosome
#' A chromosome name or index in \code{chromosomes(object)}.
#' @param size
#' Size of each point. Defaults to 2.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' data(exampleHiCDOCDataSetProcessed)
#' plotCentroids(exampleHiCDOCDataSetProcessed, chromosome = 1)
#'
#' @export
plotCentroids <- function(object, chromosome, size = 2) {
    .validateSlots(object, slots = "centroids")
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")

    df <-
        object@centroids %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::select(-chromosome) %>%
        tidyr::unite(name, c(condition, compartment))

    if (nrow(df) == 0) {
        message("No centroids for chromosome ", chromosomeName, ".")
        return(NULL)
    }

    names <- df$name
    df %<>%
        tidyr::spread(name, centroid) %>%
        tidyr::unnest(cols = dplyr::all_of(names)) %>%
        t()

    pca <- stats::prcomp(df)
    varpca <- pca$sdev^2
    propvar <- varpca / sum(varpca)
    propvar <- paste(round(100 * propvar, 2), "%")

    pca <- as.data.frame(pca$x)
    pca$group <- row.names(df)
    pca %<>%
        dplyr::mutate(group = strsplit(group, "_", fixed = TRUE)) %>%
        dplyr::mutate(condition = purrr::map_chr(group, 1),
                      compartment = purrr::map_chr(group, 2))

    plot <-
        ggplot(
            pca,
            aes(
                x = PC1,
                y = PC2,
                color = compartment,
                shape = condition
            )
        ) +
        geom_point(size = size) +
        labs(
            title = paste0("PCA on centroids of chromosome ", chromosomeName),
            x = paste("PC1 ", propvar[1]),
            y = paste("PC2 ", propvar[2])
        )
    return(plot)
}
