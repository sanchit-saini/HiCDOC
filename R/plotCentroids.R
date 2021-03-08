#' Plot the centroids PCA for a given chromosome.
#'
#' @param object an \code{HiCDOCDataSet} object
#' @param chromosomeId Character or numeric value. Name or number of
#' the chromosome, like in chromosomes(object)
#' @param size Numeric. Size of the points in geom_point.
#'
#' @return A \code{ggplot} object
#' @examples
#' object <- HiCDOCDataSetExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotCentroids(object, 1)
#' @export
plotCentroids <- function(object, chromosome, size = 2) {
    .validateSlots(object, slots = c("centroids"))
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")

    df <-
        object@centroids %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::select(-chromosome) %>%
        tidyr::unite(name, c(condition, compartment))
    names <- df$name
    df %<>%
        tidyr::spread(name, centroid) %>%
        tidyr::unnest(cols = names) %>%
        t()

    pca <- stats::prcomp(df)
    varpca <- pca$sdev^2
    propvar <- varpca / sum(varpca)
    propvar <- paste(round(100 * propvar, 2), "%")

    pca <- as.data.frame(pca$x)
    pca$group <- row.names(df)
    plot <-
        ggplot(
            pca,
            aes(
                x = PC1,
                y = PC2,
                color = group,
                shape = group
            )
        ) +
        geom_point(size = size) +
        labs(
            title = paste0("Centroids of chromosome ", chromosomeName),
            x = paste("PC1 ", propvar[1]),
            y = paste("PC2 ", propvar[2])
        )
    return(plot)
}
