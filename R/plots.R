#' Plot the interaction matrix (as heatmap).
#'
#' @param object an \code{HiCDOCDataSet} object
#' @param chromosomeId The name or number of the chromosome to plot.
#' If number, will be taken in \code{object@chromosomes[chromosomeId]}
#' @param trans character: transformation of the color scale. Default to "log2".
#' See \code{\link[ggplot2]{scale_fill_gradientn}} for other accepted values.
#' Set to NULL for no transformation.
#' @param colours character vector, vector of colours to use for n-colour
#' gradient.
#' Default to \code{c("#000066", "#ffffbf", "#990000")}.
#'
#' @return A \code{ggplot} object.
#' @examples
#' object <- HiCDOCExample()
#' p <- plotInteractionMatrix(object, chromosomeId = 1, trans = "log2")
#' @export
plotInteractionMatrix <-
    function(object,
    chromosomeId,
    trans = "log2",
    colours = c("#000066", "#ffffbf", "#990000")) {
        # Parameters
        testSlotsHiCDOC(object,
            slots = c(
                "interactions",
                "conditions",
                "totalBins",
                "binSize",
                "positions"
            )
        )
        chr <- testValidId(object, chromosomeId, "chromosomes")
        if (is.null(trans)) {
              trans <- "Identity"
          }

        posChr <- object@positions %>%
            dplyr::filter(chromosome == chr)
        # Prepare data
        interactionsChr <- object@interactions %>%
            dplyr::filter(chromosome == chr & value > 0) %>%
            dplyr::left_join(posChr %>%
                dplyr::select(
                    bin.1 = bin,
                    position.1 = start
                ),
            by = "bin.1"
            ) %>%
            dplyr::left_join(posChr %>%
                dplyr::select(
                    bin.2 = bin,
                    position.2 = start
                ),
            by = "bin.2"
            )
        nblevels <- table(object@conditions)
        nbrows <- 1
        if (max(nblevels) == min(nblevels)) {
              nbrows <- length(unique(object@conditions))
          }
        xylim <- c(min(posChr$start), max(posChr$start))

        if (nrow(interactionsChr) > 0) {
            p <-
                ggplot(
                    data = interactionsChr,
                    aes(x = position.1, y = position.2, z = value)
                ) +
                geom_raster(aes(fill = value), na.rm = TRUE) +
                geom_raster(
                    data = interactionsChr[interactionsChr$bin.1 !=
                        interactionsChr$bin.2, ],
                    aes(x = position.2, y = position.1, fill = value),
                    na.rm = TRUE
                ) +
                coord_fixed(ratio = 1) +
                theme_bw() +
                xlim(xylim) +
                scale_y_reverse(limits = rev(xylim)) +
                facet_wrap(condition ~ replicate,
                    nrow = nbrows,
                    labeller = label_wrap_gen(multi_line = FALSE)
                ) +
                labs(title = paste("Chromosome:", chr), x = "", y = "")
            p <-
                p + scale_fill_gradientn(
                    colours = colours,
                    trans = trans,
                    name = "Intensity",
                    na.value = "transparent"
                )
        } else {
            message("No interaction data, with positive value to plot")
            p <- NULL
        }
        return(p)
    }

#' Plot the distance vs intensity matrix.
#'
#' @param object an \code{HiCDOCDataSet} object
#' @return A \code{ggplot}.
#' @examples
#' object <- HiCDOCExample()
#' plotDistanceEffect(object)
#' @export
plotDistanceEffect <- function(object) {
    testSlotsHiCDOC(object, slots = c("interactions", "binSize"))

    dataplot <- object@interactions %>%
        dplyr::mutate(distance = (bin.2 - bin.1) * object@binSize)

    p <- ggplot(dataplot, aes(x = distance, y = value)) +
        geom_bin2d() +
        scale_fill_gradient(
            low = "white",
            high = "blue",
            trans = "log2"
        ) +
        geom_point(col = "transparent") + # necessary for geom_smooth
        geom_smooth(col = "red") +
        labs(title = "Distance effect")
    p <- ggExtra::ggMarginal(p,
        margins = "x",
        type = "histogram",
        fill = "transparent",
        lwd = 0.5
    )
    return(p)
}

#' Plot the concordance, i.e. the relative distance of the genomic positions
#' with respect to the centroids.
#'
#' @param object an \code{HiCDOCDataSet} object
#' @return A list of \code{ggplot}, one for each chromosome.
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotDiffConcordances(object)
#' @export
plotDiffConcordances <- function(object) {
    testSlotsHiCDOC(object,
        slots = c(
            "interactions",
            "differences",
            "concordances"
        )
    )

    changed <- object@differences %>%
        dplyr::select(-c(`padj`)) %>%
        dplyr::mutate(changed = "T")

    differences <- object@concordances %>%
        dplyr::group_by(.dots = c("chromosome", "bin", "condition")) %>%
        dplyr::summarise(median = stats::median(concordance)) %>%
        dplyr::ungroup() %>%
        tidyr::spread(condition, median) %>%
        dplyr::mutate(difference = `2` - `1`) %>%
        dplyr::select(-c(`1`, `2`)) %>%
        dplyr::left_join(changed, by = c("chromosome", "bin")) %>%
        dplyr::mutate(changed = tidyr::replace_na(changed, "F"))

    p <- ggplot(differences, aes(x = difference, fill = changed)) +
        geom_histogram() +
        labs(
            x = "Concordance",
            title = "Distribution of the differences of concordances"
        )
    return(p)
}

#' Plot the distribution of A/B compartments along the genomic positions.
#'
#' @param object an \code{HiCDOCDataSet} object
#' @param chromosomeId character or numeric value, name or number of chromosome
#' @return A list of \code{ggplot}, one for each chromosome.
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotAB(object, 1)
#' @export
plotAB <- function(object, chromosomeId) {
    testSlotsHiCDOC(object, slots = c("diagonalRatios", "compartments"))
    chromosomeId <- testValidId(object, chromosomeId, "chromosomes")

    data <- object@diagonalRatios %>%
        dplyr::left_join(
            object@compartments,
            by = c("chromosome", "condition", "bin")
        ) %>%
        dplyr::filter(chromosome == chromosomeId)

    ggplot(data, aes(x = compartment, y = value)) +
        geom_jitter(aes(color = compartment)) +
        geom_boxplot(
            outlier.colour = NA,
            fill = NA,
            colour = "grey20"
        ) +
        labs(
            color = "Compartment",
            x = "Compartment",
            y = "Difference of int.",
            title = paste0("Chromosome: ", chromosomeId)
        )
}


#' Plot the centroid distributions along the genomic positions
#' for a given chromosome.
#'
#' @param object an \code{HiCDOCDataSet} object
#' @param chromosomeId Character or numeric value. Name or number of
#' the chromosome, like in chromosomes(object)
#' @param size Numeric. Size of the points in geom_point.
#'
#' @return A \code{ggplot} object
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotCentroids(object, 1)
#' @export
plotCentroids <- function(object, chromosomeId, size = 2) {
    testSlotsHiCDOC(object, slots = c("centroids"))
    chr <- testValidId(object, chromosomeId, "chromosomes")

    df <- object@centroids %>%
        dplyr::filter(chromosome == chr) %>%
        dplyr::select(-chromosome) %>%
        tidyr::unite(name, c(condition, compartment))
    names <- df$name
    df <- df %>%
        tidyr::spread(name, centroid) %>%
        tidyr::unnest(cols = names) %>%
        t()

    pca <- stats::prcomp(df)
    varpca <- pca$sdev^2
    propvar <- varpca / sum(varpca)
    propvar <- paste(round(100 * propvar, 2), "%")

    pca <- as.data.frame(pca$x)
    pca$group <- row.names(df)
    p <-
        ggplot(pca, aes(
            x = PC1,
            y = PC2,
            color = group,
            shape = group
        )) +
        geom_point(size = size) +
        labs(
            title = paste0("Centroids of chromosome ", chr),
            x = paste("PC1 ", propvar[1]),
            y = paste("PC2 ", propvar[2])
        )
    return(p)
}
