#' Plot the interaction matrix (as heatmap).
#'
#' @param object an \code{HiCDOCDataSet} object
#' @param chromosome The name or number of the chromosome to plot.
#' If number, will be taken in \code{object@chromosomes[chromosome]}
#' @param transform character: transformation of the color scale. Default to "log2".
#' See \code{\link[ggplot2]{scale_fill_gradientn}} for other accepted values.
#' Set to NULL for no transformation.
#' @param colours character vector, vector of colours to use for n-colour
#' gradient.
#' Default to \code{c("#000066", "#ffffbf", "#990000")}.
#'
#' @return A \code{ggplot} object.
#' @examples
#' object <- HiCDOCExample()
#' p <- plotInteractions(object, chromosome = 1, transform = "log2")
#' @export
plotInteractions <- function(
    object,
    chromosome,
    transform = "log2",
    colours = c("#000066", "#ffffbf", "#990000")
) {

    validateSlots(
        object,
        slots = c(
            "interactions",
            "conditions",
            "totalBins",
            "binSize",
            "positions"
        )
    )
    chromosomeName <- validateNameOrId(object, chromosome, "chromosomes")
    if (is.null(transform)) transform <- "Identity"

    positions <-
        object@positions %>%
        dplyr::filter(chromosome == chromosomeName)
    interactions <-
        object@interactions %>%
        dplyr::filter(chromosome == chromosomeName & interaction > 0) %>%
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
        labs(title = paste("Chromosome:", chromosomeName), x = "", y = "") +
        scale_fill_gradientn(
            colours = colours,
            trans = transform,
            name = "Intensity",
            na.value = "transparent"
        )
    return(plot)
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
    validateSlots(object, slots = c("interactions", "binSize"))
    data <-
        object@interactions %>%
        dplyr::mutate(distance = (bin.2 - bin.1) * object@binSize)
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
#' plotConcordanceDistribution(object)
#' @export
plotConcordanceDistribution <- function(object) {
    validateSlots(
        object,
        slots = c(
            "interactions",
            "differences",
            "concordances"
        )
    )

    changed <-
        object@differences %>%
        dplyr::select(-pvalue.adjusted) %>%
        dplyr::mutate(changed = "T")

    differences <-
        object@concordances %>%
        dplyr::group_by(.dots = c("chromosome", "bin", "condition")) %>%
        dplyr::summarise(median = stats::median(concordance)) %>%
        dplyr::ungroup() %>%
        tidyr::spread(condition, median) %>%
        dplyr::mutate(difference = `2` - `1`) %>%
        dplyr::select(-c(`1`, `2`)) %>%
        dplyr::left_join(changed, by = c("chromosome", "bin")) %>%
        dplyr::mutate(changed = tidyr::replace_na(changed, "F"))

    p <-
        ggplot(differences, aes(x = difference, fill = changed)) +
        geom_histogram() +
        labs(
            x = "Concordance",
            title = "Distribution of the differences of concordances"
        )
    return(p)
}

#' Plot boxplots of the diagonal ratios in each compartment,
#' for a given chromosome.
#'
#' @param object an \code{HiCDOCDataSet} object
#' @param chromosome character or numeric value, name or number of chromosome
#' @return A list of \code{ggplot}, one for each chromosome.
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotSelfInteractionRatios(object, 1)
#' @export
plotSelfInteractionRatios <- function(object, chromosome) {
    validateSlots(object, slots = c("selfInteractionRatios", "compartments"))
    chromosomeName <- validateNameOrId(object, chromosome, "chromosomes")

    data <-
        object@selfInteractionRatios %>%
        dplyr::left_join(
            object@compartments,
            by = c("chromosome", "condition", "bin")
        ) %>%
        dplyr::filter(chromosome == chromosomeName)

    plot <-
        ggplot(data, aes(x = compartment, y = ratio)) +
        geom_jitter(aes(color = compartment)) +
        geom_boxplot(
            outlier.colour = NA,
            fill = NA,
            colour = "grey20"
        ) +
        labs(
            color = "Compartment",
            x = "Compartment",
            y = "Interaction difference",
            title = paste0(
                "Differences between self-interactions ",
                "and other interactions in chromosome ",
                chromosomeName
            )
        )
    return(plot)
}


#' Plot the centroids PCA for a given chromosome.
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
plotCentroids <- function(object, chromosome, size = 2) {
    validateSlots(object, slots = c("centroids"))
    chromosomeName <- validateNameOrId(object, chromosome, "chromosomes")

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
