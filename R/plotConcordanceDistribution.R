#' Plot the concordance, i.e. the relative distance of the genomic positions
#' with respect to the centroids.
#'
#' @param object an \code{HiCDOCDataSet} object
#' @return A list of \code{ggplot}, one for each chromosome.
#' @examples
#' object <- HiCDOCDataSetExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotConcordanceDistribution(object)
#' @export
plotConcordanceDistribution <- function(object) {
    .validateSlots(
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
