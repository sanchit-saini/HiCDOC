#' @title
#' Plot the distribution of concordance differences.
#'
#' @description
#' Plots the distribution of concordance differences, which are the differences
#' between concordances of each pair of replicates from different conditions. A
#' concordance can be understood as a confidence in a genomic position's
#' assigned compartment. Mathematically, it is the log ratio of a genomic
#' position's distance to each compartment's centroid, normalized by the
#' distance between both centroids, and min-maxed to a [-1,1] interval.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' data(exampleHiCDOCDataSet)
#' object <- HiCDOC(exampleHiCDOCDataSet)
#' plotConcordanceDifferences(object)
#'
#' @export
plotConcordanceDifferences <- function(object) {
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
        dplyr::mutate(changed = tidyr::replace_na(changed, "F")) %>%
        dplyr::filter(!is.na(difference))

    plot <-
        ggplot(differences, aes(x = difference, fill = changed)) +
        geom_histogram() +
        labs(
            x = "Concordance",
            title = "Distribution of concordance differences"
        )
    return(plot)
}
