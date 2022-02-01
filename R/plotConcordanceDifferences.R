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
#' data(exampleHiCDOCDataSetProcessed)
#' plotConcordanceDifferences(exampleHiCDOCDataSetProcessed)
#'
#' @export
plotConcordanceDifferences <- function(object) {
    .validateSlots(object,
                   slots = c("comparisons"))
    
    differences <-
        object@comparisons[, changed := data.table::fifelse(compartment.1 == compartment.2, "FALSE", "TRUE")]
    
    plot <-
        ggplot(differences, aes(x = difference, fill = changed)) +
        geom_histogram() +
        labs(x = "Concordance",
             title = "Distribution of concordance differences")
    return(plot)
}
