#' Suppress the small chromosomes
#'
#' The function return an HiCDOCDataSet object, with only the big chromosomes.
#' The big chromosomes must have a totalBins >= threshold.
#' The interactions are reduced to keep only those corresponding
#' to the remaining chromosomes.
#'
#' @param object A HiCDOCDataSet object
#' @param threshold Numeric value. The minimum chromosome
#' size (in number of bins), to be kept. If NULL, default to the first not
#' NULL of \code{object$smallChromosomeThreshold} and
#' \code{HiCDOCDefaultParameters$smallChromosomeThreshold}.
#'
#' @return A HiCDOCDataSet object
#' @export
#' @seealso \code{\link[HiCDOC]{filterWeakPositions}},
#' \code{\link[HiCDOC]{filterSparseChromosomes}} and
#' \code{\link[HiCDOC]{HiCDOC}} for the recommended pipeline.
#'
#' @examples
#' object <- HiCDOCExample()
#' chromosomes(object)
#' object <- filterSmallChromosomes(object)
#' chromosomes(object)
filterSmallChromosomes <- function(object, threshold = NULL) {
    if (!is.null(threshold)) {
        object@parameters$smallChromosomeThreshold <- threshold
    }
    object@parameters <- validateParameters(object@parameters)

    message(
        "Keeping chromosomes with at least ",
        object@parameters$smallChromosomeThreshold,
        " position",
        if (object@parameters$smallChromosomeThreshold != 1) "s",
        "."
    )

    bigChromosomes <-
        vapply(
            object@totalBins,
            function(x) x >= object@parameters$smallChromosomeThreshold,
            FUN.VALUE = TRUE
        )
    bigChromosomeNames <- names(bigChromosomes)[bigChromosomes]
    bigChromosomeNames <- gtools::mixedsort(bigChromosomeNames)
    smallChromosomeNames <-
        object@chromosomes[!(object@chromosomes %in% bigChromosomeNames)]

    object <- reduceHiCDOCDataSet(object, chromosomes = bigChromosomeNames)

    message(
        "Kept ",
        length(bigChromosomeNames),
        " chromosome",
        if (length(bigChromosomeNames) != 1) "s",
        if (length(bigChromosomeNames) > 0) ": " else ".",
        paste(bigChromosomeNames, collapse = ", ")
    )
    message(
        "Removed ",
        length(smallChromosomeNames),
        " chromosome",
        if (length(smallChromosomeNames) != 1) "s",
        if (length(smallChromosomeNames) > 0) ": " else ".",
        paste(smallChromosomeNames, collapse = ", ")
    )

    if (length(bigChromosomeNames) == 0) message("No data left!")

    return(object)
}
