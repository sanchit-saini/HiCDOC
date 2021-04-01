#' @title
#' Filter small chromosomes.
#'
#' @description
#' Removes chromosomes whose length (in number of positions) is smaller than the
#' threshold.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param threshold
#' The minimum length (number of positions) for a chromosome to be kept.
#' Defaults to \code{object$smallChromosomeThreshold} which is originally set to
#' \code{defaultHiCDOCParameters$smallChromosomeThreshold} = 100.
#'
#' @return
#' A filtered \code{\link{HiCDOCDataSet}}.
#'
#' @seealso
#' \code{\link{filterSparseReplicates}},
#' \code{\link{filterWeakPositions}},
#' \code{\link{HiCDOC}}
#'
#' @examples
#' data(HiCDOCDataSetExample)
#' chromosomes(HiCDOCDataSetExample)
#' object <- filterSmallChromosomes(HiCDOCDataSetExample)
#' chromosomes(object)
#'
#' @export
filterSmallChromosomes <- function(object, threshold = NULL) {

    .validateSlots(
        object,
        slots = c(
            "interactions",
            "chromosomes",
            "totalBins",
            "parameters"
        )
    )

    if (!is.null(threshold)) {
        object@parameters$smallChromosomeThreshold <- threshold
    }
    object@parameters <- .validateParameters(object@parameters)
    threshold <- object@parameters$smallChromosomeThreshold

    message(
        "Keeping chromosomes with at least ",
        threshold,
        " position",
        if (threshold != 1) "s",
        "."
    )

    bigChromosomes <-
        vapply(
            object@totalBins,
            function(totalBins) totalBins >= threshold,
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

    if (length(bigChromosomeNames) == 0) {
        warning("No data left!", call. = FALSE)
    }

    return(object)
}
