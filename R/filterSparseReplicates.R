#' @description
#' Removes sparse replicates of a given chromosome.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' The name of a chromosome.
#' @param threshold
#' The minimum percentage of non-zero interactions for a replicate to be kept.
#'
#' @return
#' A list of:
#' - The sparse condition names repeated along the sparse replicates.
#' - The sparse replicate names repeated along the sparse conditions.
#' - The filtered interactions.
#'
#' @keywords internal
#' @noRd
.filterSparseReplicatesOfChromosome <- function(
    assay,
    diagonal,
    chromosomeName,
    totalBins,
    threshold,
    validAssay,
    conditions,
    replicates
) {

    
    filledAssay <- diagonal*(!is.na(assay) &  assay > 0)
    filledPct <- colSums(filledAssay) / (totalBins * totalBins)
    toRemove <- which(filledPct < threshold)
    toRemove <- toRemove[toRemove %in% validAssay]
    
    if (length(toRemove) > 0) {
        message(
            paste(
                "Removed interactions matrix of chromosome ",
                chromosomeName,
                ", condition ",
                conditions[toRemove],
                ", replicate ",
                replicates[toRemove],
                " filled at ",
                round(filledPct[toRemove], digits = 5) * 100,
                "%.",
                collapse = "\n", sep = ""
            )
        )
    }
    
    assay[,toRemove] <- NA
    return(assay)
}

#' @title
#' Filter sparse replicates.
#'
#' @description
#' Removes chromosome replicates whose percentage of non-zero interactions is
#' smaller than the threshold.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param threshold
#' The minimum percentage of non-zero interactions for a chromosome replicate to
#' be kept. If a chromosome replicate's percentage of non-zero interactions is
#' lower than this value, it is removed. Defaults to
#' \code{object$smallChromosomeThreshold} which is originally set to
#' \code{defaultHiCDOCParameters$smallChromosomeThreshold = 30\%}.
#'
#' @return
#' A filtered \code{\link{HiCDOCDataSet}}.
#'
#' @seealso
#' \code{\link{filterSmallChromosomes}},
#' \code{\link{filterWeakPositions}},
#' \code{\link{HiCDOC}}
#'
#' @examples
#' data(exampleHiCDOCDataSet)
#' object <- exampleHiCDOCDataSet
#'
#' object <- filterSparseReplicates(object)
#'
#' @export
filterSparseReplicates <- function(object, threshold = NULL) {

    .validateSlots(
        object,
        slots = c(
            "chromosomes",
            "parameters"
        )
    )

    if (!is.null(threshold)) {
        object@parameters$sparseReplicateThreshold <- threshold
    }
    object@parameters <- .validateParameters(object@parameters)
    threshold <- object@parameters$sparseReplicateThreshold

    message(
        "Keeping replicates filled with at least ",
        threshold * 100,
        "% non-zero interactions."
    )
    
    diagonal <- InteractionSet::anchors(object)
    diagonal <- diagonal$first == diagonal$second
    diagonal <- 2 - 1*diagonal
    diagonals <- split(diagonal, 
                       SummarizedExperiment::mcols(object)$Chr, drop=FALSE)
    assayChromosomes <- S4Vectors::split(
        SummarizedExperiment::assay(object), 
        SummarizedExperiment::mcols(object)$Chr, drop=FALSE)
    
    resultAssay <-
        pbapply::pbmapply(function(a, d, c, t, v)
            .filterSparseReplicatesOfChromosome(a, d, c, t, threshold, v, 
                                                object$condition, 
                                                object$replicat),
            assayChromosomes,
            diagonals,
            object@chromosomes,
            object@totalBins,
            object@validAssay
        )
    
    assay <- do.call("rbind",resultAssay)
    if(nrow(assay)!=nrow(object)){
        stop("Something went wrong")
    }
    
    SummarizedExperiment::assay(object)  <- assay
    newValidAssay <- .determineValids(object)
    badChromosomes <-
        vapply(
            newValidAssay,
            function(x) length(x) == 0,
            FUN.VALUE = TRUE
        )
    newValidAssay[badChromosomes] <- list(NULL)
    deleted <- length(unlist(object@validAssay)) - length(unlist(newValidAssay))
    object@validAssay <- newValidAssay
    
    rowsTosuppress <- (rowSums(SummarizedExperiment::assay(object), na.rm=TRUE) == 0)
    if(sum(rowsTosuppress)>0){
        object <- object[!rowsTosuppress,]
        object <- InteractionSet::reduceRegions(object)
    }
    message(
        "Removed ",
        deleted,
        " replicate",
        if (deleted != 1) "s",
        " in total."
    )
    
    if (nrow(object) == 0) {
        warning("No data left!", call. = FALSE)
    }

    return(object)
}
