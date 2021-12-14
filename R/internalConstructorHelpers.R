#' @description
#' Determines the number of bins per chromosome.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A list of the number of bins per chromosome.
#'
#' @keywords internal
#' @noRd
.determineChromosomeSizes <- function(object) {
    totalBins <- S4Vectors::runLength(
        GenomeInfoDb::seqnames(InteractionSet::regions(object)))
    names(totalBins) <- object@chromosomes
    return(totalBins)
}

#' @description
#' Determines the valid conditions and replicates by chromosomes (not empty)
#'
#' @param interactions
#' An InteractionSet object
#'
#' @return
#' A list of length 2, validConditions and validreplicates.
#'
#' @keywords internal
#' @noRd
.determineValids <- function(object) {
    valids <- S4Vectors::split(SummarizedExperiment::assay(object), 
                               S4Vectors::mcols(object)$Chr, drop=FALSE)
    valids <- lapply(valids, colSums, na.rm=TRUE)
    valids <- lapply(valids, function(x) which(x>0 & !is.na(x)))
    return(valids)
}
#' @description
#' Fills parameters and slots describing the data. Called by a
#' \code{\link{HiCDOCDataSet}} constructor.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}} ready for analysis.
#'
#' @keywords internal
#' @noRd
.fillHiCDOCDataSet <- function(object) {

    # Fill all other slots than interactionSet part
    
    # Chromosomes and their size (max bin)
    object@chromosomes <- gtools::mixedsort(
        GenomeInfoDb::seqlevels(object@interactions))
    object@totalBins <- .determineChromosomeSizes(object)
    object@parameters <- defaultHiCDOCParameters
    
    # Valid conditions and replicats by chromosome (==not empty)
    # maybe do a function for valid conditions and replicats ?
    valids <- .determineValids(object)
    object@validAssay <-valids

    # Weakbins
    object@weakBins <- vector("list", length(object@chromosomes))
    names(object@weakBins) <- object@chromosomes

    return(object)
}
