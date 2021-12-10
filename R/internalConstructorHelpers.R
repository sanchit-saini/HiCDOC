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
        GenomeInfoDb::seqnames(InteractionSet::regions(object@interactions)))
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
.determineValids <- function(interactions) {
    valids <- S4Vectors::split(SummarizedExperiment::assay(interactions), 
                               S4Vectors::mcols(interactions)$Chr, drop=FALSE)
    valids <- lapply(valids, colSums, na.rm=TRUE)
    valids <- lapply(valids, function(x) (x>0 & !is.na(x)))
    validConditions <-
        lapply(valids, function(x) interactions$condition[x])
    validReplicates <-
        lapply(valids, function(x) interactions$replicat[x])
    return(list("validConditions" = validConditions, 
                "validReplicates" =validReplicates))
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

    # Fill all other slots than interactions
    .validateSlots(object, "interactions")

    # Chromosomes and their size (max bin)
    object@chromosomes <- GenomeInfoDb::seqlevels(object@interactions)
    object@totalBins <- .determineChromosomeSizes(object)
    object@parameters <- defaultHiCDOCParameters
    
    # Valid conditions and replicats by chromosome (==not empty)
    # maybe do a function for valid conditions and replicats ?
    valids <- .determineValids(object@interactions)
    object@validConditions <-valids$validConditions
    object@validReplicates <-valids$validReplicates

    # Weakbins
    object@weakBins <- vector("list", length(object@chromosomes))
    names(object@weakBins) <- object@chromosomes

    return(object)
}
