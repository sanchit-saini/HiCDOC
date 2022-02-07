#' modeVector Extract the mode of vector.
#'
#' @param x A vector
#'
#' @return The mode of a vector
#' @export
#'
#' @examples
#' modeVector(c(1, 2, 2, 2, 4))
#' @noRd
modeVector <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

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
    tabChromosomes <- as.data.table(InteractionSet::regions(object))
    tabChromosomes[, minIndex := min(index), by = .(seqnames)]
    # minStart to correct chromosomes not starting at the 0 position
    tabChromosomes[index == minIndex, minStart := round(start / width), by =
                       .(seqnames)]
    # Comuting chromosome entire size
    tabChromosomes <-
        tabChromosomes[, .(binSize = max(index) - min(index) + 1 + max(minStart, na.rm =
                                                                           TRUE)),
                       by = .(seqnames)]
    totalBins <- tabChromosomes$binSize
    names(totalBins) <- tabChromosomes$seqnames
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
                               S4Vectors::mcols(object)$Chr,
                               drop = FALSE)
    valids <- lapply(valids, colSums, na.rm = TRUE)
    valids <- lapply(valids, function(x)
        which(x > 0 & !is.na(x)))
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
    # Reduce the levels in interaction part
    object <- InteractionSet::reduceRegions(object)
    objectRegions <- InteractionSet::regions(object)
    ChrNames <-
        unique(as.character(GenomeInfoDb::seqnames(objectRegions)))
    ChrNames <- gtools::mixedsort(ChrNames)
    GenomeInfoDb::seqlevels(InteractionSet::regions(object),
                            pruning.mode = "coarse") <- ChrNames
    
    # Add chromosome column for split purpose
    Chr <-
        GenomeInfoDb::seqnames(InteractionSet::anchors(object, "first"))
    Chr <- S4Vectors::Rle(factor(Chr, levels = ChrNames))
    S4Vectors::mcols(object) <-  S4Vectors::DataFrame("Chr" = Chr)
    
    # Sorting interactions and assay
    ids <- InteractionSet::anchors(object, id = TRUE)
    neworder <- order(Chr, ids$first, ids$second)
    object <- object[neworder, ]
    
    # Fill all other slots than interactionSet part
    # Chromosomes and their size (max bin)
    object@chromosomes <- ChrNames
    object@totalBins <- .determineChromosomeSizes(object)
    object@parameters <- defaultHiCDOCParameters
    
    # Valid conditions and replicats by chromosome (==not empty)
    # maybe do a function for valid conditions and replicats ?
    valids <- .determineValids(object)
    object@validAssay <- valids
    
    # Weakbins
    object@weakBins <- vector("list", length(object@chromosomes))
    names(object@weakBins) <- object@chromosomes
    
    return(object)
}
