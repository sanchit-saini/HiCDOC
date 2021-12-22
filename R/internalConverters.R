#' @description
#' Fill \code{\link{InteractionSet}} with possibly missing values
#'
#' @param iset
#' An \code{\link{InteractionSet}}.
#' @param isetUnion
#' The full \code{\link{InteractionSet}}.
#' @param fill
#' Fill missing values with this.
#'
#' @return
#' The full \code{\link{InteractionSet}}.
#'
#' @keywords internal
#' @noRd
.fillInteractionSet <- function(iset, isetUnion, fill=NA){
    over <- GenomicRanges::match(iset, isetUnion)
    nc <- ncol(iset)
    newassays <- matrix(rep(fill, length(isetUnion)*nc), ncol=nc)
    newassays[over,] <- SummarizedExperiment::assay(iset)
    return(InteractionSet::InteractionSet(
        newassays, 
        isetUnion, 
        colData = SummarizedExperiment::colData(iset)))
}

#' @description
#' Merge two different \code{\link{InteractionSet}}.
#'
#' @param iset1
#' The first \code{\link{InteractionSet}}.
#' @param iset2
#' The second \code{\link{InteractionSet}}.
#' @param fill
#' Fill missing values with this.
#'
#' @return
#' The merged \code{\link{InteractionSet}}.
#'
#' @keywords internal
#' @noRd
.mergeInteractionSet <- function(iset1, iset2, fill=NA){
    unionInteractions <- GenomicRanges::union(InteractionSet::interactions(iset1),
                                              InteractionSet::interactions(iset2))
    # Complete InteractionSets
    iset1 <- .fillInteractionSet(iset1, unionInteractions, fill)
    iset2 <- .fillInteractionSet(iset2, unionInteractions, fill)

    # Merge
    newiset <- BiocGenerics::cbind(iset1, iset2)
    return(newiset)
}
