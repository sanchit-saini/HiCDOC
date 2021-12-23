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

#' @description
#' Format the outputs produced by \code{detectCompartements}.
#' @param object
#' a HiCDOCDataSet object
#' @return
#' a HiCDOCDataSet object
#'
#' @keywords internal
#' @noRd
.formatDetectCompartment <- function(object){
    chr <- object@chromosomes
    cond <- sort(unique(object$condition))
    rep <- sort(unique(object$replicate))
    
    all.regions <- InteractionSet::regions(object)
    
    # Concordances
    object@concordances[,`:=`(chromosome = factor(chromosome, levels =  chr),
                              replicate = factor(replicate, levels = rep))]
    concordances <- object@concordances
    object@concordances <- all.regions[match(concordances$index, 
                                             S4Vectors::mcols(all.regions)$index)]
    S4Vectors::mcols(object@concordances) <- 
        S4Vectors::DataFrame(concordances[,.(condition, replicate, compartment, concordance)])
    
    # Centroids
    object@centroids[,`:=`(chromosome = factor(chromosome, levels =  chr),
                           condition = factor(condition, levels = cond))]
    
    # Differences
    object@differences[, chromosome := factor(chromosome, levels =  chr)]
    object@differences[, significance := ""]
    object@differences[pvalue.adjusted <= 0.05, significance := "*"]
    object@differences[pvalue.adjusted <= 0.01, significance := "**"]
    object@differences[pvalue.adjusted <= 0.001, significance := "***"]
    object@differences[pvalue.adjusted <= 0.0001, significance := "****"]

    differences <- object@differences
    object@differences <- all.regions[match(differences$index, 
                                            S4Vectors::mcols(all.regions)$index)]
    S4Vectors::mcols(object@differences) <- 
        S4Vectors::DataFrame(differences[,.(condition.1, condition.2, 
                                            pvalue, pvalue.adjusted, 
                                            direction, significance)])
    
    # Compartments
    object@compartments[, chromosome := factor(chromosome, levels =  chr)]
    compartments <- object@compartments
    object@compartments <- all.regions[match(compartments$index, 
                                             S4Vectors::mcols(all.regions)$index)]
    S4Vectors::mcols(object@compartments) <- 
        S4Vectors::DataFrame(compartments[,.(condition, compartment)])
    
    # Distances
    object@distances[,`:=`(chromosome = factor(chromosome, levels =  chr),
                           condition = factor(condition, levels = cond),
                           replicate = factor(replicate, levels = rep))]
    
    # Comparisons
    object@comparisons[,chromosome := factor(chromosome, levels =  chr)]
    
    # selfInteractionRatios
    object@selfInteractionRatios[,`:=`(chromosome = factor(chromosome, levels =  chr),
                                       condition = factor(condition, levels = cond),
                                       replicate = factor(replicate, levels = rep))]
    
    return(object)
}
