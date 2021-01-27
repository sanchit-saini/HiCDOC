#' Suppress the small chromosomes
#'
#' The function return an HiCDOCDataSet object, with only the big chromosomes.
#' The big chromosomes must have a totalBins >= minLength.
#' The interactions are reduced to keep only those corresponding
#' to the remaining chromosomes.
#'
#' @param object A HiCDOCDataSet object
#' @param minLength Numeric value, default to
#' HiCDOCDefaultParameters$minLengthChr. The minimum chromosome
#' size (in number of bins), to be kept.
#'
#' @return A HiCDOCDataSet object
#' @export
#'
#' @examples
#' object <- HiCDOCExample()
#' chromosomes(object)
#' object <- filterSmallChromosomes(object)
#' chromosomes(object)

filterSmallChromosomes <- function(object, minLength = NULL) {
    if(! is.null(minLength) ) object@parameters$minLengthChr <- minLength
    object@parameters <- checkParameters(object@parameters,
                                         c("minLengthChr"))
    message("Keeping only the chromosomes with ",
            object@parameters$minLengthChr, " bins or more")
    bigChromosomes <- vapply(object@totalBins,
                             function(x)
                                 x >= object@parameters$minLengthChr,
                             FUN.VALUE = TRUE)
    bigChromosomes <- names(bigChromosomes)[bigChromosomes == TRUE]
    bigChromosomes <- gtools::mixedsort(bigChromosomes)
    smallChr <- object@chromosomes[!(object@chromosomes %in% bigChromosomes)]
    
    object <- reduceHiCDOCDataSet(object, chromosomes = bigChromosomes)
    
    message("Kept ",
            length(bigChromosomes),
            " chromosome",
            if (length(bigChromosomes) > 1) "s"
    )
    message("Removed ",
            length(smallChr),
            " chromosome",
            if (length(smallChr) > 1) "s :",
            paste(smallChr, collapse = ", ")
    )
    
    return (object)
}
