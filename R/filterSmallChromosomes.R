#' Suppress the small chromosomes
#'
#' The function return an HiCDOCExp object, with only the big chromosomes.
#' The big chromosomes must have a totalBins > object@minLength.
#' The interactions are reduced to keep only those corresponding 
#' to the remaining chromosomes 
#' 
#' @param object A HiCDOC object
#'
#' @return A HiCDOC object
#' @export
#'
#' @examples 
#' object <- HiCDOCExample()
#' object@chromosomes
#' object <- filterSmallChromosomes(object)
#' object@chromosomes

filterSmallChromosomes <- function(object) {

  bigChromosomes <- vapply(object@totalBins, function(x) x >= object@minLength, FUN.VALUE = TRUE)
  bigChromosomes <- names(bigChromosomes)[bigChromosomes == TRUE]
  bigChromosomes <- mixedsort(bigChromosomes)

  object@interactions %<>%
    filter(chromosome %in% bigChromosomes) %>%
    mutate(chromosome = factor(chromosome))

  object@chromosomes <- bigChromosomes
  object@totalBins <- object@totalBins[object@chromosomes]
  object@weakBins  <- object@weakBins[object@chromosomes]

  message(
    "Kept ",
    length(bigChromosomes),
    " chromosome",
    if (length(bigChromosomes) != 1) "s"
  )

  return (object)
}
