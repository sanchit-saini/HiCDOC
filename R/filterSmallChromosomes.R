#' Suppress the small chromosomes
#'
#' The function return an HiCDOCExp object, with only the big chromosomes.
#' The big chromosomes must have a totalBins >= minLength.
#' The interactions are reduced to keep only those corresponding
#' to the remaining chromosomes.
#'
#' @param object A HiCDOC object
#' @param minLength Numeric value, default to 100. The minimum chromosome
#' size (in number of bins), to be kept.
#'
#' @return A HiCDOC object
#' @export
#'
#' @examples
#' object <- HiCDOCExample()
#' chromosomes(object)
#' object <- filterSmallChromosomes(object)
#' chromosomes(object)

filterSmallChromosomes <- function(object, minLength = 100) {
    message("Keeping only the chromosomes with ", minLength, " bins or more")
    bigChromosomes <- vapply(object@totalBins,
                             function(x)
                                 x >= minLength,
                             FUN.VALUE = TRUE)
    bigChromosomes <- names(bigChromosomes)[bigChromosomes == TRUE]
    bigChromosomes <- mixedsort(bigChromosomes)

    object@interactions %<>%
        dplyr::filter(chromosome %in% bigChromosomes) %>%
        dplyr::mutate(chromosome = factor(chromosome))

    object@chromosomes <- bigChromosomes
    object@totalBins <- object@totalBins[object@chromosomes]
    object@weakBins  <- object@weakBins[object@chromosomes]

    message("Kept ",
            length(bigChromosomes),
            " chromosome",
            if (length(bigChromosomes) != 1) "s"
    )

    return (object)
}
