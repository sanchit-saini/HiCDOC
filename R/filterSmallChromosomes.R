#' @export
filterSmallChromosomes <- function(object) {

  smallChromosomes <- names(which(object@totalBins < object@minLength))

  object@interactions %<>%
    filter(!(chromosome %in% smallChromosomes)) %>%
    mutate(chromosome = factor(chromosome))

  object@chromosomes <- object@chromosomes[!object@chromosomes %in% smallChromosomes]
  object@totalBins[which(names(object@totalBins) %in% smallChromosomes)] <- NULL
  object@weakBins[which(names(object@weakBins) %in% smallChromosomes)] <- NULL

  message(
    "Removed ",
    length(smallChromosomes),
    " chromosome",
    if (length(smallChromosomes) != 1) "s"
  )

  return (object)
}
