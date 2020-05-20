#' @export
filterSmallChromosomes <- function(object) {

  bigChromosomes <- makeSymmetric(object@interactions) %>%
    select(chromosome, `position.1`) %>%
    distinct() %>%
    group_by(chromosome) %>%
    summarize(nBins = n()) %>%
    ungroup() %>%
    filter(nBins >= object@minLength) %>%
    pull(chromosome) %>%
    as.character()

  #smallChromosomes <- names(which(object@totalBins < object@minLength))
  smallChromosomes <- object@chromosomes[!object@chromosomes %in% bigChromosomes]

  object@interactions %<>%
    filter(chromosome %in% bigChromosomes) %>%
    mutate(chromosome = factor(chromosome))

  object@chromosomes <- bigChromosomes
  object@totalBins <- object@totalBins[bigChromosomes]
  object@weakBins  <- object@weakBins[bigChromosomes]

  message(
    "Removed ",
    length(smallChromosomes),
    " chromosome",
    if (length(smallChromosomes) != 1) "s"
  )

  return (object)
}
