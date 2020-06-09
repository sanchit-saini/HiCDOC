#' @export
filterSmallChromosomes <- function(object) {

  bigChromosomes <- fullInteractions(object) %>%
    select(chromosome, `position.1`) %>%
    distinct() %>%
    group_by(chromosome) %>%
    summarize(nBins = n()) %>%
    ungroup() %>%
    filter(nBins >= object@minLength) %>%
    pull(chromosome) %>%
    as.character()

  object@interactions %<>%
    filter(chromosome %in% bigChromosomes) %>%
    mutate(chromosome = factor(chromosome))

  object@chromosomes <- bigChromosomes[order(nchar(bigChromosomes), bigChromosomes)]
  object@totalBins <- object@totalBins[bigChromosomes]
  object@weakBins  <- object@weakBins[bigChromosomes]

  message(
    "Kept ",
    length(bigChromosomes),
    " chromosome",
    if (length(bigChromosomes) != 1) "s"
  )

  return (object)
}
