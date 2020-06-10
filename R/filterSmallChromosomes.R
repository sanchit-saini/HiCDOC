#' @export
filterSmallChromosomes <- function(object) {

  bigChromosomes <- fullInteractions(object) %>%
    select(chromosome, `position.1`) %>%
    distinct() %>%
    group_by(chromosome) %>%
    summarize(totalBins = n()) %>%
    ungroup() %>%
    filter(totalBins >= object@minLength) %>%
    pull(chromosome) %>%
    as.character()

  object@interactions %<>%
    filter(chromosome %in% bigChromosomes) %>%
    mutate(chromosome = factor(chromosome))

  object@chromosomes <- mixedsort(bigChromosomes)
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
