checkParameters <- function(value) {
    # TODO
}

makeFullMatrix <- function(data) {
  data %>% filter(`position 1` != `position 2`) %>%
    rename(tmp = `position 1`, `position 1` = `position 2`) %>%
    rename(`position 2` = tmp) %>%
    bind_rows(data)
}

#' @export
removeSmallChromosomes <- function(object) {
  bigChromosomes <- object@interactions %>%
    group_by(chromosome, replicate) %>%
    summarise(totalBins = n()) %>%
    summarise(length = sqrt(mean(totalBins))) %>%
    filter(length >= object@minLength) %>%
    pull(chromosome)
  message("Keeping ", length(bigChromosomes), " chromosomes.")
  bigChromosomes <- sort(bigChromosomes)
  object@chromosomes <- droplevels(bigChromosomes)
  object@interactions %<>%
    filter(chromosome %in% bigChromosomes) %>%
    mutate(chromosome = droplevels(chromosome))
  return (object)
}

sparseInteractionsToMatrix <- function(interactions, totalBins) {
  result <- matrix(0, nrow = totalBins, ncol = totalBins)
  data <- as.matrix(sapply(interactions, as.numeric))
  result[ data[, 1:2] ] <- data[, 3]
  result <- result + t(result) - diag(diag(result))
  if (!isSymmetric(result)) {
    stop("Matrix is not symmetric.")
  }
  return (result)
}
