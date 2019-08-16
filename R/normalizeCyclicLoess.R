#' @export
normalizeCyclicLoess <- function(object) {
  matrix <- object@interactionMatrix %>%
    spread(replicate, value)
  replicates <- as.vector(unique(object@interactionMatrix$replicate))
  matrices <- object@interactionMatrix %>%
    split(.$replicate) %>%
    map(~ select(.x, -replicate))
  hicexp <- make_hicexp(
    data_list = matrices,
    groups = object@conditions,
    remove.regions = NULL
  )
  normalized <- cyclic_loess(hicexp, parallel = TRUE)
  output <- as_tibble(hic_table(normalized)) %>%
    select(-D)
  colnames(output) <- c("chromosome",
                        "position 1",
                        "position 2",
                        replicates)
  object@interactionMatrix <- output %>%
    mutate(chromosome = object@chromosomes[chromosome]) %>%
    gather(replicates, key = "replicate", value = "value")
  return(object)
}
