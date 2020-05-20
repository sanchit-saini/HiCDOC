#' @export
normalizeTechnicalBiases <- function(object) {

  chromosome_as_list <- levels(object@interactions$chromosome)

  matrices <- object@interactions %>%
    mutate(chromosome = as.integer(chromosome)) %>%
    group_split(condition, replicate) %>%
    map(function(x) select(x, -c(condition, replicate)))


  hicexp <- make_hicexp(
    data_list = matrices,
    groups = object@conditions,
    remove.regions = NULL,
    # Default filtering parameters explicitely specified here
    remove_zeros = FALSE, zero.p = 0.8, A.min = 5, filter = TRUE
  )

  normalized <- cyclic_loess(hicexp, parallel = FALSE)
  #normalized <- cyclic_loess(hicexp, parallel = TRUE)
  output <- hic_table(normalized) %>% as_tibble() %>% select(-D)

  colnames(output) <- c(
    "chromosome", "position.1", "position.2", seq_along(object@replicates)
  )

  object@interactions <- output %>%
    gather(
      as.character(seq_along(object@replicates)),
      key = "i",
      value = "value"
    ) %>%
    mutate(i = factor(as.integer(i))) %>%
    mutate(condition = factor(object@conditions[i])) %>%
    mutate(replicate = factor(object@replicates[i])) %>%
    mutate(chromosome = factor(chromosome_as_list[chromosome])) %>%
    select(chromosome, position.1, position.2, condition, replicate, value)

  return(object)
}
