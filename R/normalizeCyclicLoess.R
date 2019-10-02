#' @export
normalizeCyclicLoess <- function(object) {
  matrices <- object@interactionMatrix %>%
    mutate(chromosome = as.integer(chromosome)) %>%
    unite("condRep", c(condition, replicate)) %>%
    mutate(condRep = factor(condRep, levels = paste(object@conditions, object@replicates, sep = "_"))) %>%
    split(.$condRep) %>%
    map(~ select(.x, -condRep))
  names(matrices) <- paste0("Sample", seq_along(object@conditions))
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
                        seq_along(object@replicates))
  object@interactionMatrix <- output %>%
    mutate(chromosome = factor(object@chromosomes[chromosome])) %>%
    gather(as.character(seq_along(object@replicates)), key = "condRep", value = "value") %>%
    mutate(condRep = factor(as.integer(condRep))) %>%
    mutate(condition = factor(object@conditions[condRep])) %>%
    mutate(replicate = object@replicates[condRep]) %>%
    select(chromosome, `position 1`, `position 2`, condition, replicate, value)
  return(object)
}
