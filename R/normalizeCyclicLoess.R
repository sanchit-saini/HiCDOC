#' @export
normalizeCyclicLoess <- function(object) {
  matrices <- object@interactionMatrix %>%
    mutate(chromosome = as.integer(chromosome)) %>%
    unite("condRep", c(condition, replicate)) %>%
    split(.$condRep) %>%
    map(~ select(.x, -condRep))
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
                        names(matrices))
  object@interactionMatrix <- output %>%
    mutate(chromosome = factor(object@chromosomes[chromosome])) %>%
    gather(names(matrices), key = "condRep", value = "value") %>%
    separate(condRep, c("condition", "replicate")) %>%
    mutate(condition = factor(condition)) %>%
    mutate(replicate = factor(replicate)) %>%
    select(chromosome, `position 1`, `position 2`, condition, replicate, value)
  return(object)
}
