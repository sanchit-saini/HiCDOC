#' @export
findPValues <- function(object) {

  if (is.null(object@interactionMatrix)) {
    stop(paste0("Interaction matrix is not loaded yet.  ",
                "Please provide a matrix first."))
  }
  if (is.null(object@concordances)) {
    stop(paste0("Concordance is not computed.  ",
                "Please run 'detectConstrainedKMeans' first."))
  }
  if (is.null(object@compartments)) {
    stop(paste0("Compartments are not computed.  ",
                "Please run 'detectConstrainedKMeans' first."))
  }

  # Compute differences of concordances
  differences <- object@concordances %>%
    group_by(.dots = c("chromosome", "position", "condition")) %>%
    summarize(value = median(value)) %>%
    spread(condition, value) %>%
    mutate(value = `2` - `1`) %>%
    select(-c(`1`, `2`))

  # Compute number of compartments
  object@DIR <- object@compartments %>%
    group_by(.dots = c("chromosome", "position")) %>%
    summarize(nCompartments = n_distinct(value),
              compartment = first(value)) %>%
    ungroup() %>%
    left_join(differences, by = c("chromosome", "position")) %>%
    rename(start = position) %>%
    mutate(end = start + object@binSize) %>%
    mutate(percentRank = percent_rank(value)) %>%
    filter(nCompartments == 2) %>%
    mutate(pvalue = if_else(percentRank < 0.5,
                           2 * percentRank,
                           2 * (1 - percentRank))) %>%
    mutate(pvalue = if_else(pvalue < 0, 0, pvalue)) %>%
    mutate(pvalue = if_else(pvalue > 1, 1, pvalue)) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    mutate(padj = if_else(compartment == 1, -padj, padj)) %>%
    select(-c(nCompartments, compartment, value, percentRank, pvalue))

  return(object)
}
