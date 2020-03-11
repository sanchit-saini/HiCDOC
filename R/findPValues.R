#' @export
findPValues <- function(object) {

  if (is.null(object@interactions)) {
    stop(paste0(
      "Interaction matrix is not loaded yet.  ",
      "Please provide a matrix first."
    ))
  }
  if (is.null(object@concordances)) {
    stop(paste0(
      "Concordance is not computed.  ",
      "Please run 'detectConstrainedKMeans' first."
    ))
  }
  if (is.null(object@compartments)) {
    stop(paste0(
      "Compartments are not computed.  ",
      "Please run 'detectConstrainedKMeans' first."
    ))
  }


  # Compute median of differences between pairs of concordances
  # N.b. median of differences != difference of medians
  totalReplicates = length(object@replicates)

  concordanceDifferences <- object@concordances %>%
    arrange(chromosome, position, condition, replicate)

  concordanceDifferences %<>%
    # Duplicate each row "totalReplicates" times
    uncount(totalReplicates) %>%
    bind_cols(
      # Duplicate the table "totalReplicates" times
      concordanceDifferences[rep(
        seq_len(nrow(concordanceDifferences)),
        totalReplicates
      ),] %>%
      # Arrange by chromosome and position to join with the duplicated rows
      arrange(chromosome, position)
    ) %>%
    select(-chromosome1, -position1) %>%
    rename(condition.1 = condition, condition.2 = condition1) %>%
    filter(as.numeric(condition.1) < as.numeric(condition.2)) %>%
    group_by(chromosome, position, condition.1, condition.2) %>%
    summarise(value = median(abs(value - value1))) %>%
    ungroup()


  # Format compartments per pair of conditions
  totalConditions = length(unique(object@conditions))

  compartmentComparisons <- object@compartments %>%
    arrange(chromosome, position, condition)

  compartmentComparisons %<>%
    # Duplicate each row "totalConditions" times
    uncount(totalConditions) %>%
    bind_cols(
      # Duplicate the table "totalConditions" times
      compartmentComparisons[rep(
        seq_len(nrow(compartmentComparisons)),
        totalConditions
      ),] %>%
      # Arrange by chromosome and position to join with the duplicated rows
      arrange(chromosome, position)
    ) %>%
    select(-chromosome1, -position1) %>%
    rename(
      condition.1 = condition,
      condition.2 = condition1,
      compartment.1 = value,
      compartment.2 = value1
    ) %>%
    filter(as.numeric(condition.1) < as.numeric(condition.2))


  # Compute p-values for switching positions
  # P-values for a condition pair computed from the whole genome distribution
  object@differences <- compartmentComparisons %>%
    left_join(
      concordanceDifferences,
      by = c("chromosome", "position", "condition.1", "condition.2")
    ) %>%
    mutate(H0_value = if_else(
      compartment.1 == compartment.2,
      value,
      NA_real_
    )) %>%
    group_by(condition.1, condition.2) %>%
    mutate(quantile = ecdf(H0_value)(value)) %>%
    filter(compartment.1 != compartment.2) %>%
    mutate(pvalue = if_else(
      quantile < 0.5,
      2 * quantile,
      2 * (1 - quantile)
    )) %>%
    mutate(pvalue = if_else(pvalue < 0, 0, pvalue)) %>%
    mutate(pvalue = if_else(pvalue > 1, 1, pvalue)) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    ungroup() %>%
    rename(start = position) %>%
    mutate(end = start + object@binSize) %>%
    mutate(direction = factor(if_else(compartment.1 == "A", "A->B", "B->A"))) %>%
    select(chromosome, start, end, condition.1, condition.2, pvalue, padj, direction)

  return(object)
}
