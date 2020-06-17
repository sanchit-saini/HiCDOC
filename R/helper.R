checkParameters <- function(value) {
    # TODO
}

fullInteractions <- function(object) {

  # Duplicate positions
  interactions <- object@interactions %>%
    filter(position.1 != position.2) %>%
    rename(temp = position.1, position.1 = position.2) %>%
    rename(position.2 = temp) %>%
    bind_rows(object@interactions)

  # Fill missing positions with zeros
  interactions <- tibble(
    chromosome = rep(object@chromosomes, lapply(object@totalBins, function(x) {
      x*x*length(object@replicates)
    })),
    condition = as.vector(unlist(map(object@totalBins, function(x) {
      rep(object@conditions, each = x*x)
    }))),
    replicate = as.vector(unlist(map(object@totalBins, function(x) {
      rep(object@replicates, each = x*x)
    }))),
    position.1 = as.vector(unlist(map(object@totalBins, function(x) {
      rep(0:(x-1), length(object@replicates)*x) * object@binSize
    }))),
    position.2 = as.vector(unlist(map(object@totalBins, function(x) {
      rep(0:(x-1), each = x, times = length(object@replicates)) * object@binSize
    }))),
    value      = 0.0,
  ) %>% mutate(
    chromosome = factor(chromosome),
    condition = factor(condition),
    replicate = factor(replicate),
    position.1 = as.integer(position.1),
    position.2 = as.integer(position.2)
  ) %>% anti_join(
    interactions,
    by = c("chromosome", "condition", "replicate", "position.1", "position.2")
  ) %>%
    bind_rows(interactions) %>%
    arrange(chromosome, condition, replicate, position.1, position.2)

  return (interactions)
}

sparseInteractionsToMatrix <- function(
  object,
  chromosomeId,
  conditionId,
  replicateId,
  filter = FALSE
) {

  totalBins <- object@totalBins[[chromosomeId]]

  interactions <- object@interactions %>%
    filter(chromosome == chromosomeId) %>%
    filter(condition == conditionId) %>%
    filter(replicate == replicateId) %>%
    mutate(
      bin.1 = position.1 / object@binSize + 1,
      bin.2 = position.2 / object@binSize + 1
    ) %>%
    select(bin.1, bin.2, value)
  if (nrow(interactions) == 0) {
    message("Warning: interaction matrix is empty")
    return(matrix(, nrow = 0, ncol = 0))
  }
  result <- matrix(0, nrow = totalBins, ncol = totalBins)
  data <- as.matrix(sapply(interactions, as.numeric))
  result[ data[, 1:2] ] <- data[, 3]
  result <- result + t(result) - diag(diag(result))
  if (!isSymmetric(result)) {
    stop("Matrix is not symmetric.")
  }

  if (filter && length(object@weakBins[[chromosomeId]]) > 0) {
    result = result[
      -object@weakBins[[chromosomeId]],
      -object@weakBins[[chromosomeId]]
    ]
  }

  return (result)
}

matrixToSparseInteractions <- function(
  m,
  object,
  chromosomeId,
  conditionId,
  replicateId
) {

  totalBins <- object@totalBins[[chromosomeId]]

  if (nrow(m) < totalBins) {
    refilled <- matrix(0, nrow=totalBins, ncol=totalBins)
    refilled[
      -object@weakBins[[chromosomeId]],
      -object@weakBins[[chromosomeId]]
    ] <- m
    m <- refilled
  }

  return (tibble(
    chromosome = chromosomeId,
    position.1 = (rep(seq(totalBins), each = totalBins) - 1) * object@binSize,
    position.2 = (rep(seq(totalBins), times = totalBins) - 1) * object@binSize,
    condition = conditionId,
    replicate = replicateId,
    value = as.vector(t(m))
  ) %>% filter(position.1 <= position.2 & value != 0))
}

buildABComparison <- function(object) {
  diagonal <- object@interactions %>%
    filter(position.1 == position.2) %>%
    rename(position = position.1) %>%
    select(-position.2) %>%
    group_by(chromosome, position, condition) %>%
    summarise(diagValue = median(value)) %>%
    ungroup()
  offDiagonal <- fullInteractions(object) %>%
    filter(position.1 != position.2) %>%
    group_by(chromosome, position.1, condition, replicate) %>%
    summarise(value = mean(value)) %>%
    summarise(offDiagValue = median(value)) %>%
    ungroup() %>%
    rename(position = position.1)
  full_join(diagonal, offDiagonal, by = c("chromosome", "position", "condition")) %>%
    right_join( object@compartments, by = c("chromosome", "position", "condition")) %>%
    rename(compartment = value) %>%
    replace_na(list(diagValue = 0, offDiagValue = 0)) %>%
    mutate(diffValue = diagValue - offDiagValue) %>%
    select(-c(position, diagValue, offDiagValue)) %>%
    mutate(chromosome = factor(chromosome, levels = object@chromosomes))
}

predictABCompartments <- function(object) {
  compartments <- buildABComparison(object) %>%
    group_by(chromosome, condition, compartment) %>%
    summarise(value = median(diffValue)) %>%
    ungroup() %>%
    spread(key = compartment, value = value, fill = 0) %>%
    mutate(A = if_else(`1` >= `2`, 1, 2)) %>%
    select(-c(`1`, `2`))
  object@compartments %<>%
    left_join(compartments, by = c("chromosome", "condition")) %>%
    mutate(value = factor(if_else(value == A, "A", "B"))) %>%
    select(-c(A))
  object@concordances %<>%
    left_join(compartments, by = c("chromosome", "condition")) %>%
    mutate(change = if_else(A == 1, 1, -1)) %>%
    mutate(value = change * value) %>%
    select(-c(A, change))
  object@distances %<>%
    left_join(compartments, by = c("chromosome", "condition")) %>%
    mutate(cluster = factor(if_else(cluster == A, "A", "B"))) %>%
    select(-c(A))
  object@centroids %<>%
    left_join(compartments, by = c("chromosome", "condition")) %>%
    mutate(compartment = factor(if_else(compartment == A, "A", "B"))) %>%
    select(-c(A))
  return(object)
}

computePValues <- function(object) {

  # Compute median of differences between pairs of concordances
  # N.b. median of differences != difference of medians
  totalReplicates = length(object@replicates)

  concordanceDifferences <- object@concordances %>%
    arrange(chromosome, position, condition, replicate)

  concordanceDifferences %<>%
    # Duplicate each row "totalReplicates" times
    uncount(totalReplicates) %>%
    cbind(
      # Duplicate the table "totalReplicates" times
      concordanceDifferences[rep(
        seq_len(nrow(concordanceDifferences)),
        totalReplicates
      ),] %>%
      # Arrange by chromosome and position to join with the duplicated rows
      arrange(chromosome, position) %>%
      set_colnames(paste(colnames(concordanceDifferences), 2, sep="."))
    ) %>%
    select(-chromosome.2, -position.2) %>%
    rename(condition.1 = condition, value.1 = value) %>%
    filter(as.numeric(condition.1) < as.numeric(condition.2)) %>%
    group_by(chromosome, position, condition.1, condition.2) %>%
    summarise(value = median(abs(value.1 - value.2))) %>%
    ungroup()

  # Format compartments per pair of conditions
  totalConditions = length(unique(object@conditions))

  compartmentComparisons <- object@compartments %>%
    arrange(chromosome, position, condition)

  compartmentComparisons %<>%
    # Duplicate each row "totalConditions" times
    uncount(totalConditions) %>%
    cbind(
      # Duplicate the table "totalConditions" times
      compartmentComparisons[rep(
        seq_len(nrow(compartmentComparisons)),
        totalConditions
      ),] %>%
      # Arrange by chromosome and position to join with the duplicated rows
      arrange(chromosome, position) %>%
      set_colnames(paste(colnames(compartmentComparisons), 2, sep="."))
    ) %>%
    select(-chromosome.2, -position.2) %>%
    rename(
      condition.1 = condition,
      compartment.1 = value,
      compartment.2 = value.2
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
    mutate(pvalue = 1 - quantile) %>%
    mutate(pvalue = if_else(pvalue < 0, 0, pvalue)) %>%
    mutate(pvalue = if_else(pvalue > 1, 1, pvalue)) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    ungroup() %>%
    rename(start = position) %>%
    mutate(end = start + object@binSize) %>%
    mutate(direction = factor(if_else(compartment.1 == "A", "A->B", "B->A"))) %>%
    select(chromosome, start, end, condition.1, condition.2, pvalue, padj, direction) %>%
    arrange(order(mixedsort(chromosome)), start, end, condition.1, condition.2)

  return(object)
}
