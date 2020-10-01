checkParameters <- function(value) {
  # TODO
}

##- sparseInteractionsToFullInteractions -------------------------------------#
##----------------------------------------------------------------------------#
#' Build the full interactions tibble filled with zeros for a chromosome in a
#' condition and replicate.
#'
#' @param object        A \code{HiCDOCExp} object.
#' @param chromosomeId  A chromosome.
#' @param conditionId   A condition.
#' @param replicateId   A replicate.
#' @param filter        Whether to remove weak positions. Defaults to TRUE.
#'
#' @return An interactions tibble.
sparseInteractionsToFullInteractions <- function(
  object,
  chromosomeId,
  conditionId,
  replicateId,
  filter = TRUE
) {
  totalBins <- object@totalBins[[chromosomeId]]
  positions <- (seq_len(totalBins) - 1) * object@binSize

  interactions = object@interactions %>%
    filter(chromosome == chromosomeId) %>%
    filter(condition == conditionId) %>%
    filter(replicate == replicateId)

  mirrorInteractions = interactions %>%
    filter(position.1 != position.2) %>%
    rename(position.1 = position.2, position.2 = position.1)

  fullInteractions <- tibble(
      chromosome = chromosomeId,
      condition = conditionId,
      replicate = replicateId,
      position.1 = rep(positions, totalBins),
      position.2 = rep(positions, each = totalBins),
      value = 0
    ) %>%
    left_join(
      bind_rows(interactions, mirrorInteractions),
      by = c("chromosome", "condition", "replicate", "position.1", "position.2")
    ) %>%
    mutate(value = coalesce(value.y, value.x)) %>%
    select(-c(value.x, value.y))

  if (filter) {
    weakBins <- object@weakBins[[chromosomeId]]
    weakPositions <- (weakBins - 1) * object@binSize

    fullInteractions %<>%
     filter(!(position.1 %in% weakPositions)) %>%
     filter(!(position.2 %in% weakPositions))
  }

  return (fullInteractions)
}

##- sparseInteractionsToMatrix -----------------------------------------------#
##----------------------------------------------------------------------------#
#' Build the interaction matrix for a chromosome in a condition and replicate.
#'
#' @param object        A \code{HiCDOCExp} object.
#' @param chromosomeId  A chromosome.
#' @param conditionId   A condition.
#' @param replicateId   A replicate.
#' @param filter        Shrink the matrix by removing weak rows/columns.
#'
#' @return A matrix.
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
    select(bin.1, bin.2, value) %>%
    filter(value != 0) %>%
    as.matrix()

  if (nrow(interactions) == 0) {
    return(matrix(0, nrow = 0, ncol = 0))
  }
  if (!is.numeric(interactions)) {
    stop("Error: non numeric matrix of interactions.", call. = TRUE)
  }
  result <- matrix(0, nrow = totalBins, ncol = totalBins)
  result[interactions[, c(1, 2)]] <- interactions[, 3]
  result <- result + t(result) - diag(diag(result))
  if (!isSymmetric(result)) {
    stop("Matrix is not symmetric.")
  }

  if (filter && length(object@weakBins[[chromosomeId]]) > 0) {
    result <- result[
      -object@weakBins[[chromosomeId]],
      -object@weakBins[[chromosomeId]]
    ]
  }

  return (result)
}

##- matrixToSparseInteractions -----------------------------------------------#
##----------------------------------------------------------------------------#
#' Build the interactions tibble for a chromosome in a condition and replicate.
#'
#' @param m             A matrix.
#' @param object        A \code{HiCDOCExp} object.
#' @param chromosomeId  A chromosome.
#' @param conditionId   A condition.
#' @param replicateId   A replicate.
#'
#' @return An interactions tibble.
matrixToSparseInteractions <- function(
  m,
  object,
  chromosomeId,
  conditionId,
  replicateId
) {
  totalBins <- object@totalBins[[chromosomeId]]

  if (nrow(m) < totalBins) {
    refilled <- matrix(0, nrow = totalBins, ncol = totalBins)
    refilled[
      -object@weakBins[[chromosomeId]],
      -object@weakBins[[chromosomeId]]
    ] <- m
    m <- refilled
  }

  return (
    tibble(
      chromosome = chromosomeId,
      position.1 = (rep(seq(totalBins), each = totalBins) - 1) * object@binSize,
      position.2 = (rep(seq(totalBins), times = totalBins) - 1) * object@binSize,
      condition = conditionId,
      replicate = replicateId,
      value = as.vector(t(m))
    ) %>% filter(position.1 <= position.2 & value != 0)
  )
}

##- predictAB ----------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Use ratio between diagonal and off-diagonal interactions to determine which
#' clusters correspond to compartments A and B.
#'
#' @param object A \code{HiCDOCExp} object.
#'
#' @return A \code{HiCDOCExp} object, with diagonalRatios, and with A and B
#' labels replacing cluster numbers in centroids, compartments, distances and
#' concordances.
predictAB <- function(object) {

  chromosomeIds = rep(object@chromosomes, each = length(object@replicates))
  conditionIds = rep(object@conditions, length(object@chromosomes))
  replicateIds = rep(object@replicates, length(object@chromosomes))

  groups = cbind(chromosomeIds, conditionIds, replicateIds)

  object@diagonalRatios <- bind_rows(apply(groups, 1, function(group) {
    fullInteractions <- sparseInteractionsToFullInteractions(
      object,
      group[[1]],
      group[[2]],
      group[[3]]
    )

    diagonal <- fullInteractions %>%
      filter(position.1 == position.2) %>%
      select(-position.2) %>%
      rename(position = position.1) %>%
      rename(diagonal = value)

    offDiagonal <- fullInteractions %>%
      filter(position.1 != position.2) %>%
      select(-position.2) %>%
      rename(position = position.1) %>%
      group_by(position) %>%
      mutate(offDiagonal = median(value)) %>%
      ungroup() %>%
      select(-value) %>%
      distinct()

    diagonalRatios <- diagonal %>%
      left_join(
        offDiagonal,
        by = c("chromosome", "condition", "replicate", "position")
      ) %>%
      mutate(value = diagonal - offDiagonal) %>%
      select(-c(diagonal, offDiagonal))

    return (diagonalRatios)
  }))

  compartments <- object@compartments %>%
    rename(compartment = value) %>%
    left_join(
      object@diagonalRatios,
      by = c("chromosome", "condition", "position")
    ) %>%
    group_by(chromosome, condition, compartment) %>%
    summarize(value = median(value)) %>%
    ungroup() %>%
    spread(
      key = compartment,
      value = value,
      fill = 0
    ) %>%
    mutate(A = if_else(`1` >= `2`, 2, 1)) %>%
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

  return (object)
}

##- computePValues -----------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Get p-values for genomic positions whose assigned compartment switches
#' between two conditions:
#' 1. For each pair of replicates in different conditions, for each genomic
#'    position, compute the absolute difference between its concordances.
#' 2. For each pair of conditions, for each genomic position, compute the
#'    median of its concordance differences.
#' 3. For each pair of conditions, for each genomic position whose assigned
#'    compartment switches, rank its median against the empirical cumulative
#'    distribution of medians of all non-switching positions in that condition
#'    pair. Adjust the resulting p-value with the Benjamini–Hochberg procedure.
#'
#' @param object A \code{HiCDOCExp} object.
#'
#' @return A \code{HiCDOCExp} object, with differences and their p-values.
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
      rename_with(function(old_colnames) paste(old_colnames, 2, sep = "."))
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
      rename_with(function(old_colnames) paste(old_colnames, 2, sep = "."))
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
    mutate(
      H0_value = if_else(compartment.1 == compartment.2, value, NA_real_)
    ) %>%
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
    mutate(
      direction = factor(if_else(compartment.1 == "A", "A->B", "B->A"))
    ) %>%
    select(
      chromosome,
      start,
      end,
      condition.1,
      condition.2,
      pvalue,
      padj,
      direction
    ) %>%
    arrange(order(mixedsort(chromosome)), start, end, condition.1, condition.2)

  return (object)
}
