checkParameters <- function(value) {
  # TODO
}

#' Title
#'
#' @param chr current chromosome
#' @param interactions_chr interaction matrix for the chromosome
#' @param totalBins_chr totalBins for the chromosome
#' @param binSize binSize of the object (object@binSize)
#' @param replicates replicates of the object (as given by object@replicates)
#' @param conditions conditions of the object, with repetitions (as given by object@conditions)
#'
#' @return the full interaction matrix (filled with 0 values) for the chromosome
#' @export
#'
#' @examples
fullInteractionsChr <- function(chr, interactionsChr, totalBinsChr, binSize, replicates, conditions){
  positionsChr <- (seq_len(totalBinsChr) - 1) * binSize

  # Duplicate positions
  interactionsChr <- interactionsChr %>% filter(value>0)
  interactionsChr <- interactionsChr %>%
    filter(position.1 != position.2) %>%
    rename(position.2 = position.1, position.1 = position.2) %>%
    bind_rows(interactionsChr) %>%
    group_by(condition, replicate, position.1, position.2) %>%
    mutate(value = mean(value))

  fullmatrixChr <- tibble("condition" = factor(conditions), "replicate"= factor(replicates)) %>%
    expand(nesting(condition, replicate),
           "position.1" = positionsChr,
           "position.2" = positionsChr)

  fullmatrixChr <- fullmatrixChr %>%
    dplyr::left_join(interactionsChr, by = c("condition", "replicate", "position.1", "position.2")) %>%
    dplyr::mutate(value = ifelse(is.na(value)==T, 0, value)) %>%
    select(chromosome, condition, replicate, position.1, position.2, value)
  return(fullmatrixChr)
}


#' Full Interactions Matrix
#'
#' From an interactions matrix, that can be sparse and an upper of lower matrix, 
#' generate the full interactions matrix, symetric and filling gaps with 0.
#' @param object A \code{HiCDOCExp} object
#'
#' @return a tibble, long version of the full interactions matrix.
fullInteractions <- function(object) {
  interactions <- purrr::map_dfr(object@chromosomes,
                               function(x)
                                 fullInteractionsChr(
                                   x,
                                   object@interactions[object@interactions$chromosome == x,],
                                   object@totalBins[[x]],
                                   object@binSize,
                                   object@replicates,
                                   object@conditions
                                 )) %>%
  mutate(chromosome = factor(chromosome, levels = object@chromosomes))
  
  # # Duplicate positions
  # interactions <- object@interactions %>%
  #   filter(position.1 != position.2) %>%
  #   rename(temp = position.1, position.1 = position.2) %>%
  #   rename(position.2 = temp) %>%
  #   bind_rows(object@interactions)
  # 
  # # Fill missing positions with zeros
  # interactions <- tibble(
  #   chromosome = rep(
  #     as.factor(object@chromosomes),
  #     lapply(object@totalBins, function(x) {
  #       x * x * length(object@replicates)
  #     })
  #   ),
  #   condition = as.factor(flatten_chr(map(object@totalBins, function(x) {
  #     rep(object@conditions, each = x * x)
  #   }))),
  #   replicate = as.factor(flatten_chr(map(object@totalBins, function(x) {
  #     rep(object@replicates, each = x * x)
  #   }))),
  #   position.1 = flatten_int(map(object@totalBins, function(x) {
  #     rep(0:(x - 1), length(object@replicates) * x) * object@binSize
  #   })),
  #   position.2 = flatten_int(map(object@totalBins, function(x) {
  #     rep(0:(x - 1),
  #         each = x,
  #         times = length(object@replicates)) * object@binSize
  #   })),
  #   value      = 0.0,
  # ) %>% anti_join(
  #   interactions,
  #   by = c(
  #     "chromosome",
  #     "condition",
  #     "replicate",
  #     "position.1",
  #     "position.2"
  #   )
  # ) %>%
  #   bind_rows(interactions) %>%
  #   arrange(chromosome, condition, replicate, position.1, position.2)
  
  return (interactions)
}

sparseInteractionsToMatrix <- function(object,
                                       chromosomeId,
                                       conditionId,
                                       replicateId,
                                       filter = FALSE) {
  totalBins <- object@totalBins[[chromosomeId]]
  
  interactions <- object@interactions %>%
    filter(chromosome == chromosomeId) %>%
    filter(condition == conditionId) %>%
    filter(replicate == replicateId) %>%
    mutate(bin.1 = position.1 / object@binSize + 1,
           bin.2 = position.2 / object@binSize + 1) %>%
    select(bin.1, bin.2, value) %>%
    filter(value != 0) %>%
    as.matrix()
  
  if (nrow(interactions) == 0) {
    message("Warning: interaction matrix is empty")
    return(matrix(0, nrow = 0, ncol = 0))
  }
  if (is.numeric(interactions) == F) {
    stop("Error: non numeric matrix of interactions", call. = TRUE)
  }
  result <- matrix(0, nrow = totalBins, ncol = totalBins)
  result[interactions[, c(1, 2)]] <- interactions[, 3]
  result <- result + t(result) - diag(diag(result))
  if (!isSymmetric(result)) {
    stop("Matrix is not symmetric.")
  }
  
  if (filter && length(object@weakBins[[chromosomeId]]) > 0) {
    result = result[-object@weakBins[[chromosomeId]], -object@weakBins[[chromosomeId]]]
  }
  
  return (result)
}

matrixToSparseInteractions <- function(m,
                                       object,
                                       chromosomeId,
                                       conditionId,
                                       replicateId) {
  totalBins <- object@totalBins[[chromosomeId]]
  
  if (nrow(m) < totalBins) {
    refilled <- matrix(0, nrow = totalBins, ncol = totalBins)
    refilled[-object@weakBins[[chromosomeId]], -object@weakBins[[chromosomeId]]] <-
      m
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
  
  full_join(diagonal,
            offDiagonal,
            by = c("chromosome", "position", "condition")) %>%
    right_join(object@compartments,
               by = c("chromosome", "position", "condition")) %>%
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
    spread(key = compartment,
           value = value,
           fill = 0) %>%
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
      concordanceDifferences[rep(seq_len(nrow(concordanceDifferences)),
                                 totalReplicates),] %>%
        # Arrange by chromosome and position to join with the duplicated rows
        arrange(chromosome, position) %>%
        rename_with(function(old_colnames)
          paste(old_colnames, 2, sep = "."))
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
      compartmentComparisons[rep(seq_len(nrow(compartmentComparisons)),
                                 totalConditions),] %>%
        # Arrange by chromosome and position to join with the duplicated rows
        arrange(chromosome, position) %>%
        rename_with(function(old_colnames)
          paste(old_colnames, 2, sep = "."))
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
    mutate(H0_value = if_else(compartment.1 == compartment.2,
                              value,
                              NA_real_)) %>%
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
    select(chromosome,
           start,
           end,
           condition.1,
           condition.2,
           pvalue,
           padj,
           direction) %>%
    arrange(order(mixedsort(chromosome)), start, end, condition.1, condition.2)
  
  return(object)
}

#' Balanced condition levels, pasted with replicates
#'
#' @param conditions conditions vector, with repetitions, as in object@conditions
#'
#' @return Character vector filled to get the same number of replicates by condition, 
#' to use in \code{plotInteractionsMatrix()}
#'
#' @examples
CompleteCondLevels <- function(conditions, replicates){
  repcond <- paste(conditions, replicates)
  cond <- unique(conditions)
  lc <- purrr::map_int(cond, function(x) length(conditions[conditions==x]))
  lctofill <- max(lc) - lc
  completecond <- c(mapply(function(x, y, z) c(rep(x, y), rep(x, z)), cond, lc, lctofill, SIMPLIFY = T))
  return(completecond)
}



