euclideanDistance <- function(x, y) {
  sqrt(sum((x - y)^2))
}

distanceRatio <- function(x, centroids, eps=1e-10) {
  return (log(
    (euclideanDistance(x, centroids[[1]]) + eps)
    / (euclideanDistance(x, centroids[[2]]) + eps)
  ))
}

#' @export
detectConstrainedKMeans <- function(object) {

  input <- object@interactions %>% mutate(
    bin1 = `position 1` / object@binSize + 1,
    bin2 = `position 2` / object@binSize + 1
  )

  object@compartments <- tibble()
  object@concordances <- tibble()
  object@distances    <- tibble()
  object@centroids    <- tibble()

  for (chromosome in object@chromosomes) {

    message("Chromosome: ", chromosome)
    chromosomeInteractions <- input[input$chromosome == chromosome,]
    totalBins <- max(chromosomeInteractions$bin1, chromosomeInteractions$bin2)
    if (totalBins == -Inf) next

    positions <- (seq.int(totalBins) - 1) * object@binSize

    for (conditionId in unique(object@conditions)) {

      interactions <- matrix(nrow = 0, ncol = totalBins)
      replicates <- object@replicates[which(object@conditions == conditionId)]

      for (replicateId in replicates) {
        message("Replicate: ", conditionId, ".", replicateId)
        replicateInteractions <- chromosomeInteractions %>%
          filter(condition == conditionId) %>%
          filter(replicate == replicateId) %>%
          select(bin1, bin2, value)
        interactions <- rbind(
          interactions,
          sparseInteractionsToMatrix(replicateInteractions, totalBins)
        )
      }

      mustLink <- matrix(
        rep(0:(length(replicates) - 1), totalBins)*totalBins
        + rep(0:(totalBins - 1), each = length(replicates)),
        nrow = totalBins,
        byrow = TRUE
      )

      clusteringOutput <- constrainedClustering(
        interactions,
        mustLink,
        object@kMeansDelta,
        object@kMeansIterations,
        object@kMeansRestarts
      )

      clusters <- clusteringOutput[["clusters"]][0:totalBins] + 1
      centroids <- clusteringOutput[["centroids"]]

      min <- distanceRatio(
        centroids[[1]],
        centroids
      )

      max <- distanceRatio(
        centroids[[2]],
        centroids
      )

      concordances <- apply(interactions, 1, function(row) {
        2 * (distanceRatio(row, centroids) - min) / (max - min) - 1
      })

      distances <- apply(interactions, 1, function(row) {
        c(
          euclideanDistance(row, centroids[[1]]),
          euclideanDistance(row, centroids[[2]])
        )
      })

      object@compartments %<>% bind_rows(tibble(
        chromosome = chromosome,
        position = positions,
        condition = conditionId,
        value = clusters
      ))

      object@concordances %<>% bind_rows(tibble(
        chromosome = chromosome,
        position = rep(positions, length(replicates)),
        condition = conditionId,
        replicate = rep(replicates, each = ncol(interactions)),
        value = concordances
      ))

      object@distances %<>% bind_rows(tibble(
        chromosome = chromosome,
        position = rep(rep(positions, length(replicates)), 2),
        condition = conditionId,
        replicate = rep(rep(replicates, each = ncol(interactions)), 2),
        cluster = c(
          rep(1, length(replicates)*ncol(interactions)),
          rep(2, length(replicates)*ncol(interactions))
        ),
        value = c(distances[1,], distances[2,])
      ))

      object@centroids %<>% bind_rows(tibble(
        chromosome = chromosome,
        condition = conditionId,
        compartment = c(1, 2),
        centroid = centroids
      ))
    }
  }

  object@compartments %<>% mutate(
    chromosome = factor(chromosome),
    condition = factor(condition),
    value = factor(value)
  )

  object@concordances %<>% mutate(
    chromosome = factor(chromosome),
    condition = factor(condition),
    replicate = factor(replicate)
  )

  object@distances %<>% mutate(
    chromosome = factor(chromosome),
    condition = factor(condition),
    replicate = factor(replicate)
  )

  object@centroids %<>% mutate(
    chromosome = factor(chromosome),
    condition = factor(condition),
    compartment = factor(compartment)
  )

  return(object)
}
