euclideanDistance <- function(x, y) {
  sqrt(sum((x - y)^2))
}

distanceRatio <- function(x, centroids, eps=1e-10) {
  return (log(
    (euclideanDistance(x, centroids[[1]]) + eps)
    / (euclideanDistance(x, centroids[[2]]) + eps)
  ))
}


clusterize <- function(object) {
  object@compartments <- tibble()
  object@concordances <- tibble()
  object@distances    <- tibble()
  object@centroids    <- tibble()
  progress            <- progress_estimated(length(object@chromosomes) * length(unique(object@conditions)))

  for (chromosomeId in object@chromosomes) {

    #message("Chromosome: ", chromosomeId)
    totalBins <- object@totalBins[[chromosomeId]]
    if (totalBins == -Inf) next

    positions <- (seq.int(totalBins) - 1) * object@binSize

    # Correct for filtered bins
    if (!is.null(object@weakBins[[chromosomeId]])) {
      positions <- positions[-object@weakBins[[chromosomeId]]]
      if (totalBins < length(object@weakBins[[chromosomeId]])) {
        message("Problem while filtering bins")
        message("Chr: ", chromosomeId)
        message("# bins: ", totalBins)
        message("max bins", max(object@interactions$position.1, object@interactions$position.2))
        message("bin size", object@binSize)
        message("# weak bins: ", length(object@weakBins[[chromosomeId]]))
        message("weak bins: ", paste(object@weakBins[[chromosomeId]]), " ")
      }
      totalBins <- totalBins - length(object@weakBins[[chromosomeId]])
    }

    for (conditionId in unique(object@conditions)) {

      interactions <- matrix(nrow = 0, ncol = totalBins)
      replicates <- object@replicates[which(object@conditions == conditionId)]

      for (replicateId in replicates) {
        #message("Replicate: ", conditionId, ".", replicateId)
        newInteractions <- sparseInteractionsToMatrix(
          object,
          chromosomeId,
          conditionId,
          replicateId,
          filter = TRUE
        )
        interactions <- rbind(interactions, newInteractions)
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
        chromosome = chromosomeId,
        position = positions,
        condition = conditionId,
        value = clusters
      ))

      object@concordances %<>% bind_rows(tibble(
        chromosome = chromosomeId,
        position = rep(positions, length(replicates)),
        condition = conditionId,
        replicate = rep(replicates, each = ncol(interactions)),
        value = concordances
      ))

      object@distances %<>% bind_rows(tibble(
        chromosome = chromosomeId,
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
        chromosome = chromosomeId,
        condition = conditionId,
        compartment = c(1, 2),
        centroid = centroids
      ))

      progress$tick()$print()
    }
  }
  progress$stop()

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

#' @export
detectCompartments <- function(object) {

  # Refiltering, since previous normalizations  may have introduced empty
  # rows/cols.
  message("Refiltering data...")
  #object  <- filterWeakPositions(object)
  object <- filterSmallChromosomes(object)

  message("Clustering...")
  object <- clusterize(object)
  message("Predicting compartments...")
  object <- predictABCompartments(object)
  message("Computing p-values...")
  object <- computePValues(object)
  message("...done.")

  return(object)
}
