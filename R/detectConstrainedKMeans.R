euclideanDistance <- function(p1, p2) {
   sqrt(sum((p1 - p2)^2))
}

concordanceFunction <- function(x, p1, p2) {
  eps <- 1e-10
  log2((euclideanDistance(x, p1) + eps) / (euclideanDistance(x, p2) + eps))
}


detectConstrainedKMeans <- function(object) {
  input <- object@interactionMatrix %>%
    mutate(bin1 = `position 1` / object@binSize + 1) %>%
    mutate(bin2 = `position 2` / object@binSize + 1)

  mustLinkSeed <- list(rep1 = which(object@conditions == 1) - 1,
                       rep2 = which(object@conditions == 2) - 1)

  object@distances    <- tibble()
  object@compartments <- tibble()
  object@concordances <- tibble()

  for (chr in object@chromosomes) {
    inputChromosome <- filter(input, chromosome == chr)
    positions <- unique(inputChromosome$`position 1`)
    n <- max(inputChromosome$bin1, inputChromosome$bin2)
    bigmat <- matrix(nrow = 0, ncol = n)
    message(paste0("Chromosome ", chr, ", of dim. ", n))
    for (rep in object@replicates) {
      message(paste0("  Replicate ", rep))
      inputReplicate <- inputChromosome %>%
        filter(replicate == rep) %>%
        select(-c(chromosome, replicate, `position 1`, `position 2`))
      mat <- acast(inputReplicate, bin1 ~ bin2, fill = 0, value.var = "value")
      mat <- mat + t(mat) - diag(diag(mat))
      if (!isSymmetric(mat)) {
        stop(paste0("Matrix ", chr, "/", rep, " is not symmetric"))
      }
      bigmat <- rbind(bigmat, mat)
    }
    mustLink <- do.call("rbind",
                        lapply(mustLinkSeed,
                               function(x) {
                                 matrix(rep(x, n)*n +
                                          rep(0:(n - 1), each = length(x)),
                                        nrow = n, byrow = TRUE) } ))
    #mustLink <- do.call("rbind", lapply(mustLinkSeed, function(x) { matrix(rep(x, n)*n + rep(1:n, each=length(x)), nrow=n, byrow=TRUE) } ))
    #clusters <- ckmeans(bigmat, 2, mustLink, c())
    clusters <- constrainedClustering(bigmat, mustLink, 0.001, 10) + 1
    centroids <- lapply(1:2,
                        function(i) {
                          colMeans(bigmat[ which(clusters == i),  ])
                        }
                       )
    object@compartments %<>% bind_rows(tibble(
                                  chromosome = chr,
                                  position   = positions,
                                  condition  = 1,
                                  value = head(clusters, n))) %>%
                      bind_rows(tibble(
                                  chromosome = chr,
                                  position   = positions,
                                  condition  = 2,
                                  value = tail(clusters, n)))

    minMaxDistances <- lapply(centroids,
                              function(x) {
                                concordanceFunction(x, centroids[[1]], centroids[[2]])
                              }
                             )

    object@concordances <- bind_rows(object@concordances,
                              tibble(chromosome = chr,
                                     position = rep(positions,
                                                    object@nReplicates),
                                     replicate = rep(object@replicates,
                                                     each = ncol(bigmat)),
                                     value =
                                       apply(bigmat, 1, function(x) {
                                           (2 * concordanceFunction(x,
                                                                    centroids[[1]],
                                                                    centroids[[2]]) -
                                             minMaxDistances[[1]]) /
                                             (minMaxDistances[[2]] - minMaxDistances[[1]]) - 0.5
                                       })
                             ))
    object@distances %<>%
                 bind_rows(tibble(chromosome = chr,
                      position = rep(positions, object@nReplicates),
                      cluster = 1,
                      replicate = rep(object@replicates, each = ncol(bigmat)),
                      distance = apply(bigmat, 1, function(x) {
                                      euclideanDistance(x, centroids[[1]])}))) %>%
                 bind_rows(tibble(chromosome = chr,
                      position = rep(positions, object@nReplicates),
                      cluster = 2,
                      replicate = rep(object@replicates, each = ncol(bigmat)),
                      distance = apply(bigmat, 1, function(x) {
                                      euclideanDistance(x, centroids[[2]])})))
  }
  return(object)
}
