#' @description
#' Computes the euclidean distance between two vectors.
#'
#' @param x
#' A vector.
#' @param y
#' A vector.
#'
#' @return
#' A float number.
#'
#' @keywords internal
#' @noRd
.euclideanDistance <- function(x, y) {
    sqrt(sum((x - y)^2))
}

#' @description
#' Computes the log ratio of the distance of a position to each centroid.
#'
#' @param x
#' The vector of a genomic position.
#' @param centroids
#' A list of two vectors.
#' @param eps
#' A small float number to avoid log(0).
#'
#' @return
#' A float number.
#'
#' @keywords internal
#' @noRd
.distanceRatio <- function(x, centroids, eps = 1e-10) {
    return(
        log(
            (.euclideanDistance(x, centroids[[1]]) + eps) /
            (.euclideanDistance(x, centroids[[2]]) + eps)
        )
    )
}

#' @description
#' Assigns correct cluster labels by comparing centroids across conditions.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}} with corrected cluster labels in compartments,
#' concordances, distances and centroids.
#'
#' @keywords internal
#' @noRd
.tieCentroids <- function(object) {

    selectedChromosomeNames <- names(base::Filter(function(x) !is.null(x),
                                                  object@validConditions))

    referenceConditionNames <-
        lapply(
            object@validConditions[selectedChromosomeNames],
            function(conditionNames) conditionNames[1]
            ) %>%
        unlist() %>%
        as.vector()

    referenceCentroids <-
        dplyr::tibble(
            chromosome = factor(
                selectedChromosomeNames,
                levels =
                    gtools::mixedsort(unique(object@chromosomes))),
            condition = factor(
                referenceConditionNames,
                levels = gtools::mixedsort(unique(object@conditions)))
        ) %>%
        dplyr::left_join(
            object@centroids,
            by = c("chromosome", "condition")
        )

    referenceCentroids <-
        referenceCentroids[!is.na(referenceCentroids$compartment),]

    clusters <-
        referenceCentroids %>%
        tidyr::pivot_wider(
            names_from = compartment,
            values_from = centroid,
            names_prefix = "reference."
        ) %>%
        dplyr::select(-condition) %>%
        dplyr::full_join(
            object@centroids %>%
            dplyr::filter(chromosome %in% selectedChromosomeNames) %>%
            tidyr::pivot_wider(
                names_from = compartment,
                values_from = centroid,
                names_prefix = "centroid."
            ),
            by = "chromosome"
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            cluster.1 = dplyr::if_else(
                .euclideanDistance(centroid.1, reference.1) *
                .euclideanDistance(centroid.2, reference.2) <
                .euclideanDistance(centroid.1, reference.2) *
                .euclideanDistance(centroid.2, reference.1),
                1,
                2
            )
        ) %>%
        dplyr::mutate(cluster.2 = dplyr::if_else(cluster.1 == 1, 2, 1)) %>%
        dplyr::select(-centroid.1, -centroid.2, -reference.1, -reference.2)

    object@compartments %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::mutate(
            compartment = dplyr::if_else(compartment == 1, cluster.1, cluster.2)
        ) %>%
        dplyr::select(-cluster.1, -cluster.2)

    object@concordances %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::mutate(
            change = dplyr::if_else(
                compartment == 1,
                dplyr::if_else(compartment == cluster.1, 1, -1),
                dplyr::if_else(compartment == cluster.2, 1, -1)
            )
        ) %>%
        dplyr::mutate(concordance = change * concordance) %>%
        dplyr::mutate(
            compartment = dplyr::if_else(compartment == 1, cluster.1, cluster.2)
        ) %>%
        dplyr::select(-cluster.1, -cluster.2, -change)

    object@distances %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::mutate(
            compartment = dplyr::if_else(compartment == 1, cluster.1, cluster.2)
        ) %>%
        dplyr::select(-cluster.1, -cluster.2) %>%
        dplyr::mutate(condition = factor(
            condition,
            levels = gtools::mixedsort(unique(object@conditions))))

    object@centroids %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::mutate(
            compartment = dplyr::if_else(compartment == 1, cluster.1, cluster.2)
        ) %>%
        dplyr::select(-cluster.1, -cluster.2)

    return(object)
}

#' @description
#' Constructs a link matrix of interaction rows to be clustered together.
#'
#' @param totalReplicates
#' The number of replicates.
#' @param totalBins
#' The number of bins.
#'
#' @return
#' A matrix, each row holding the row indices of interactions to be clustered
#' together.
#'
#' @keywords internal
#' @noRd
.constructLinkMatrix <- function(totalReplicates, totalBins) {
    return(
        matrix(
            rep(0:(totalReplicates - 1), totalBins) * totalBins +
            rep(0:(totalBins - 1), each = totalReplicates),
            nrow = totalBins,
            byrow = TRUE
        )
    )
}

#' @description
#' Segregates positions of a chromosome into two clusters using constrained
#' K-means.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' The name of a chromosome.
#' @param conditionName
#' The name of a condition.
#'
#' @return
#' A list of:
#' - The compartment (cluster number) of each position.
#' - The concordance (float) of each genomic position in each replicate.
#' - The distances to centroids (float) of each position in each replicate.
#' - The centroid (vector) of each cluster.
#'
#' @md
#' @keywords internal
#' @noRd
.clusterizeChromosomeCondition <- function(
    object,
    chromosomeName,
    conditionName
) {
    totalBins <- object@totalBins[[chromosomeName]]
    bins <- seq.int(totalBins)
    if (!is.null(object@weakBins[[chromosomeName]])) {
        bins <- bins[-object@weakBins[[chromosomeName]]]
        totalBins <- totalBins - length(object@weakBins[[chromosomeName]])
    }
    if (totalBins <= 1) return(NULL)

    replicateNames <-
        object@validReplicates[[chromosomeName]][
            which(object@validConditions[[chromosomeName]] == conditionName)
        ]
    if(length(replicateNames) == 0) return(NULL)

    interactions <-
        purrr::map(
            replicateNames,
            function(replicateName) {
                .sparseInteractionsToMatrix(
                    object,
                    chromosomeName,
                    conditionName,
                    replicateName,
                    filter = TRUE
                )
            }
        )
    interactions <- do.call("rbind", interactions)

    mustLink <- .constructLinkMatrix(length(replicateNames), totalBins)
    clusteringOutput <-
        constrainedClustering(
            interactions,
            mustLink,
            object@parameters$kMeansDelta,
            object@parameters$kMeansIterations,
            object@parameters$kMeansRestarts
        )
    clusters <- clusteringOutput[["clusters"]][0:totalBins] + 1
    centroids <- clusteringOutput[["centroids"]]

    min <- .distanceRatio(centroids[[1]], centroids)
    max <- .distanceRatio(centroids[[2]], centroids)

    concordances <-
        apply(interactions, 1, function(row) {
            2 * (.distanceRatio(row, centroids) - min) / (max - min) - 1
        })

    distances <-
        apply(interactions, 1, function(row) {
            c(
                .euclideanDistance(row, centroids[[1]]),
                .euclideanDistance(row, centroids[[2]])
            )
        })

    dfCompartments <-
        dplyr::tibble(
            chromosome = factor(
                chromosomeName,
                levels = gtools::mixedsort(unique(object@chromosomes))
            ),
            condition = factor(
                conditionName,
                levels = gtools::mixedsort(unique(object@conditions))
            ),
            bin = bins,
            compartment = factor(clusters, levels = c(1, 2))
        )

    dfConcordances <-
        purrr::map_dfr(seq_len(length(replicateNames)), ~dfCompartments) %>%
        dplyr::mutate(
            replicate = rep(
                factor(
                    replicateNames,
                    levels = gtools::mixedsort(unique(object@replicates))
                ),
                each = ncol(interactions)
            ),
            concordance = concordances
        ) %>%
        dplyr::select(
            chromosome, condition, replicate, bin, compartment, concordance
        )

    dfDistances <-
        purrr::map_dfr(seq_len(2), ~dfConcordances) %>%
        dplyr::select(-concordance) %>%
        dplyr::mutate(
            compartment = rep(
                factor(c(1, 2), levels = c(1, 2)),
                each = length(replicateNames) * ncol(interactions)
            ),
            distance = c(t(distances))
        )

    dfCentroids <-
        dplyr::tibble(
            chromosome = factor(
                chromosomeName,
                levels = gtools::mixedsort(unique(object@chromosomes))
            ),
            condition = factor(
                conditionName,
                levels = gtools::mixedsort(unique((object@conditions)))
            ),
            compartment = factor(c(1, 2), levels = c(1, 2)),
            centroid = centroids
        )

    return(list(
        "compartments" = dfCompartments,
        "concordances" = dfConcordances,
        "distances" = dfDistances,
        "centroids" = dfCentroids
    ))
}

#' @description
#' Runs the clustering to detect compartments in each chromosome and condition.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param parallel
#' Whether or not to parallelize the processing. Defaults to TRUE.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}} with compartments, concordances, distances
#' and centroids.
#'
#' @keywords internal
#' @noRd
.clusterize <- function(object, parallel = TRUE) {

    chromosomeNames <-
        rep(object@chromosomes, each = length(unique(object@conditions)))
    conditionNames <-
            rep(gtools::mixedsort(unique(object@conditions)),
                length(object@chromosomes))

    if (!parallel) {
        result <-
            pbapply::pbmapply(
                FUN = function(chromosomeName, conditionName) {
                    .clusterizeChromosomeCondition(
                        object,
                        chromosomeName,
                        conditionName
                    )
                },
                chromosomeNames,
                conditionNames,
                SIMPLIFY = FALSE
            )
    } else {
        reducedObjects <-
            purrr::map2(
                chromosomeNames,
                conditionNames,
                function(chromosomeName, conditionName) {
                    reduceHiCDOCDataSet(
                        object,
                        chromosomes = chromosomeName,
                        conditions = conditionName,
                        dropLevels = FALSE
                    )
                }
            )
        if(!is.null(BiocParallel::bpparam()$RNGseed))
            BiocParallel::bpstart(BiocParallel::bpparam())
        # bpstart to address this issue (and ensure reproducibility):
        # https://github.com/Bioconductor/BiocParallel/issues/122
        result <-
            BiocParallel::bpmapply(
                FUN = .clusterizeChromosomeCondition,
                reducedObjects,
                chromosomeNames,
                conditionNames,
                SIMPLIFY = FALSE,
                BPPARAM = BiocParallel::bpparam()
            )
        if(!is.null(BiocParallel::bpparam()$RNGseed))
            BiocParallel::bpstop(BiocParallel::bpparam())
        # Idem than bpstart
    }

    object@compartments <-
        purrr::map_dfr(result, "compartments") %>%
        dplyr::arrange(
            chromosome,
            condition,
            bin
        )

    object@concordances <-
        purrr::map_dfr(result, "concordances") %>%
        dplyr::arrange(
            chromosome,
            condition,
            replicate,
            bin
        )

    object@distances <-
        purrr::map_dfr(result, "distances") %>%
        dplyr::arrange(
            chromosome,
            condition,
            replicate,
            bin,
            compartment
        )

    object@centroids <-
        purrr::map_dfr(result, "centroids") %>%
        dplyr::arrange(
            chromosome,
            condition,
            compartment
        )

    object <- .tieCentroids(object)

    return(object)
}

#' @description
#' Computes the ratio of self interaction vs median of other interactions for
#' each position.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' The name of a chromosome.
#' @param conditionName
#' The name of a condition.
#' @param replicateName
#' The name of a replicate.
#'
#' @return
#' A tibble.
#'
#' @keywords internal
#' @noRd
.computeSelfInteractionRatios <- function(
    object,
    chromosomeName,
    conditionName,
    replicateName
) {

    interactions <-
        object@interactions %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::filter(condition == conditionName) %>%
        dplyr::filter(replicate == replicateName)

    diagonal <-
        interactions %>%
        dplyr::filter(bin.1 == bin.2) %>%
        dplyr::select(-bin.2) %>%
        dplyr::rename(bin = bin.1) %>%
        dplyr::rename(diagonal = interaction)

    offDiagonal <-
        interactions %>%
        dplyr::filter(bin.1 != bin.2) %>%
        tidyr::pivot_longer(
            cols = tidyr::starts_with("bin"),
            names_to = "namepos",
            values_to = "bin"
        ) %>%
        dplyr::select(-namepos) %>%
        dplyr::group_by(chromosome, condition, replicate, bin) %>%
        dplyr::summarise(offDiagonal = stats::median(interaction)) %>%
        dplyr::ungroup()

    selfInteractionRatios <-
        diagonal %>%
        dplyr::full_join(
            offDiagonal,
            by = c(
                "chromosome",
                "condition",
                "replicate",
                "bin"
            )
        ) %>%
        tidyr::replace_na(list(diagonal = 0, offDiagonal = 0)) %>%
        dplyr::mutate(ratio = diagonal - offDiagonal) %>%
        dplyr::select(-diagonal, -offDiagonal)

    return(selfInteractionRatios)
}

#' @description
#' Uses ratio between self interactions and other interactions to determine
#' which clusters correspond to compartments A and B.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param parallel
#' Whether or not to parallelize the processing. Defaults to TRUE.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}, with selfInteractionRatios, and with A and B
#' labels replacing cluster numbers in compartments, concordances, distances,
#' and centroids.
#'
#' @keywords internal
#' @noRd
.predictCompartmentsAB <- function(object, parallel = TRUE) {

    chromosomeNames <- rep(object@chromosomes, each = length(object@replicates))
    conditionNames <- rep(object@conditions, length(object@chromosomes))
    replicateNames <- rep(object@replicates, length(object@chromosomes))

    if (!parallel) {
        ratios <-
            pbapply::pbmapply(
                function(chromosomeName, conditionName, replicateName) {
                    .computeSelfInteractionRatios(
                        object,
                        chromosomeName,
                        conditionName,
                        replicateName
                    )
                },
                chromosomeNames,
                conditionNames,
                replicateNames,
                SIMPLIFY = FALSE
            )
        object@selfInteractionRatios <- do.call("rbind", ratios)
    } else {
        reducedObjects <-
            purrr::pmap(
                list(chromosomeNames, conditionNames, replicateNames),
                function(chromosomeName, conditionName, replicateName) {
                    reduceHiCDOCDataSet(
                        object,
                        chromosomes = chromosomeName,
                        conditions = conditionName,
                        replicates = replicateName,
                        dropLevels = FALSE
                    )
                }
            )
        object@selfInteractionRatios <-
            BiocParallel::bpmapply(
                FUN = .computeSelfInteractionRatios,
                reducedObjects,
                chromosomeNames,
                conditionNames,
                replicateNames,
                SIMPLIFY = FALSE,
                BPPARAM = BiocParallel::bpparam()
            )
        object@selfInteractionRatios <-
            do.call("rbind", object@selfInteractionRatios)
    }

    compartments <-
        object@compartments %>%
        dplyr::left_join(
            object@selfInteractionRatios,
            by = c("chromosome", "condition", "bin")
        ) %>%
        dplyr::group_by(chromosome, compartment) %>%
        dplyr::summarize(ratio = stats::median(ratio)) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(
            names_from = compartment,
            values_from = ratio,
            values_fill = list(ratio = 0),
            names_prefix = "ratio."
        ) %>%
        dplyr::mutate(A = dplyr::if_else(ratio.1 >= ratio.2, 2, 1)) %>%
        dplyr::select(-ratio.1, -ratio.2)

    object@compartments %<>%
        dplyr::left_join(compartments, by = "chromosome") %>%
        dplyr::mutate(
            compartment = factor(
                dplyr::if_else(compartment == A, "A", "B"),
                levels = c("A", "B")
            )
        ) %>%
        dplyr::select(-A)

    object@concordances %<>%
        dplyr::left_join(compartments, by = "chromosome") %>%
        dplyr::mutate(change = dplyr::if_else(A == 1, 1, -1)) %>%
        dplyr::mutate(concordance = change * concordance) %>%
        dplyr::mutate(
            compartment = factor(
                dplyr::if_else(compartment == A, "A", "B"),
                levels = c("A", "B")
            )
        ) %>%
        dplyr::select(-A, -change)

    object@distances %<>%
        dplyr::left_join(compartments, by = "chromosome") %>%
        dplyr::mutate(
            compartment = factor(
                dplyr::if_else(compartment == A, "A", "B"),
                levels = c("A", "B")
            )
        ) %>%
        dplyr::select(-A)

    object@centroids %<>%
        dplyr::left_join(compartments, by = "chromosome") %>%
        dplyr::mutate(
            compartment = factor(
                dplyr::if_else(compartment == A, "A", "B"),
                levels = c("A", "B")
            )
        ) %>%
        dplyr::select(-A)

    return(object)
}

#' @description
#' Computes p-values for genomic positions whose assigned compartment switches
#' between two conditions.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}, with differences and their p-values.
#'
#' @keywords internal
#' @noRd
.computePValues <- function(object) {

    # Compute median of differences between pairs of concordances
    # N.b. median of differences != difference of medians
    totalReplicates <- length(object@replicates)

    concordanceDifferences <-
        object@concordances %>%
        dplyr::arrange(chromosome, bin, condition, replicate)

    concordanceDifferences %<>%
        # Duplicate each row "totalReplicates" times
        tidyr::uncount(totalReplicates) %>%
        cbind(
            # Duplicate the table "totalReplicates" times
            concordanceDifferences[rep(
                seq_len(nrow(concordanceDifferences)),
                totalReplicates
            ), ] %>%
            # Arrange by chromosome and position
            # to join with the duplicated rows
            dplyr::arrange(chromosome, bin) %>%
            dplyr::rename_with(function(column) {
                paste(column, 2, sep = ".")
            })
        ) %>%
        dplyr::select(-chromosome.2, -bin.2) %>%
        dplyr::rename(condition.1 = condition, concordance.1 = concordance) %>%
        dplyr::filter(as.numeric(condition.1) < as.numeric(condition.2)) %>%
        dplyr::group_by(chromosome, bin, condition.1, condition.2) %>%
        dplyr::summarise(
            value = stats::median(abs(concordance.1 - concordance.2))
        ) %>%
        dplyr::ungroup()

    # Format compartments per pair of conditions
    totalConditions <- length(unique(object@conditions))

    compartmentComparisons <-
        object@compartments %>%
        dplyr::arrange(chromosome, bin, condition)

    compartmentComparisons %<>%
        # Duplicate each row "totalConditions" times
        tidyr::uncount(totalConditions) %>%
        cbind(
            # Duplicate the table "totalConditions" times
            compartmentComparisons[rep(
                seq_len(nrow(compartmentComparisons)),
                totalConditions
            ), ] %>%
            # Arrange by chromosome and position
            # to join with the duplicated rows
            dplyr::arrange(chromosome, bin) %>%
            dplyr::rename_with(function(column) {
                paste(column, 2, sep = ".")
            })
        ) %>%
        dplyr::select(-chromosome.2, -bin.2) %>%
        dplyr::rename(
            condition.1 = condition,
            compartment.1 = compartment
        ) %>%
        dplyr::filter(as.numeric(condition.1) < as.numeric(condition.2))

    # Compute p-values for switching positions
    # P-values for a condition pair computed from the whole genome distribution
    object@differences <-
        compartmentComparisons %>%
        dplyr::left_join(
            concordanceDifferences,
            by = c("chromosome", "bin", "condition.1", "condition.2")
        ) %>%
        dplyr::mutate(
            H0_value =
                dplyr::if_else(compartment.1 == compartment.2, value, NA_real_)
        ) %>%
        dplyr::group_by(condition.1, condition.2) %>%
        dplyr::mutate(quantile = ifelse(value > 0,
                                        stats::ecdf(H0_value)(value),
                                        NA)) %>%
        dplyr::filter(compartment.1 != compartment.2) %>%
        dplyr::mutate(pvalue = 1 - quantile) %>%
        dplyr::mutate(pvalue = dplyr::if_else(pvalue < 0, 0, pvalue)) %>%
        dplyr::mutate(pvalue = dplyr::if_else(pvalue > 1, 1, pvalue)) %>%
        dplyr::mutate(
            pvalue.adjusted = stats::p.adjust(pvalue, method = "BH")
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            direction = factor(
                dplyr::if_else(compartment.1 == "A", "A->B", "B->A"),
                levels = c("A->B", "B->A")
            )
        ) %>%
        dplyr::select(
            chromosome,
            bin,
            condition.1,
            condition.2,
            pvalue,
            pvalue.adjusted,
            direction
        ) %>%
        dplyr::arrange(
            chromosome,
            bin,
            condition.1,
            condition.2
        )

    return(object)
}

#' @title
#' A and B compartments detection and differences across conditions.
#'
#' @description
#' Detects compartments for each genomic position in each condition, and
#' computes p-values for compartment differences between conditions.
#'
#' @details
#' \subsection{Genomic positions clustering}{
#' To clusterize genomic positions, the algorithm follows these steps:
#'     \enumerate{
#'         \item{
#'             For each chromosome and condition, get the interaction vectors of
#'             each genomic position. Each genomic position can have multiple
#'             interaction vectors, corresponding to the multiple replicates in
#'             that condition.
#'         }
#'         \item{
#'             For each chromosome and condition, use constrained K-means to
#'             clusterize the interaction vectors, forcing replicate interaction
#'             vectors into the same cluster. The euclidean distance between
#'             interaction vectors determines their similarity.
#'         }
#'         \item{
#'             For each interaction vector, compute its concordance, which is
#'             the confidence in its assigned cluster. Mathematically, it is the
#'             log ratio of its distance to each centroid, normalized by the
#'             distance between both centroids, and min-maxed to a [-1,1]
#'             interval.
#'         }
#'         \item{
#'             For each chromosome, compute the distance between all centroids
#'             and the centroids of the first condition. The cross-condition
#'             clusters whose centroids are closest are given the same cluster
#'             label. This results in two clusters per chromosome, spanning all
#'             conditions.
#'         }
#'     }
#' }
#' \subsection{A/B compartments prediction}{
#' To match each cluster with an A or B compartment, the algorithm follows these
#' steps:
#'     \enumerate{
#'         \item{
#'             For each genomic position, compute its self interaction ratio,
#'             which is the difference between its self interaction and the
#'             median of its other interactions.
#'         }
#'         \item{
#'             For each chromosome, for each cluster, get the median self
#'             interaction ratio of the genomic positions in that cluster.
#'         }
#'         \item{
#'             For each chromosome, the cluster with the smallest median self
#'             interaction ratio is matched with compartment A, and the cluster
#'             with the greatest median self interaction ratio is matched with
#'             compartment B. Compartment A being open, there are more overall
#'             interactions between distant genomic positions, so it is assumed
#'             that the difference between self interactions and other
#'             interactions is lower than in compartment B.
#'         }
#'     }
#' }
#' \subsection{Significant differences detection}{
#' To find significant compartment differences across conditions, and compute
#' their p-values, the algorithm follows three steps:
#'     \enumerate{
#'         \item{
#'             For each pair of replicates in different conditions, for each
#'             genomic position, compute the absolute difference between its
#'             concordances.
#'         }
#'         \item{
#'             For each pair of conditions, for each genomic position, compute
#'             the median of its concordance differences.
#'         }
#'         \item{
#'             For each pair of conditions, for each genomic position whose
#'             assigned compartment switches, rank its median against the
#'             empirical cumulative distribution of medians of all non-switching
#'             positions in that condition pair. Adjust the resulting p-value
#'             with the Benjamini–Hochberg procedure.
#'         }
#'     }
#' }
#' \subsection{Parallel processing}{
#' The parallel version of detectCompartments uses the
#' \code{\link{BiocParallel}{bpmapply}} function. Before to call the
#' function in parallel you should specify the parallel parameters such as:
#'     \itemize{
#'         \item{On Linux:
#'
#'              \code{multiParam <- BiocParallel::MulticoreParam(workers = 10)}
#'              \code{BiocParallel::register(multiParam, default = TRUE)}
#'          }
#'          \item{On Windows:
#'
#'              \code{multiParam <- BiocParallel::SnowParam(workers = 10)}
#'              \code{BiocParallel::register(multiParam, default = TRUE)}
#'         }
#'     }
#' }
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @param parallel
#' Whether or not to parallelize the processing. Defaults to TRUE.
#' See 'Details'.
#'
#' @param kMeansDelta
#' The convergence stop criterion for the clustering. When the centroids'
#' distances between two iterations is lower than this value, the clustering
#' stops. Defaults to \code{object$kMeansDelta} which is originally set to
#' \code{defaultHiCDOCParameters$kMeansDelta} = 0.0001.
#'
#' @param kMeansIterations
#' The maximum number of iterations during clustering. Defaults to
#' \code{object$kMeansIterations} which is originally set to
#' \code{defaultHiCDOCParameters$kMeansIterations} = 50.
#'
#' @param kMeansRestarts
#' The amount of times the clustering is restarted. For each restart, the
#' clustering iterates until convergence or reaching the maximum number of
#' iterations. The clustering that minimizes inner-cluster variance is selected.
#' Defaults to \code{object$kMeansRestarts} which is originally set to
#' \code{defaultHiCDOCParameters$kMeansRestarts} = 20.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}, with compartments, concordances, distances,
#' centroids, and differences.
#'
#' @examples
#' data(HiCDOCDataSetExample)
#' # Run all steps f=of filter and normalization
#' object <- filterSmallChromosomes(HiCDOCDataSetExample)
#' object <- filterSparseReplicates(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#'
#' @usage
#' detectCompartments(
#'     object,
#'     parallel = TRUE,
#'     kMeansDelta = NULL,
#'     kMeansIterations = NULL,
#'     kMeansRestarts = NULL
#' )
#'
#' @export
detectCompartments <- function(
    object,
    parallel = TRUE,
    kMeansDelta = NULL,
    kMeansIterations = NULL,
    kMeansRestarts = NULL
) {

    .validateSlots(
        object,
        slots = c(
            "interactions",
            "chromosomes",
            "conditions",
            "replicates",
            "totalBins",
            "resolution",
            "weakBins",
            "parameters"
        )
    )
    if (is.null(object@validConditions) | is.null(object@validReplicates)) {
        valid <- lapply(object@chromosomes,
                        FUN = function(x)
                            object@interactions %>%
                            dplyr::filter(chromosome == x &
                                              interaction > 0) %>%
                            dplyr::select(condition, replicate) %>%
                            unique())
        object@validConditions <-
            lapply(valid, function(x) x$condition)
        object@validReplicates <-
            lapply(valid, function(x) x$replicate)
        names(object@validConditions) <- as.vector(object@chromosomes)
        names(object@validReplicates) <- as.vector(object@chromosomes)
    }
    if (!is.null(kMeansDelta)) {
        object@parameters$kMeansDelta <- kMeansDelta
    }
    if (!is.null(kMeansIterations)) {
        object@parameters$kMeansIterations <- kMeansIterations
    }
    if (!is.null(kMeansRestarts)) {
        object@parameters$kMeansRestarts <- kMeansRestarts
    }
    object@parameters <- .validateParameters(object@parameters)

    message("Clustering genomic positions.")
    object <- .clusterize(object, parallel)
    message("Predicting A/B compartments.")
    object <- .predictCompartmentsAB(object, parallel)
    message("Detecting significant differences.")
    object <- .computePValues(object)

    return(object)
}
