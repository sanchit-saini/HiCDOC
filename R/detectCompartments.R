## - .euclideanDistance ------------------------------------------------------#
## --------------------------------------------------------------------------#
#' Compute the euclidean distance between two vectors.
#'
#' @param x  A vector.
#' @param y  A vector.
#'
#' @return A float number.
#' @keywords internal
#' @noRd
.euclideanDistance <- function(x, y) {
    sqrt(sum((x - y)^2))
}

## - .distanceRatio ----------------------------------------------------------#
## --------------------------------------------------------------------------#
#' Compute the log ratio of the distance of a position to each centroid.
#'
#' @param x          A vector of a genomic position.
#' @param centroids  A list of 2 vectors.
#' @param eps        A small epsilon to avoid log(0).
#'
#' @return A float number.
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

## - .tieCentroids -----------------------------------------------------------#
## --------------------------------------------------------------------------#
#' Assign correct cluster labels by comparing centroids across conditions.
#'
#' @param object A \code{HiCDOCDataSet} object.
#'
#' @return A \code{HiCDOCDataSet} object, with selfInteractionRatios, and with
#' corrected cluster labels in centroids, compartments, distances and
#' concordances.
#' @keywords internal
#' @noRd
.tieCentroids <- function(object) {
    totalConditions <- length(unique(object@conditions))
    referenceCondition <- object@conditions[[1]]

    clusters <-
        object@centroids %>%
        dplyr::filter(condition == referenceCondition) %>%
        tidyr::pivot_wider(
            names_from = compartment,
            values_from = centroid,
            names_prefix = "reference."
        ) %>%
        dplyr::select(-c(condition)) %>%
        tidyr::uncount(totalConditions) %>%
        dplyr::bind_cols(
            object@centroids %>%
            tidyr::pivot_wider(
                names_from = compartment,
                values_from = centroid,
                names_prefix = "centroid."
            ) %>%
            dplyr::select(-c(chromosome))
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
        dplyr::select(-c(centroid.1, centroid.2, reference.1, reference.2))

    object@compartments %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::mutate(
            compartment = dplyr::if_else(compartment == 1, cluster.1, cluster.2)
        ) %>%
        dplyr::select(-c(cluster.1, cluster.2))

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
        dplyr::select(-c(cluster.1, cluster.2, change))

    object@distances %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::mutate(
            compartment = dplyr::if_else(compartment == 1, cluster.1, cluster.2)
        ) %>%
        dplyr::select(-c(cluster.1, cluster.2))

    object@centroids %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::mutate(
            compartment = dplyr::if_else(compartment == 1, cluster.1, cluster.2)
        ) %>%
        dplyr::select(-c(cluster.1, cluster.2))

    return(object)
}

#' .constructLinkMatrix
#'
#' @param nbReplicates number of replicates under a condition
#' @param nbBins number of bins of a chromosome
#'
#' @return an matrix filled with 0 to use in constraned clustering
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

## - clusterize for 1 chromosome and 1 condition-----------------------------#
## --------------------------------------------------------------------------#
#' Segregate genomic positions into 2 clusters using constrained k-means.
#'
#' @param object  A \code{HiCDOCDataSet} object.
#' @param chromosomeName A character or numeric value.
#' Name or number of the chromosome
#' @param conditionName A character or numeric value.
#' Name or number of the condition
#'
#' @return A list of length 4:
#' * centroid (vector) for each cluster
#' * compartment (cluster number) for each genomic position
#' * distances to centroids (float) for each genomic position in
#' each replicate
#' * concordance (float between -1 and  1) for each genomic position in
#' each replicate
#' @md
#' @keywords internal
#' @noRd
.clusterizeChromosomeCondition <- function(
    object,
    chromosomeName,
    conditionName
) {

    totalBins <- object@totalBins[[chromosomeName]]
    if (totalBins == -Inf) return(NULL)
    positions <- seq.int(totalBins)
    # Substract filtered positions
    if (!is.null(object@weakBins[[chromosomeName]])) {
        positions <- positions[-object@weakBins[[chromosomeName]]]
        totalBins <- totalBins - length(object@weakBins[[chromosomeName]])
    }
    replicateNames <-
        object@replicates[which(object@conditions == conditionName)]
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
            chromosome = factor(chromosomeName, levels = object@chromosomes),
            bin = positions,
            condition = factor(
                conditionName,
                levels = unique(object@conditions)
            ),
            compartment = factor(clusters, levels = c(1, 2))
        )

    dfConcordances <-
        purrr::map_dfr(seq_len(length(replicateNames)), ~dfCompartments) %>%
        dplyr::mutate(
            replicate = rep(
                factor(replicateNames, levels = unique(object@replicates)),
                each = ncol(interactions)
            ),
            concordance = concordances
        ) %>%
        dplyr::select(
            chromosome, bin, condition, replicate, compartment, concordance
        )

    dfDistances <-
        purrr::map_dfr(seq_len(2), ~dfConcordances) %>%
        dplyr::select(-concordance) %>%
        dplyr::mutate(
            compartment = rep(
                factor(c(1, 2)),
                each = length(replicateNames) * ncol(interactions)
            ),
            distance = c(t(distances))
        )

    dfCentroids <-
        dplyr::tibble(
            chromosome = factor(chromosomeName, levels = object@chromosomes),
            condition = factor(
                conditionName,
                levels = unique(object@conditions)
            ),
            compartment = factor(c(1, 2)),
            centroid = centroids
        )

    return(list(
        "compartments" = dfCompartments,
        "concordances" = dfConcordances,
        "distances" = dfDistances,
        "centroids" = dfCentroids
    ))
}


## - clusterize -------------------------------------------------------------#
## --------------------------------------------------------------------------#
#' Segregate genomic positions into 2 clusters using constrained k-means.
#'
#' @param object  A \code{HiCDOCDataSet} object.
#' @param parallel Logical. Should parallel processing be used?
#'
#' @return A \code{HiCDOCDataSet} object, with:
#' - centroid (vector) for each cluster
#' - compartment (cluster number) for each genomic position
#' - distances to centroids (float) for each genomic position in each replicate
#' - concordance (float in [-1, 1]) for each genomic position in each replicate
#' @keywords internal
#' @noRd
.clusterize <- function(object, parallel = FALSE) {
    chromosomeNames <-
        rep(object@chromosomes, each = length(unique(object@conditions)))
    conditionNames <-
        rep(unique(object@conditions), length(object@chromosomes))

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
        BiocParallel::bpstart(BiocParallel::bpparam())
        # bpstart to address this issue (and ensure reproductibility):
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
        BiocParallel::bpstop(BiocParallel::bpparam())
        # Idem than bpstart
    }

    object@compartments <- purrr::map_dfr(result, "compartments")
    object@concordances <- purrr::map_dfr(result, "concordances")
    object@distances <- purrr::map_dfr(result, "distances")
    object@centroids <- purrr::map_dfr(result, "centroids")

    object <- .tieCentroids(object)
    return(object)
}

## - .computeSelfInteractionRatios ---------------------------------------------------------#
## --------------------------------------------------------------------------#
#' Compute the ratio of the reads which are on the diagonal vs those off-
#'   diagonal, for each bin, and each matrix.
#'
#' @param object an HiCDOCDataSet object
#' @param chromosomeName A character or numeric value.
#' Name or number of the chromosome
#' @param conditionName A character or numeric value.
#' Name or number of the condition
#' @param replicateName A character or numeric value.
#' Name or number of the replicate
#'
#' @return a \code{tibble}.
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
        dplyr::select(-c(diagonal, offDiagonal))

    return(selfInteractionRatios)
}

## - .predictCompartmentsAB --------------------------------------------------------------#
## --------------------------------------------------------------------------#
#' Use ratio between diagonal and off-diagonal interactions to determine which
#' clusters correspond to compartments A and B.
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param parallel Logical. Should parallel processing be used?
#'
#' @return A \code{HiCDOCDataSet} object, with selfInteractionRatios, and with A and B
#' labels replacing cluster numbers in centroids, compartments, distances and
#' concordances.
#' @keywords internal
#' @noRd
.predictCompartmentsAB <- function(
    object,
    parallel = FALSE
) {
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
        object@selfInteractionRatios <- do.call("rbind", object@selfInteractionRatios)
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
        dplyr::select(-c(ratio.1, ratio.2))

    object@compartments %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(
            compartment = factor(dplyr::if_else(compartment == A, "A", "B"))
        ) %>%
        dplyr::select(-c(A))

    object@concordances %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(change = dplyr::if_else(A == 1, 1, -1)) %>%
        dplyr::mutate(concordance = change * concordance) %>%
        dplyr::select(-c(A, change))

    object@distances %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(
            compartment = factor(dplyr::if_else(compartment == A, "A", "B"))
        ) %>%
        dplyr::select(-c(A))

    object@centroids %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(
            compartment = factor(dplyr::if_else(compartment == A, "A", "B"))
        ) %>%
        dplyr::select(-c(A))

    return(object)
}

## - .computePValues ---------------------------------------------------------#
## --------------------------------------------------------------------------#
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
#' @param object A \code{HiCDOCDataSet} object.
#'
#' @return A \code{HiCDOCDataSet} object, with differences and their p-values.
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
            dplyr::rename_with(function(old_colnames) {
                paste(old_colnames, 2, sep = ".")
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
            dplyr::rename_with(function(old_colnames) {
                paste(old_colnames, 2, sep = ".")
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
        dplyr::mutate(quantile = stats::ecdf(H0_value)(value)) %>%
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
                dplyr::if_else(compartment.1 == "A", "A->B", "B->A")
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
            order(gtools::mixedsort(chromosome)),
            bin,
            condition.1,
            condition.2
        )

    return(object)
}

## - detectCompartments -----------------------------------------------------#
## --------------------------------------------------------------------------#
#' Detect compartments for each genomic position in each condition, and compute
#' p-values for compartment differences between conditions.
#'
#' @param object  A \code{HiCDOCDataSet} object.
#' @param parallel Logical. Should parallel processing be used?
#' Default to FALSE.
#' @param kMeansDelta A numerical value. The maximum number of 2-means
#' iterations. If NULL, default to the first not NULL of
#' \code{object$kMeansDelta} and \code{HiCDOCDefaultParameters$kMeansDelta}.
#' @param kMeansIterations A numerical value. The stop criterion of
#' convergence of the 2-means method. If NULL, default to the first not
#' NULL of \code{object$kMeansIterations} and
#' \code{HiCDOCDefaultParameters$kMeansIterations}.
#' @param kMeansRestarts A numerical value. The maximum number of restarts
#' for the 2-means. If NULL, default to the first not NULL of
#' \code{object$kMeansRestarts} and
#' \code{HiCDOCDefaultParameters$kMeansRestarts}.
#'
#' @return A \code{HiCDOCDataSet} object, with centroids, compartments,
#' distances,concordances and differences.
#'
#' @examples
#' object <- HiCDOCExample()
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' @export
detectCompartments <- function(
    object,
    parallel = FALSE,
    kMeansDelta = NULL,
    kMeansIterations = NULL,
    kMeansRestarts = NULL
) {
    .validateSlots(
        object,
        slots = c(
            "chromosomes",
            "conditions",
            "totalBins",
            "binSize",
            "weakBins",
            "interactions",
            "parameters"
        )
    )
    # Parameters
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

    message("Clustering.")
    object <- .clusterize(object, parallel)
    message("Predicting compartments.")
    object <- .predictCompartmentsAB(object, parallel)
    message("Computing p-values.")
    object <- .computePValues(object)

    return(object)
}
