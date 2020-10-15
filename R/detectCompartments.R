##- euclideanDistance --------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Compute the euclidean distance between two vectors.
#'
#' @param x  A vector.
#' @param y  A vector.
#'
#' @return A float number.
euclideanDistance <- function(x, y) {
    sqrt(sum((x - y) ^ 2))
}

##- distanceRatio ------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Compute the log ratio of the distance of a position to each centroid.
#'
#' @param x          A vector of a genomic position.
#' @param centroids  A list of 2 vectors.
#' @param eps        A small epsilon to avoid log(0).
#'
#' @return A float number.
distanceRatio <- function(x, centroids, eps = 1e-10) {
    return (log(
        (euclideanDistance(x, centroids[[1]]) + eps)
        / (euclideanDistance(x, centroids[[2]]) + eps)
    ))
}

##- clusterize ---------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Segregate genomic positions into 2 clusters using constrained k-means.
#'
#' @param object  A \code{HiCDOCExp} object.
#'
#' @return A \code{HiCDOCExp} object, with:
#' - centroid (vector) for each cluster
#' - compartment (cluster number) for each genomic position
#' - distances to centroids (float) for each genomic position in each replicate
#' - concordance (float in [-1, 1]) for each genomic position in each replicate
clusterize <- function(object) {
    object@compartments <- tibble()
    object@concordances <- tibble()
    object@distances    <- tibble()
    object@centroids    <- tibble()

    progress <- progress_estimated(length(object@chromosomes) * length(unique(object@conditions)))

    for (chromosomeId in object@chromosomes) {
        totalBins <- object@totalBins[[chromosomeId]]
        if (totalBins == -Inf) next

        positions <- (seq.int(totalBins) - 1) * object@binSize

        # Correct for filtered bins
        if (!is.null(object@weakBins[[chromosomeId]])) {
            positions <- positions[-object@weakBins[[chromosomeId]]]
            totalBins <- totalBins - length(object@weakBins[[chromosomeId]])
        }

        for (conditionId in unique(object@conditions)) {
            interactions <- matrix(nrow = 0, ncol = totalBins)
            replicates <- object@replicates[which(object@conditions == conditionId)]

            for (replicateId in replicates) {
                replicateInteractions <- sparseInteractionsToMatrix(
                    object,
                    chromosomeId,
                    conditionId,
                    replicateId,
                    filter = TRUE
                )
                interactions <- rbind(interactions, replicateInteractions)
            }

            mustLink <- matrix(
                rep(0:(length(replicates) - 1), totalBins)
                * totalBins
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

            min <- distanceRatio(centroids[[1]], centroids)
            max <- distanceRatio(centroids[[2]], centroids)

            concordances <- apply(interactions, 1, function(row) {
                2 * (distanceRatio(row, centroids) - min) / (max - min) - 1
            })

            distances <- apply(interactions, 1, function(row) {
                c(
                    euclideanDistance(row, centroids[[1]]),
                    euclideanDistance(row, centroids[[2]])
                )
            })

            object@compartments %<>% bind_rows(
                tibble(
                    chromosome = chromosomeId,
                    position = positions,
                    condition = conditionId,
                    compartment = clusters
                )
            )

            object@concordances %<>% bind_rows(
                tibble(
                    chromosome = chromosomeId,
                    position = rep(positions, length(replicates)),
                    condition = conditionId,
                    replicate = rep(replicates, each = ncol(interactions)),
                    compartment = rep(clusters, length(replicates)),
                    concordance = concordances
                )
            )

            object@distances %<>% bind_rows(
                tibble(
                    chromosome = chromosomeId,
                    position = rep(rep(positions, length(replicates)), 2),
                    condition = conditionId,
                    replicate = rep(rep(replicates, each = ncol(interactions)), 2),
                    compartment = c(
                        rep(1, length(replicates) * ncol(interactions)),
                        rep(2, length(replicates) * ncol(interactions))
                    ),
                    distance = c(t(distances))
                )
            )

            object@centroids %<>% bind_rows(
                tibble(
                    chromosome = chromosomeId,
                    condition = conditionId,
                    compartment = c(1, 2),
                    centroid = centroids
                )
            )

            progress$tick()$print()
        }
    }

    progress$stop()

    object@compartments %<>% dplyr::mutate(
        chromosome = factor(chromosome),
        condition = factor(condition),
        compartment = factor(compartment)
    )

    object@concordances %<>% dplyr::mutate(
        chromosome = factor(chromosome),
        condition = factor(condition),
        replicate = factor(replicate)
    )

    object@distances %<>% dplyr::mutate(
        chromosome = factor(chromosome),
        condition = factor(condition),
        replicate = factor(replicate)
    )

    object@centroids %<>% dplyr::mutate(
        chromosome = factor(chromosome),
        condition = factor(condition),
        compartment = factor(compartment)
    )

    return (tieCentroids(object))
}

##- tieCentroids -------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Assign correct cluster labels by comparing centroids across conditions.
#'
#' @param object A \code{HiCDOCExp} object.
#'
#' @return A \code{HiCDOCExp} object, with diagonalRatios, and with corrected
#' cluster labels in centroids, compartments, distances and concordances.
tieCentroids <- function(object) {
    totalConditions <- length(unique(object@conditions))
    referenceCondition <- object@conditions[[1]]

    clusters <- object@centroids %>%
        dplyr::filter(condition == referenceCondition) %>%
        dplyr::group_by(chromosome, condition) %>%
        spread(key = compartment, value = centroid) %>%
        dplyr::ungroup() %>%
        uncount(totalConditions) %>%
        dplyr::rename(reference.1 = `1`) %>%
        dplyr::rename(reference.2 = `2`) %>%
        dplyr::select(-c(condition)) %>%
        dplyr::bind_cols(
            object@centroids %>%
                dplyr::group_by(chromosome, condition) %>%
                spread(key = compartment, value = centroid) %>%
                dplyr::ungroup() %>%
                dplyr::rename(centroid.1 = `1`) %>%
                dplyr::rename(centroid.2 = `2`) %>%
                dplyr::select(-c(chromosome))
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(`1` = if_else(
            euclideanDistance(centroid.1, reference.1) *
            euclideanDistance(centroid.2, reference.2) <
            euclideanDistance(centroid.1, reference.2) *
            euclideanDistance(centroid.2, reference.1), 1, 2
        )) %>%
        dplyr::mutate(`2` = if_else(`1` == 1, 2, 1)) %>%
        dplyr::select(-c(centroid.1, centroid.2, reference.1, reference.2))

    object@compartments %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::mutate(compartment = if_else(compartment == 1, `1`, `2`)) %>%
        dplyr::select(-c(`1`, `2`))

    object@concordances %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(change = if_else(
            compartment == 1,
            if_else(compartment == `1`, 1, -1),
            if_else(compartment == `2`, 1, -1)
        )) %>%
        dplyr::mutate(concordance = change * concordance) %>%
        dplyr::mutate(compartment = if_else(compartment == 1, `1`, `2`)) %>%
        dplyr::select(-c(`1`, `2`, change))

    object@distances %<>%
        dplyr::left_join(clusters, by=c("chromosome", "condition")) %>%
        dplyr::mutate(compartment = if_else(compartment == 1, `1`, `2`)) %>%
        dplyr::select(-c(`1`, `2`))

    object@centroids %<>%
        dplyr::left_join(clusters, by=c("chromosome", "condition")) %>%
        dplyr::mutate(compartment = if_else(compartment == 1, `1`, `2`)) %>%
        dplyr::select(-c(`1`, `2`))

    return (object)
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
            dplyr::filter(position.1 == position.2) %>%
            dplyr::select(-position.2) %>%
            dplyr::rename(position = position.1) %>%
            dplyr::rename(diagonal = value)

        offDiagonal <- fullInteractions %>%
            dplyr::filter(position.1 != position.2) %>%
            dplyr::select(-position.2) %>%
            dplyr::rename(position = position.1) %>%
            dplyr::group_by(position) %>%
            dplyr::mutate(offDiagonal = median(value)) %>%
            dplyr::ungroup() %>%
            dplyr::select(-value) %>%
            dplyr::distinct()

        diagonalRatios <- diagonal %>%
            dplyr::left_join(
                offDiagonal,
                by = c("chromosome", "condition", "replicate", "position")
            ) %>%
            dplyr::mutate(value = diagonal - offDiagonal) %>%
            dplyr::select(-c(diagonal, offDiagonal))

        return (diagonalRatios)
    }))

    compartments <- object@compartments %>%
        dplyr::left_join(
            object@diagonalRatios,
            by = c("chromosome", "condition", "position")
        ) %>%
        dplyr::group_by(chromosome, compartment) %>%
        dplyr::summarize(value = median(value)) %>%
        dplyr::ungroup() %>%
        tidyr::spread(
            key = compartment,
            value = value,
            fill = 0
        ) %>%
        dplyr::mutate(A = if_else(`1` >= `2`, 2, 1)) %>%
        dplyr::select(-c(`1`, `2`))

    object@compartments %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(compartment = factor(if_else(compartment == A, "A", "B"))) %>%
        dplyr::select(-c(A))

    object@concordances %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(change = if_else(A == 1, 1, -1)) %>%
        dplyr::mutate(concordance = change * concordance) %>%
        dplyr::select(-c(A, change))

    object@distances %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(compartment = factor(if_else(compartment == A, "A", "B"))) %>%
        dplyr::select(-c(A))

    object@centroids %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(compartment = factor(if_else(compartment == A, "A", "B"))) %>%
        dplyr::select(-c(A))

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
        dplyr::arrange(chromosome, position, condition, replicate)

    concordanceDifferences %<>%
        # Duplicate each row "totalReplicates" times
        uncount(totalReplicates) %>%
        cbind(
            # Duplicate the table "totalReplicates" times
            concordanceDifferences[rep(seq_len(nrow(concordanceDifferences)),
                                       totalReplicates), ] %>%
                # Arrange by chromosome and position to join with the duplicated rows
                dplyr::arrange(chromosome, position) %>%
                dplyr::rename_with(function(old_colnames)
                    paste(old_colnames, 2, sep = "."))
        ) %>%
        dplyr::select(-chromosome.2,-position.2) %>%
        dplyr::rename(condition.1 = condition, concordance.1 = concordance) %>%
        dplyr::filter(as.numeric(condition.1) < as.numeric(condition.2)) %>%
        dplyr::group_by(chromosome, position, condition.1, condition.2) %>%
        dplyr::summarise(value = median(abs(concordance.1 - concordance.2))) %>%
        dplyr::ungroup()

    # Format compartments per pair of conditions
    totalConditions = length(unique(object@conditions))

    compartmentComparisons <- object@compartments %>%
        dplyr::arrange(chromosome, position, condition)

    compartmentComparisons %<>%
        # Duplicate each row "totalConditions" times
        uncount(totalConditions) %>%
        cbind(
            # Duplicate the table "totalConditions" times
            compartmentComparisons[rep(seq_len(nrow(compartmentComparisons)),
                                       totalConditions), ] %>%
                # Arrange by chromosome and position to join with the duplicated rows
                dplyr::arrange(chromosome, position) %>%
                dplyr::rename_with(function(old_colnames)
                    paste(old_colnames, 2, sep = "."))
        ) %>%
        dplyr::select(-chromosome.2,-position.2) %>%
        dplyr::rename(
            condition.1 = condition,
            compartment.1 = compartment
        ) %>%
        dplyr::filter(as.numeric(condition.1) < as.numeric(condition.2))

    # Compute p-values for switching positions
    # P-values for a condition pair computed from the whole genome distribution
    object@differences <- compartmentComparisons %>%
        dplyr::left_join(
            concordanceDifferences,
            by = c("chromosome", "position", "condition.1", "condition.2")
        ) %>%
        dplyr::mutate(H0_value = if_else(compartment.1 == compartment.2, value, NA_real_)) %>%
        dplyr::group_by(condition.1, condition.2) %>%
        dplyr::mutate(quantile = ecdf(H0_value)(value)) %>%
        dplyr::filter(compartment.1 != compartment.2) %>%
        dplyr::mutate(pvalue = 1 - quantile) %>%
        dplyr::mutate(pvalue = if_else(pvalue < 0, 0, pvalue)) %>%
        dplyr::mutate(pvalue = if_else(pvalue > 1, 1, pvalue)) %>%
        dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) %>%
        dplyr::ungroup() %>%
        dplyr::rename(start = position) %>%
        dplyr::mutate(end = start + object@binSize) %>%
        dplyr::mutate(direction = factor(if_else(compartment.1 == "A", "A->B", "B->A"))) %>%
        dplyr::select(
            chromosome,
            start,
            end,
            condition.1,
            condition.2,
            pvalue,
            padj,
            direction
        ) %>%
        dplyr::arrange(order(mixedsort(chromosome)), start, end, condition.1, condition.2)

    return (object)
}

##- detectCompartments -------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Detect compartments for each genomic position in each condition, and compute
#' p-values for compartment differences between conditions.
#'
#' @param object  A \code{HiCDOCExp} object.
#'
#' @return A \code{HiCDOCExp} object, with centroids, compartments, distances,
#' concordances and differences.
#'
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' @export
detectCompartments <- function(object) {
    message("Clustering...")
    object <- clusterize(object)
    message("Predicting compartments...")
    object <- predictAB(object)
    message("Computing p-values...")
    object <- computePValues(object)
    message("Done.")

    return (object)
}
