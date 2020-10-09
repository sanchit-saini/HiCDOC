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
    return (log((
        euclideanDistance(x, centroids[[1]]) + eps
    )
    / (
        euclideanDistance(x, centroids[[2]]) + eps
    )))
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
        if (totalBins == -Inf)
            next
        
        positions <- (seq.int(totalBins) - 1) * object@binSize
        
        # Correct for filtered bins
        if (!is.null(object@weakBins[[chromosomeId]])) {
            positions <- positions[-object@weakBins[[chromosomeId]]]
            if (totalBins < length(object@weakBins[[chromosomeId]])) {
                message("Problem while filtering bins")
                message("Chr: ", chromosomeId)
                message("# bins: ", totalBins)
                message(
                    "max bins: ",
                    max(
                        object@interactions$position.1,
                        object@interactions$position.2
                    )
                )
                message("bin size: ", object@binSize)
                message("# weak bins: ", length(object@weakBins[[chromosomeId]]))
                message(
                    "weak bins: ",
                    paste(object@weakBins[[chromosomeId]]),
                    sep = " ",
                    collapse = " "
                )
            }
            totalBins <-
                totalBins - length(object@weakBins[[chromosomeId]])
        }
        
        for (conditionId in unique(object@conditions)) {
            interactions <- matrix(nrow = 0, ncol = totalBins)
            replicates <-
                object@replicates[which(object@conditions == conditionId)]
            
            for (replicateId in replicates) {
                replicateInteractions <- sparseInteractionsToMatrix(object,
                                                                    chromosomeId,
                                                                    conditionId,
                                                                    replicateId,
                                                                    filter = TRUE)
                interactions <- rbind(interactions, replicateInteractions)
            }
            
            mustLink <- matrix(
                rep(0:(length(
                    replicates
                ) - 1), totalBins) * totalBins
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
            
            min <- distanceRatio(centroids[[1]],
                                 centroids)
            
            max <- distanceRatio(centroids[[2]],
                                 centroids)
            
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
                    value = clusters
                )
            )
            
            object@concordances %<>% bind_rows(
                tibble(
                    chromosome = chromosomeId,
                    position = rep(positions, length(replicates)),
                    condition = conditionId,
                    replicate = rep(replicates, each = ncol(interactions)),
                    value = concordances
                )
            )
            
            object@distances %<>% bind_rows(
                tibble(
                    chromosome = chromosomeId,
                    position = rep(rep(
                        positions, length(replicates)
                    ), 2),
                    condition = conditionId,
                    replicate = rep(rep(
                        replicates, each = ncol(interactions)
                    ), 2),
                    cluster = c(rep(
                        1, length(replicates) * ncol(interactions)
                    ),
                    rep(
                        2, length(replicates) * ncol(interactions)
                    )),
                    value = c(t(distances))
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
    
    object@diagonalRatios <-
        bind_rows(apply(groups, 1, function(group) {
            fullInteractions <- sparseInteractionsToFullInteractions(object,
                                                                     group[[1]],
                                                                     group[[2]],
                                                                     group[[3]])
            
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
                left_join(offDiagonal,
                          by = c("chromosome", "condition", "replicate", "position")) %>%
                mutate(value = diagonal - offDiagonal) %>%
                select(-c(diagonal, offDiagonal))
            
            return (diagonalRatios)
        }))
    
    compartments <- object@compartments %>%
        rename(compartment = value) %>%
        left_join(object@diagonalRatios,
                  by = c("chromosome", "condition", "position")) %>%
        group_by(chromosome, condition, compartment) %>%
        summarize(value = median(value)) %>%
        ungroup() %>%
        spread(key = compartment,
               value = value,
               fill = 0) %>%
        mutate(A = if_else(`1` >= `2`, 2, 1)) %>%
        select(-c(`1`, `2`))
    
    object@compartments %<>%
        left_join(compartments, by = c("chromosome", "condition")) %>%
        mutate(value = factor(if_else(value == A, "A", "B"))) %>%
        select(-c(A))
    
    object@concordances %<>%
        left_join(compartments, by = c("chromosome", "condition")) %>%
        mutate(change = if_else(A == 1, 1,-1)) %>%
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
            concordanceDifferences[rep(seq_len(nrow(concordanceDifferences)),
                                       totalReplicates), ] %>%
                # Arrange by chromosome and position to join with the duplicated rows
                arrange(chromosome, position) %>%
                rename_with(function(old_colnames)
                    paste(old_colnames, 2, sep = "."))
        ) %>%
        select(-chromosome.2,-position.2) %>%
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
                                       totalConditions), ] %>%
                # Arrange by chromosome and position to join with the duplicated rows
                arrange(chromosome, position) %>%
                rename_with(function(old_colnames)
                    paste(old_colnames, 2, sep = "."))
        ) %>%
        select(-chromosome.2,-position.2) %>%
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
        mutate(H0_value = if_else(compartment.1 == compartment.2, value, NA_real_)) %>%
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
