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

##- tieCentroids -------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Assign correct cluster labels by comparing centroids across conditions.
#'
#' @param object A \code{HiCDOCDataSet} object.
#'
#' @return A \code{HiCDOCDataSet} object, with diagonalRatios, and with corrected
#' cluster labels in centroids, compartments, distances and concordances.
tieCentroids <- function(object) {
    totalConditions <- length(unique(object@conditions))
    referenceCondition <- object@conditions[[1]]

    clusters <- object@centroids %>%
        dplyr::filter(condition == referenceCondition) %>%
        tidyr::pivot_wider(names_from = compartment,
                           values_from = centroid,
                           names_prefix = "reference.") %>%
        dplyr::select(-c(condition)) %>%
        tidyr::uncount(totalConditions) %>%
        dplyr::bind_cols(
            object@centroids %>%
                tidyr::pivot_wider(names_from = compartment,
                                   values_from = centroid,
                                   names_prefix = "centroid.") %>%
                dplyr::select(-c(chromosome))
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(cl1 = dplyr::if_else(
            euclideanDistance(centroid.1, reference.1) *
                euclideanDistance(centroid.2, reference.2) <
                euclideanDistance(centroid.1, reference.2) *
                euclideanDistance(centroid.2, reference.1), 1, 2
        )) %>%
        dplyr::mutate(cl2 = dplyr::if_else(cl1 == 1, 2, 1)) %>%
        dplyr::select(-c(centroid.1, centroid.2, reference.1, reference.2))

    object@compartments %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::mutate(compartment = dplyr::if_else(compartment == 1, cl1, cl2)) %>%
        dplyr::select(-c(cl1, cl2))

    object@concordances %<>%
        dplyr::left_join(clusters, by = c("chromosome", "condition")) %>%
        dplyr::mutate(change = dplyr::if_else(
            compartment == 1,
            dplyr::if_else(compartment == cl1, 1, -1),
            dplyr::if_else(compartment == cl2, 1, -1)
        )) %>%
        dplyr::mutate(concordance = change * concordance) %>%
        dplyr::mutate(compartment = dplyr::if_else(compartment == 1, cl1, cl2)) %>%
        dplyr::select(-c(cl1, cl2, change))

    object@distances %<>%
        dplyr::left_join(clusters, by=c("chromosome", "condition")) %>%
        dplyr::mutate(compartment = dplyr::if_else(compartment == 1, cl1, cl2)) %>%
        dplyr::select(-c(cl1, cl2))

    object@centroids %<>%
        dplyr::left_join(clusters, by=c("chromosome", "condition")) %>%
        dplyr::mutate(compartment = dplyr::if_else(compartment == "1", cl1, cl2)) %>%
        dplyr::select(-c(cl1, cl2))

    return (object)
}

##- clusterize for 1 chromosome and 1 condition-------------------------------#
##----------------------------------------------------------------------------#
#' Segregate genomic positions into 2 clusters using constrained k-means.
#'
#' @param object  A \code{HiCDOCDataSet} object.
#' @param chromosomeId A character or numeric value.
#' Name or number of the chromosome
#' @param conditionId A character or numeric value.
#' Name or number of the condition
#'
#' @return A \code{HiCDOCDataSet} object, with:
#' - centroid (vector) for each cluster
#' - compartment (cluster number) for each genomic position
#' - distances to centroids (float) for each genomic position in each replicate
#' - concordance (float in [-1, 1]) for each genomic position in each replicate
clusterizeChrCond <- function(object, chromosomeId, conditionId) {
    chr <- testChromosome(object, chromosomeId)
    cond <- testCondition(object, conditionId)

    totalBinsChr <- object@totalBins[[chr]]
    if (totalBinsChr == -Inf)
        return(NULL)

    positions <- (seq.int(totalBinsChr) - 1) * object@binSize

    # Correct for filtered bins
    if (!is.null(object@weakBins[[chromosomeId]])) {
        positions <- positions[-object@weakBins[[chromosomeId]]]
        totalBinsChr <- totalBinsChr - length(object@weakBins[[chromosomeId]])
    }

    replicates <- object@replicates[which(object@conditions == cond)]

    replicateInteractions <- purrr::map(replicates, function(x)
        sparseInteractionsToMatrix(object,
                                   chr,
                                   cond,
                                   x,filter = TRUE))
    interactions <- do.call("rbind", replicateInteractions)

    mustLink <- matrix(
        rep(0:(length(replicates) - 1), totalBinsChr) * totalBinsChr +
            rep(0:(totalBinsChr - 1), each = length(replicates)),
        nrow = totalBinsChr,
        byrow = TRUE
    )

    clusteringOutput <- constrainedClustering(
        interactions,
        mustLink,
        object@parameters$kMeansDelta,
        object@parameters$kMeansIterations,
        object@parameters$kMeansRestarts
    )

    clusters <- clusteringOutput[["clusters"]][0:totalBinsChr] + 1
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

    dfCompartments <- tibble(
        chromosome = factor(chr, levels = object@chromosomes),
        position = positions,
        condition = factor(cond, levels = unique(object@conditions)),
        compartment = factor(clusters, levels = c(1,2))
    )

    dfConcordances <-
        purrr::map_dfr(seq_len(length(replicates)), ~dfCompartments) %>%
        dplyr::mutate(replicate = rep(factor(replicates, levels = unique(object@replicates)),
                                      each = ncol(interactions)),
                      concordance = concordances) %>%
        dplyr::select(chromosome, position, condition, replicate, compartment, concordance)

    dfDistances <-
        purrr::map_dfr(seq_len(2), ~dfConcordances) %>%
        dplyr::select(-concordance) %>%
        dplyr::mutate(
            compartment = rep(factor(c(1,2)),
                              each = length(replicates) * ncol(interactions)),
            distance = c(t(distances))
        )

    dfCentroids <- tibble(
        chromosome = factor(chr, levels = object@chromosomes),
        condition = factor(cond, levels = unique(object@conditions)),
        compartment = factor(c(1,2)),
        centroid = centroids
    )

    return(list(
        "compartments" = dfCompartments,
        "concordances" = dfConcordances,
        "distances" = dfDistances,
        "centroids" = dfCentroids
    ))
}


##- clusterize ---------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Segregate genomic positions into 2 clusters using constrained k-means.
#'
#' @param object  A \code{HiCDOCDataSet} object.
#'
#' @return A \code{HiCDOCDataSet} object, with:
#' - centroid (vector) for each cluster
#' - compartment (cluster number) for each genomic position
#' - distances to centroids (float) for each genomic position in each replicate
#' - concordance (float in [-1, 1]) for each genomic position in each replicate
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- clusterize(object)
clusterize <- function(object) {
    object@parameters <- checkParameters(object@parameters,
                                         c("kMeansDelta",
                                           "kMeansIterations",
                                           "kMeansRestarts"))
    vectChr <- rep(object@chromosomes, each = length(unique(object@conditions)))
    vectCond <- rep(unique(object@conditions), length(object@chromosomes))

    clusterRes <- purrr::map2(vectChr, vectCond,
                              function(.x, .y) clusterizeChrCond(object, .x, .y))
    object@compartments <- purrr::map_dfr(clusterRes, "compartments")
    object@concordances <- purrr::map_dfr(clusterRes, "concordances")
    object@distances <- purrr::map_dfr(clusterRes, "distances")
    object@centroids <- purrr::map_dfr(clusterRes, "centroids")

    object <- tieCentroids(object)
    return(object)
}

##- diagonalRatios -----------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Compute the ratio of the reads which are on the diagonal vs those off-
#'   diagonal, for each bin, and each matrix.
#'
#' @param object an HiCDOCDataSet object
#' @param chromosomeId A character or numeric value.
#' Name or number of the chromosome
#' @param conditionId A character or numeric value.
#' Name or number of the condition
#' @param replicateId A character or numeric value.
#' Name or number of the replicate
#'
#' @return a \code{tibble}.
diagonalRatios <- function(object, chromosomeId, conditionId, replicateId){
    interactions = object@interactions %>%
        dplyr::filter(chromosome == chromosomeId) %>%
        dplyr::filter(condition == conditionId) %>%
        dplyr::filter(replicate == replicateId)
    
    diagonal <- interactions %>%
        filter(position.1 == position.2) %>%
        select(-position.2) %>%
        rename(position = position.1) %>%
        rename(diagonal = value)
    
    offDiagonal <- interactions %>%
        filter(position.1 != position.2) %>%
        pivot_longer(cols=starts_with("position"), names_to="namepos",
                     values_to="position") %>%
                     select(-namepos) %>%
                     group_by(chromosome, condition, replicate, position) %>%
                     summarise(offDiagonal = median(value)) %>%
                     ungroup()
                     
    diagonalRatios <- diagonal %>%
        left_join(offDiagonal,
                  by = c("chromosome",
                         "condition",
                         "replicate",
                         "position")) %>%
        mutate(value = diagonal - offDiagonal) %>%
        select(-c(diagonal, offDiagonal))
    
    return (diagonalRatios)
}

##- predictAB ----------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Use ratio between diagonal and off-diagonal interactions to determine which
#' clusters correspond to compartments A and B.
#'
#' @param object A \code{HiCDOCDataSet} object.
#'
#' @return A \code{HiCDOCDataSet} object, with diagonalRatios, and with A and B
#' labels replacing cluster numbers in centroids, compartments, distances and
#' concordances.
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- clusterize(object)
#' object <- predictAB(object)
predictAB <- function(object) {
    chromosomeIds = rep(object@chromosomes, each = length(object@replicates))
    conditionIds = rep(object@conditions, length(object@chromosomes))
    replicateIds = rep(object@replicates, length(object@chromosomes))

    object@diagonalRatios <- purrr::pmap_dfr(
        list(chromosomeIds, conditionIds, replicateIds),
        function(x, y, z)
            diagonalRatios(object, x, y, z))

    compartments <- object@compartments %>%
        dplyr::left_join(
            object@diagonalRatios,
            by = c("chromosome", "condition", "position")
        ) %>%
        dplyr::group_by(chromosome, compartment) %>%
        dplyr::summarize(value = median(value)) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = compartment,
                           values_from = value,
                           values_fill = 0,
                           names_prefix = "val") %>%
        dplyr::mutate(A = dplyr::if_else(val1 >= val2, 2, 1)) %>%
        dplyr::select(-c(val1, val2))

    object@compartments %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(compartment = factor(dplyr::if_else(compartment == A, "A", "B"))) %>%
        dplyr::select(-c(A))

    object@concordances %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(change = dplyr::if_else(A == 1, 1, -1)) %>%
        dplyr::mutate(concordance = change * concordance) %>%
        dplyr::select(-c(A, change))

    object@distances %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(compartment = factor(dplyr::if_else(compartment == A, "A", "B"))) %>%
        dplyr::select(-c(A))

    object@centroids %<>%
        dplyr::left_join(compartments, by = c("chromosome")) %>%
        dplyr::mutate(compartment = factor(dplyr::if_else(compartment == A, "A", "B"))) %>%
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
#' @param object A \code{HiCDOCDataSet} object.
#'
#' @return A \code{HiCDOCDataSet} object, with differences and their p-values.
computePValues <- function(object) {
    # Compute median of differences between pairs of concordances
    # N.b. median of differences != difference of medians
    totalReplicates = length(object@replicates)

    concordanceDifferences <- object@concordances %>%
        dplyr::arrange(chromosome, position, condition, replicate)

    concordanceDifferences %<>%
        # Duplicate each row "totalReplicates" times
        tidyr::uncount(totalReplicates) %>%
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
        tidyr::uncount(totalConditions) %>%
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
        dplyr::mutate(H0_value = dplyr::if_else(compartment.1 == compartment.2, value, NA_real_)) %>%
        dplyr::group_by(condition.1, condition.2) %>%
        dplyr::mutate(quantile = ecdf(H0_value)(value)) %>%
        dplyr::filter(compartment.1 != compartment.2) %>%
        dplyr::mutate(pvalue = 1 - quantile) %>%
        dplyr::mutate(pvalue = dplyr::if_else(pvalue < 0, 0, pvalue)) %>%
        dplyr::mutate(pvalue = dplyr::if_else(pvalue > 1, 1, pvalue)) %>%
        dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) %>%
        dplyr::ungroup() %>%
        dplyr::rename(start = position) %>%
        dplyr::mutate(end = start + object@binSize) %>%
        dplyr::mutate(direction = factor(dplyr::if_else(compartment.1 == "A", "A->B", "B->A"))) %>%
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
#' @param object  A \code{HiCDOCDataSet} object.
#'
#' @return A \code{HiCDOCDataSet} object, with centroids, compartments, distances,
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
    testSlotsHiCDOC(
        object,
        slots = c(
            "chromosomes",
            "conditions",
            "totalBins",
            "binSize",
            "weakBins",
            "interactions"
        )
    )
    message("Clustering...")
    object <- clusterize(object)
    message("Predicting compartments...")
    object <- predictAB(object)
    message("Computing p-values...")
    object <- computePValues(object)
    message("Done.")

    return (object)
}
