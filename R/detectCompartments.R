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

    validChromosomeNames <-
        names(base::Filter(function(x) !is.null(x), object@validAssay))
    referenceConditionNames <-
        vapply(validChromosomeNames,
               FUN = function(x) 
                   sort(object$condition[object@validAssay[[x]]])[1],
               FUN.VALUE = "")
    
    refCentroids <- data.table::merge.data.table(
        object@centroids[compartment == 1 & 
                             condition == referenceConditionNames[chromosome], 
                         .(chromosome, ref.1 = centroid)],
        object@centroids[compartment == 2 & 
                             condition == referenceConditionNames[chromosome], 
                         .(chromosome, ref.2 = centroid)],
        all=T
    )
    
    clusters <- data.table::merge.data.table(
        object@centroids[compartment == 1, 
                         .(chromosome, condition, centroid.1 = centroid)],
        object@centroids[compartment == 2, 
                         .(chromosome, condition, centroid.2 = centroid)],
        all=T
    )
    clusters <- data.table::merge.data.table(
        clusters,
        refCentroids, 
        all=TRUE)
   
    c1_r1 <- mapply(function(x, y) .euclideanDistance(unlist(x), unlist(y)), 
                    clusters$centroid.1, clusters$ref.1)
    c1_r2 <- mapply(function(x, y) .euclideanDistance(unlist(x), unlist(y)), 
                    clusters$centroid.1, clusters$ref.2)
    c2_r1 <- mapply(function(x, y) .euclideanDistance(unlist(x), unlist(y)), 
                    clusters$centroid.2, clusters$ref.1)
    c2_r2 <- mapply(function(x, y) .euclideanDistance(unlist(x), unlist(y)), 
                    clusters$centroid.2, clusters$ref.2)
    clusters[,cluster.1 := 1*((c1_r1*c2_r2) >= (c1_r2 * c2_r1)) + 1] 
    clusters[,cluster.2 := 1 + (cluster.1 == 1)]
    clusters <- clusters[,.(chromosome, condition, cluster.1, cluster.2)]

    object@compartments  <- data.table::merge.data.table(
        object@compartments,
        clusters,
        by=c("chromosome", "condition"),
        all.x=T,
        sort = FALSE
    )
    object@compartments[,compartment := ifelse(compartment == 1, cluster.1, cluster.2)]
    object@compartments[ ,`:=`(cluster.1 = NULL, cluster.2 = NULL)]
    
    object@concordances <- data.table::merge.data.table(
        object@concordances,
        clusters,
        by = c("chromosome", "condition"),
        all.x=TRUE,
        sort=FALSE
    )
    
    object@concordances[, change := -1]
    object@concordances[compartment == 1 & compartment == cluster.1, change := 1]
    object@concordances[compartment == 2 & compartment == cluster.2, change := 1]
    object@concordances[, concordance := change * concordance]
    object@concordances[, compartment := data.table::fifelse(compartment == 1, cluster.1, cluster.2)]
    object@concordances[,`:=`(cluster.1 = NULL, cluster.2 = NULL, change = NULL)]

    object@distances <- data.table::merge.data.table(
        object@distances,
        clusters,
        by = c("chromosome", "condition"),
        all.x=TRUE,
        sort=FALSE
    )
    object@distances[,compartment := 
                         data.table::fifelse(
                             compartment == 1, cluster.1, cluster.2)]
    object@distances[,`:=`(cluster.1 = NULL, cluster.2 = NULL)]
    
    object@centroids <- data.table::merge.data.table(
        object@centroids,
        clusters,
        by = c("chromosome", "condition"),
        all.x=TRUE,
        sort=FALSE
    )
    object@centroids[,compartment := 
                         data.table::fifelse(
                             compartment == 1, cluster.1, cluster.2)]
    object@centroids[,`:=`(cluster.1 = NULL, cluster.2 = NULL)]
    
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
    object
) {
    chromosomeName <- object@chromosomes
    conditionName <- object$condition[1]
    
    nbBins <- length(InteractionSet::regions(object))
    validAssay <- object@validAssay[[chromosomeName]]
    replicateNames <- object$replicate[validAssay]
    orderAssay <- validAssay[order(replicateNames)]
    
    if (length(replicateNames) == 0) return(NULL)
    isetChromosome <- InteractionSet::InteractionSet(
        SummarizedExperiment::assay(object),
        InteractionSet::interactions(object)
    )
    matAssay <- 
        lapply(orderAssay,
               FUN = function(x) {
                   InteractionSet::inflate(isetChromosome, 
                                           rows = chromosomeName,
                                           columns = chromosomeName,
                                           sample = x)
               })
    matAssay <- lapply(matAssay, function(x) x@matrix)
    matAssay <- do.call("rbind", matAssay)
    matAssay[is.na(matAssay)] <- 0
    
    mustLink <- .constructLinkMatrix(length(replicateNames), nbBins)
    clusteringOutput <-
        constrainedClustering(
            matAssay,
            mustLink,
            object@parameters$kMeansDelta,
            object@parameters$kMeansIterations,
            object@parameters$kMeansRestarts
        )
    # TODO : question : pourquoi on ne prend que les premiers ?
    # Quel est l'intérêt de retourner 2 fois ?
    clusters <- clusteringOutput[["clusters"]][0:nbBins] + 1
    centroids <- clusteringOutput[["centroids"]]

    min <- .distanceRatio(centroids[[1]], centroids)
    max <- .distanceRatio(centroids[[2]], centroids)

    concordances <-
        apply(matAssay, 1, function(row) {
            2 * (.distanceRatio(row, centroids) - min) / (max - min) - 1
        })

    distances <-
        apply(matAssay, 1, function(row) {
            c(
                .euclideanDistance(row, centroids[[1]]),
                .euclideanDistance(row, centroids[[2]])
            )
        })
    
    indexes <- S4Vectors::mcols(InteractionSet::regions(object))$index
    dfCompartments <- data.table::data.table(
        "chromosome" = chromosomeName,
        "index" = indexes,
        "condition" = conditionName,
        "compartment" = clusters
    )
    
    dfConcordances <- data.table::data.table(
        "chromosome" = chromosomeName,
        "index" = rep(indexes, length(replicateNames)),
        "condition" = conditionName,
        "replicate" = rep(sort(replicateNames), each = nbBins),
        "compartment" = rep(clusters, length(replicateNames)),
        "concordance" = concordances
    )
        
    dfDistances <- data.table::data.table(
        "chromosome" = chromosomeName,
        "index" = rep(indexes, 2 * length(replicateNames)),
        "condition" = conditionName,
        "replicate" = rep(rep(sort(replicateNames), each = nbBins), 2),
        "compartment" = rep(c(1, 2), each = length(replicateNames) * nbBins),
        "distance" = c(t(distances)))
    
    dfCentroids <- data.table::data.table(
        "chromosome" = chromosomeName,
        "condition" = conditionName,
        "compartment" = c(1, 2),
        "centroid" = centroids)
    
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
    
    condByChromosomes <- 
        lapply(object@chromosomes,
               FUN = function(x) {
                   data.frame("chr" = x, 
                              "cond" = sort(unique(
                                  object$condition[object@validAssay[[x]]])))
               })
    condByChromosomes <- Reduce(rbind, condByChromosomes)
    
    reducedObjects <- 
        mapply(function(x, y){
            reduceHiCDOCDataSet(object,
                                chromosomes = x,
                                conditions = y,
                                dropLevels = TRUE)
        }, condByChromosomes$chr, 
           condByChromosomes$cond)
    
    result <- .internalApply(
        parallel,
        reducedObjects,
        .clusterizeChromosomeCondition,
        type="lapply"
        )
    
    compartments <- lapply(result, function(x) x[["compartments"]])
    concordances <- lapply(result, function(x) x[["concordances"]])
    distances <- lapply(result, function(x) x[["distances"]])
    centroids <- lapply(result, function(x) x[["centroids"]])
    
    object@compartments <- data.table::rbindlist(compartments)
    object@concordances <- data.table::rbindlist(concordances)
    object@distances <- data.table::rbindlist(distances)
    object@centroids <- data.table::rbindlist(centroids)
    
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
    object
) {
    # TODO : j'ai essayé d'utiliser anchorIds() mais ça ne me renvoie pas les
    # bonnes valeurs d'index.
    tmp <- as.data.frame(InteractionSet::interactions(object))
    ids <- list("first" = tmp$index1, "second" = tmp$index2)
    diagonal <- ids$first == ids$second
    cn <- paste(object$condition, object$replicate)
    
    # Values on diagonal
    onDiagonal <- data.table(
        "chromosome" = tmp[diagonal,]$seqnames1,
        "index" = ids$first[diagonal],
        SummarizedExperiment::assay(object)[diagonal,])
    data.table::setnames(onDiagonal, c("chromosome", "index", cn))
    onDiagonal <- data.table::melt.data.table(onDiagonal, 
                                              id.vars=c("chromosome", "index"))
    onDiagonal <- onDiagonal[!is.na(value)]
    
    # Compute median by bin, out of diagonal
    offDiagonal <- rbind(
        SummarizedExperiment::assay(object)[!diagonal,],
        SummarizedExperiment::assay(object)[!diagonal,]
    ) 
    offDiagonal <- data.table::data.table(
        "index" = c(ids$first[!diagonal], ids$second[!diagonal]),
              offDiagonal)
    setnames(offDiagonal, c("index", cn))
    offDiagonal <- offDiagonal[, lapply(.SD, 
                                        median, 
                                        na.rm=TRUE), 
                               by=index, 
                               .SDcols=colnames(offDiagonal)[-1]] 
    offDiagonal <- data.table::melt.data.table(offDiagonal, 
                                               id.vars = "index", 
                                               value.name = "median")
    offDiagonal <- offDiagonal[!is.na(median)]
    
    onDiagonal<- data.table::merge.data.table(onDiagonal, 
                                              offDiagonal, 
                                              all=T,
                                              by=c("index", "variable"),
                                              sort=FALSE)
    onDiagonal[is.na(value), value := 0]
    onDiagonal[is.na(median), median := 0]
    onDiagonal[,ratio := value - median]
    onDiagonal[, c("condition", "replicate") := 
                   tstrsplit(variable, " ", fixed=TRUE)]
    onDiagonal <- onDiagonal[,.(chromosome, index, condition, replicate, ratio)]
    return(onDiagonal)
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
    
    ratios <-.computeSelfInteractionRatios(object)
    object@selfInteractionRatios <-ratios
    
    # TODO : est-ce utile de garder le slot selfInteractionRatios ?
    compartments <- data.table::merge.data.table(
        object@compartments, 
        object@selfInteractionRatios,
        by = c("chromosome", "index", "condition"),
        all.x=T,
        sort=FALSE
    )
    # TODO !! A vérifier ici, les valeurs ne sont pas bonnes !!!
    # La sortie de selfInteractionRatios est bonne 
    compartments <- compartments[,.(ratio = median(ratio, na.rm=T)), 
                                 by=.(chromosome, compartment)]
    compartments <- data.table::dcast(compartments, 
                                       chromosome  ~ compartment, 
                                       value.var = "ratio",
                                      fill=0)
    
    compartments[, A := data.table::fifelse(`1` >= `2`, 2, 1)]
    compartments <- compartments[, .(chromosome, A)]
    
    object@compartments <- data.table::merge.data.table(
        object@compartments, 
        compartments, 
        by="chromosome",
        all.x=T
    )
    
    object@compartments[,compartment := 
                            data.table::fifelse(compartment == A, "A", "B")]
    object@compartments[,compartment := factor(compartment, levels=c("A", "B"))]
    object@compartments[,A := NULL]
   
    object@concordances <- data.table::merge.data.table(
        object@concordances, 
        compartments, 
        by="chromosome",
        all.x=T
    )
    object@concordances[,change := 
                            data.table::fifelse(A == 1, 1, -1)]
    object@concordances[,concordance :=  change * concordance]
    object@concordances[,compartment := factor(data.table::fifelse(
        compartment == A, "A", "B"), levels=c("A", "B"))]
    object@concordances[,change := NULL]
    object@concordances[,A := NULL]
    
    object@distances <- data.table::merge.data.table(
        object@distances, 
        compartments, 
        by="chromosome",
        all.x=T
    )
    object@distances[,compartment := factor(data.table::fifelse(
        compartment == A, "A", "B"), levels=c("A", "B"))]
    object@distances[,A := NULL]
    
    object@centroids <- data.table::merge.data.table(
        object@centroids, 
        compartments, 
        by="chromosome",
        all.x=T
    )
    object@centroids[,compartment := factor(data.table::fifelse(
        compartment == A, "A", "B"), levels=c("A", "B"))]
    object@centroids[,A := NULL]
    
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
    nbReplicates <- length(object$replicate)
    concordances <- object@concordances 
    concordances[,condition := factor(condition, 
                                     levels=sort(unique(object$condition)))]
    data.table::setorder(concordances, chromosome, index, condition, replicate)
    diff1 <- concordances[rep(seq_len(nrow(concordances)), each=nbReplicates),
                         .(chromosome, 
                           index,
                           condition.1 = condition, 
                           concordance.1 = concordance)]
    diff2 <-  concordances[rep(seq_len(nrow(concordances)), nbReplicates),
                          .(chromosome = chromosome, 
                            index = index,
                            condition.2 = condition, 
                            concordance.2 = concordance)]
    data.table::setorder(diff2, chromosome, index)
    diff2[,`:=`(chromosome = NULL, index = NULL)]
    diffConcordance <- cbind(diff1, diff2)
    rm(diff1, diff2)
    diffConcordance <- diffConcordance[as.numeric(condition.1) < as.numeric(condition.2)]
    diffConcordance <- diffConcordance[,.(difference = stats::median(abs(concordance.1 - concordance.2))),
                                  by = .(chromosome, index, condition.1, condition.2)]
 
    # Format compartments per pair of conditions
    # Join medians of differences and pairs of conditions
    nbConditions <- length(unique(object$condition))
    compartments <- object@compartments 
    compartments[,condition := factor(condition, 
                                     levels=sort(unique(object$condition)))]
    data.table::setorder(compartments, chromosome, index, condition)
    comp1 <- compartments[rep(seq_len(nrow(compartments)), each=nbConditions),
                          .(chromosome, 
                            index,
                            condition.1 = condition, 
                            compartment.1 = compartment)]
    comp2 <- compartments[rep(seq_len(nrow(compartments)), nbConditions),
                          .(chromosome, 
                            index,
                            condition.2 = condition, 
                            compartment.2 = compartment)]
    data.table::setorder(comp2, chromosome, index)
    comp2[,`:=`(chromosome = NULL, index = NULL)]
    comparisons <- cbind(comp1, comp2)
    rm(comp1, comp2)
    comparisons <- comparisons[as.numeric(condition.1) < as.numeric(condition.2)]
    comparisons <- data.table::merge.data.table(comparisons, 
                                                diffConcordance,
                                                by = c("chromosome", "index", "condition.1", "condition.2"))
    data.table::setcolorder(comparisons, c("chromosome",
                               "index",
                               "condition.1",
                               "condition.2",
                               "compartment.1",
                               "compartment.2",
                               "difference"))
    object@comparisons <- comparisons
    
    # Compute p-values for switching positions
    # P-values for a condition pair computed from the whole genome distribution
    differences <- copy(comparisons)
    differences[compartment.1 == compartment.2 ,H0_value := difference]
    data.table::setorder(differences, condition.1, condition.2)
    quantiles <- split(differences, list(differences$condition.1, differences$condition.2))
    quantiles <- lapply(quantiles, function(x) x[difference > 0])
    quantiles <- lapply(quantiles, function(x) 
        if(nrow(x)>0) { return(stats::ecdf(x$H0_value)(x$difference))} else return(NULL))
    quantiles <- do.call("c", quantiles)
    differences[difference > 0, quantile := quantiles]
    
    # Pvalues
    differences <- differences[compartment.1 != compartment.2]
    differences[, pvalue := 1 - quantile]
    differences[pvalue<0, pvalue := 0]
    differences[pvalue>1, pvalue := 1]
    pvalAdjust <- split(differences, list(differences$condition.1, differences$condition.2))
    pvalAdjust <- lapply(pvalAdjust, function(x) if(nrow(x) > 0) return(stats::p.adjust(x$pvalue, method = "BH")) else return(NULL))
    pvalAdjust <- do.call("c", pvalAdjust)
    differences[, pvalue.adjusted := pvalAdjust]
    
    # Changes
    differences[,direction := data.table::fifelse(compartment.1 == "A", "A->B", "B->A")]
    differences[,direction := factor(direction, levels = c("A->B", "B->A"))]
    differences <- differences[,.(chromosome, 
                                  index, 
                                  condition.1, 
                                  condition.2, 
                                  pvalue, 
                                  pvalue.adjusted, 
                                  direction)]
    data.table::setorder(differences, chromosome, index, condition.1, condition.2)
    object@differences <- differences
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
#' \code{\link[BiocParallel]{bpmapply}} function. Before to call the
#' function in parallel you should specify the parallel parameters such as:
#'     \itemize{
#'         \item{On Linux:
#'
#'              \code{multiParam <- BiocParallel::MulticoreParam(workers = 10)}
#'          }
#'          \item{On Windows:
#'
#'              \code{multiParam <- BiocParallel::SnowParam(workers = 10)}
#'         }
#'     }
#'     And then you can register the parameters to be used by BiocParallel:
#'
#'     \code{BiocParallel::register(multiParam, default = TRUE)}
#'
#'     You should be aware that using MulticoreParam, reproducibility of the
#'     detectCompartments function using a RNGseed may not work. See this
#'     \href{https://github.com/Bioconductor/BiocParallel/issues/122}{issue}
#'     for more details.
#' }
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @param parallel
#' Whether or not to parallelize the processing. Defaults to FALSE
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
#' data(exampleHiCDOCDataSet)

#' ## Run all filtering and normalization steps (not run for timing reasons)
#' # object <- filterSmallChromosomes(exampleHiCDOCDataSet)
#' # object <- filterSparseReplicates(object)
#' # object <- filterWeakPositions(object)
#' # object <- normalizeTechnicalBiases(object)
#' # object <- normalizeBiologicalBiases(object)
#' # object <- normalizeDistanceEffect(object)
#'
#' # Detect compartments and differences across conditions
#' object <- detectCompartments(exampleHiCDOCDataSet)
#'
#' @usage
#' detectCompartments(
#'     object,
#'     parallel = FALSE,
#'     kMeansDelta = NULL,
#'     kMeansIterations = NULL,
#'     kMeansRestarts = NULL
#' )
#'
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
            "validAssay",
            "parameters"
        )
    )

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
    object <- .tieCentroids(object)
    
    message("Predicting A/B compartments.")
    object <- .predictCompartmentsAB(object, parallel)
    message("Detecting significant differences.")
    object <- .computePValues(object)
    
    # Reformating outputs
    chr <- object@chromosomes
    cond <- sort(unique(object$condition))
    rep <- sort(unique(object$replicate))
    
    object@centroids[,`:=`(chromosome = factor(chromosome, levels =  chr),
                           condition = factor(condition, levels = cond))]
    object@differences[, chromosome := factor(chromosome, levels =  chr)]
    object@compartments[, chromosome := factor(chromosome, levels =  chr)]
    object@concordances[,`:=`(chromosome = factor(chromosome, levels =  chr),
                              replicate = factor(replicate, levels = rep))]
    object@distances[,`:=`(chromosome = factor(chromosome, levels =  chr),
                           condition = factor(condition, levels = cond),
                           replicate = factor(replicate, levels = rep))]
    object@comparisons[,chromosome := factor(chromosome, levels =  chr)]
    object@selfInteractionRatios[,`:=`(chromosome = factor(chromosome, levels =  chr),
                           condition = factor(condition, levels = cond),
                           replicate = factor(replicate, levels = rep))]
    
    # TODO : objets à trier. 
    return(object)
}
