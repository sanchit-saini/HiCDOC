#' Plot the interaction matrix (as heatmap).
#'
#' @param objet an \code{HiCDOCExp} object
#' @param chromosomeId The name or number of the chromosome to plot.
#' If number, will be taken in \code{object@chromosomes[chromosomeId]}
#' @param trans character: transformation of the color scale. Default to "log2". 
#' See \code{\link[ggplot2::scale_fill_gradient]{scale_fill_gradient}} for other accepted values. 
#' Set to NULL for no transformation.
#' @return A \code{ggplot} object.
#' @examples
#' object <- HiCDOCExample()
#' p <- plotInteractionMatrix(object, chromosomeId = 1, trans = "log2")
#' @export
plotInteractionMatrix <- function(object, chromosomeId, trans = "log2") {
  testSlotsHiCDOCExp(object,
                     slots = c("interactions"))
  chr <- testchromosome(object, chromosomeId)
  
  interactionsChr <- object@interactions %>% 
    filter(chromosome == chr & value>0) 
  
  # Duplicate positions (sparse symetric matrix)
  interactionsChr <- interactionsChr %>%
    filter(position.1 != position.2) %>%
    rename(position.2 = position.1, position.1 = position.2) %>%
    bind_rows(interactionsChr) %>%
    group_by(condition, replicate, position.1, position.2) %>%
    mutate(intensity = mean(value)) %>%
    unite("rep_cond", replicate, condition, sep="_")
  
  xylim <- c(0, (object@totalBins[[chr]] - 1) * object@binSize)
 
  if (nrow(interactionsChr) > 0) {
    p <- interactionsChr %>%
      ggplot(aes(x = position.1, y = position.2, z = intensity)) +
        geom_tile(aes(fill = intensity)) +
        coord_fixed(ratio = 1) +
        theme_bw() +
        labs(x = "", y = "") +
        facet_wrap(vars(rep_cond), nrow=length(unique(object@conditions))) + 
        xlim(xylim) + ylim(xylim) + 
        labs(title = paste("chromosome:", chr))
    if ((length(unique(interactionsChr$intensity)) > 1) & is.null(trans)==F) {
      p <- p + scale_fill_gradient(low = "white", high = "blue", trans = trans)
    }
    else {
      p <- p + scale_fill_gradient(low = "white", high = "blue")
    }
  } else {
    p <- NULL
  }
  return(p)
}


#' Plot the distance vs intensity matrix.
#'
#' @param objet an \code{HiCDOCExp} object
#' @return A \code{ggplot}.
#' @examples
#' object <- HiCDOCExample()
#' plotDistanceEffect(object)
#' @export
plotDistanceEffect <- function(object) {

  if (is.null(object@interactions)) {
    stop(paste0("Interaction matrix is not loaded yet.  ",
                "Please provide a matrix first."))
  }

    p <- object@interactions %>%
    mutate(distance = position.2 - position.1) %>%
    ggplot(aes(x = distance, y = value)) +
    geom_bin2d() +
    scale_fill_gradient(low = "white", high = "blue", trans = "log2") +
    geom_point(col="transparent") + # necessary for geom_smooth
    geom_smooth(col="red")
  p <- ggMarginal(p, margins = "x", type = "histogram", fill = "transparent")
  return(p)
}

#' Plot the concordance, i.e. the relative distance of the genomic positions with respect to the centroids.
#'
#' @param objet an \code{HiCDOCExp} object
#' @return A list of \code{ggplot}, one for each chromosome.
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotConcordances(object)
#' @export
plotConcordances <- function(object) {

  if (is.null(object@interactions)) {
    stop(paste0("Interaction matrix is not loaded yet.  ",
                "Please provide a matrix first."))
  }
  if (is.null(object@differences)) {
    stop(paste0("Differentially interacting regions are not computed.  ",
                "Please run 'detectCompartmentSwitches' first."))
  }
  if (is.null(object@concordances)) {
    stop(paste0("Concordance is not computed.  ",
                "Please run 'detectCompartments' first."))
  }

  changed <- object@differences %>%
    select(-c(`end`, `padj`)) %>%
    rename(position = start) %>%
    mutate(changed = "T")
  differences <- object@concordances %>%
    group_by(.dots = c("chromosome", "position", "condition")) %>%
    summarise(value = median(value)) %>%
    ungroup() %>%
    spread(condition, value) %>%
    mutate(value = `2` - `1`) %>%
    select(-c(`1`, `2`)) %>%
    left_join(changed, by = c("chromosome", "position")) %>%
    mutate(changed = replace_na(changed, "F")) %>%
    select(-c(chromosome, position))

  p <- ggplot(differences, aes(x = value, fill = changed)) +
       geom_histogram() +
       ggtitle("Distribution of the differences of concordances") +
       xlab("Concordance")

  return(p)
}

.plotAB <- function(data) {
  ggplot(data, aes(x = data$compartment, y = data$diffValue)) +
    geom_jitter(aes(color = data$compartment)) +
    geom_boxplot(outlier.colour = NA, fill = NA, colour = "grey20") +
    labs(color = "Compartment", x = "Compartment", y = "Difference of int.")
}

#' Plot the distribution of A/B compartments along the genomic positions.
#'
#' @param objet an \code{HiCDOCExp} object
#' @return A list of \code{ggplot}, one for each chromosome.
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotAB(object)
#' @export
plotAB <- function(object) {
  p <- buildABComparison(object) %>%
    group_by(chromosome) %>%
    group_split() %>%
    map(.plotAB)
  names(p) <- object@chromosomes
  return(p)
}

.plotCentroids <- function(data) {
  df <- data %>%
    select(-chromosome) %>%
    spread(name, centroid) %>%
    unnest(cols=data$name) %>% t()

  pca <- prcomp(df)
  varpca <- pca$sdev^2
  propvar <- varpca/sum(varpca)
  propvar <- paste(round(100*propvar,2),"%")
  
  pca <- as.data.frame(pca$x)
  pca$group <- row.names(df)
  ggplot(pca, aes(x = PC1, y = PC2, color = group)) + geom_point() +
    xlab(paste("PC1 ", propvar[1])) + 
    ylab(paste("PC2 ", propvar[2]))
}

#' Plot the centroid distributions along the genomic positions.
#'
#' @param objet an \code{HiCDOCExp} object
#' @return A list of \code{ggplot}, one for each chromosome.
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' object <- detectCompartments(object)
#' plotCentroids(object)
#' @export
plotCentroids <- function(object) {

    p <- object@centroids %>%
      unite(name, c(condition, compartment)) %>%
      group_by(chromosome) %>%
      group_split() %>%
      map(.plotCentroids)
    names(p) <- object@chromosomes
    return(p)
}
