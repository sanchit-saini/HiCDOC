
#' @export
plotInteractionMatrix <- function(object, log) {

  if (is.null(object@interactions)) {
    stop(paste0("Interaction matrix is not loaded yet.  ",
                "Please provide a matrix first."))
  }
  fullMatrix <- fullInteractions(object) %>%
    mutate(value = value + 0.0001) %>%
    rename(intensity = value) %>%
    unite("rep_cond", replicate, condition) %>%
    mutate(rep_cond = factor(rep_cond, levels = paste(object@replicates, object@conditions, sep = "_")))
  plots = list()
  for (chr in object@chromosomes) {
    tmp <- fullMatrix %>%
      filter(chromosome == chr)
      # mutate(replicate = factor(mapply(function(X, Y) {
      #                             object@replicates[[X]][Y]
      #                           },
      #                           X = condition,
      #                           Y = replicate)))
    if (nrow(tmp) > 0) {
      p <- tmp %>%
        ggplot(aes(x = position.1, y = position.2, z = intensity)) +
          geom_tile(aes(fill = intensity)) +
          coord_fixed(ratio = 1) +
          theme_bw() +
          labs(x = "", y = "") +
          facet_grid(cols = vars(rep_cond))
      if ((length(unique(tmp$intensity)) > 1) & (log)) {
        p <- p + scale_fill_gradient(low = "white", high = "blue", trans = "log2")
      }
      else {
        p <- p + scale_fill_gradient(low = "white", high = "blue")
      }
      plots[[chr]] <- p
    }
  }
  return(plots)
}

#' @export
plotMD <- function(object) {

  if (is.null(object@interactions)) {
    stop(paste0("Interaction matrix is not loaded yet.  ",
                "Please provide a matrix first."))
  }

  p <- object@interactions %>%
    mutate(distance = position.2 - position.1) %>%
    ggplot(aes(x = distance, y = value)) +
    stat_bin_hex() +
    scale_fill_gradient(low = "white", high = "blue", trans = "log2") +
    geom_smooth()
  p <- ggMarginal(p, margins = "x", type = "histogram", fill = "transparent")
  return(p)
}

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
    unnest() %>% t()
  df
  pca <- prcomp(df)
  pca <- as.data.frame(pca$x)
  pca$group <- row.names(df)
  ggplot(pca, aes(x = PC1, y = PC2, color = group)) + geom_point()
}

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
