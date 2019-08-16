#' @export
plotInteractionMatrix <- function(object, log) {

  if (is.null(object@interactionMatrix)) {
    stop(paste0("Interaction matrix is not loaded yet.  ",
                "Please provide a matrix first."))
  }

  fullMatrix <- object@interactionMatrix %>%
    filter(`position 1` != `position 2`) %>%
    rename(tmp = `position 1`,
           `position 1` = `position 2`) %>%
    rename(`position 2` = tmp) %>%
    bind_rows(object@interactionMatrix) %>%
    mutate(value = value + 0.0001) %>%
    rename(intensity = value)
  plots = list()
  for (chr in object@chromosomes) {
    p <- fullMatrix %>%
      filter(chromosome == chr) %>%
      ggplot(aes(x = `position 1`, y = `position 2`, z = intensity)) +
        geom_tile(aes(fill = intensity)) +
        coord_fixed(ratio = 1) +
        theme_bw() +
        labs(x = "", y = "") +
        facet_grid(cols = vars(replicate), rows = vars(chromosome))
    if (log) {
      p <- p + scale_fill_gradient(low = "white", high = "blue", trans = "log2")
    }
    else {
      p <- p + scale_fill_gradient(low = "white", high = "blue")
    }
    plots[[chr]] <- p
  }
  return(plots)
}

#' @export
plotMD <- function(object) {

  if (is.null(object@interactionMatrix)) {
    stop(paste0("Interaction matrix is not loaded yet.  ",
                "Please provide a matrix first."))
  }

  p <- object@interactionMatrix %>%
    mutate(distance = `position 2` - `position 1`) %>%
    ggplot(aes(x = distance, y = value)) +
    stat_bin_hex() +
    scale_fill_gradient(low = "white", high = "blue", trans = "log2") +
    geom_smooth()
  p <- ggMarginal(p, margins = "x", type = "histogram", fill = "transparent")
  return(p)
}

#' @export
plotConcordances <- function(object) {

  if (is.null(object@interactionMatrix)) {
    stop(paste0("Interaction matrix is not loaded yet.  ",
                "Please provide a matrix first."))
  }
  if (is.null(object@DIR)) {
    stop(paste0("Differentially interacting regions are not computed.  ",
                "Please run 'findPValues' first."))
  }
  if (is.null(object@concordances)) {
    stop(paste0("Concordance is not computed.  ",
                "Please run 'detectConstrainedKMeans' first."))
  }

  changed <- object@DIR %>%
    select(-c(`end`, `padj`)) %>%
    rename(position = start) %>%
    mutate(changed = "T")
  differences <- object@concordances %>%
    separate("replicate", c(NA, "condition", "replicate")) %>%
    group_by(.dots = c("chromosome", "position", "condition")) %>%
    summarize(value = median(value)) %>%
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
