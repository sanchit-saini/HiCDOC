plotConcordance <- function(concordance,
                            minX    = minX,
                            maxX    = maxX) {
    concordance %>%
      ggplot(aes(x = position, y = value, color = replicate)) +
      geom_line() +
      xlim(minX, maxX) +
      geom_hline(yintercept = 0.0, size = 0.1) +
      facet_grid(rows = vars(condition), margins = FALSE) +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = "bottom",
            strip.text.y = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            panel.spacing = unit(0, "cm"))
}

plotCompartments <- function(compartments,
                             binSize = binSize,
                             minX    = minX,
                             maxX    = maxX) {
    compartments %>%
      mutate(compartment = factor(value)) %>%
      ggplot(aes(x = position, fill = compartment)) +
      geom_histogram(binwidth = binSize) +
      xlim(minX, maxX) +
      facet_grid(rows = vars(condition),
                 margins = FALSE) +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = "none",
            strip.text.y = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            panel.spacing = unit(0, "cm"))
}

plotPValue <- function(differences,
                       binSize = binSize,
                       minX    = minX,
                       maxX    = maxX) {
    if (nrow(differences) == 0) {
      return(NULL)
    }
    differences %>%
      mutate(padj = -sign(padj) * log10(abs(padj))) %>%
      ggplot(aes(x = start,
                 y = 0,
                 width = binSize / 2)) +
      geom_tile(aes(fill = padj)) +
      xlim(minX, maxX) +
      scale_fill_gradient2() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = "top",
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            panel.spacing = unit(0, "cm"))
}

plotChr <- function(chr, object = object) {
  theme_set(theme_minimal())
  maxX <- object@concordances %>%
    filter(chromosome == chr) %>%
    summarise(max = max(position)) %>%
    pull()
  minX <- -maxX / 100.0
  maxX <- 101 / 100 * maxX

  pConcordance <- plotConcordance(object@concordances %>%
                                    filter(chromosome == chr),
                                  minX = minX,
                                  maxX = maxX)
  pCompartment <- plotCompartments(object@compartments %>%
                                     filter(chromosome == chr),
                                   binSize = object@binSize,
                                   minX    = minX,
                                   maxX    = maxX)
  pPValue <- plotPValue(object@differences %>%
                          filter(chromosome == chr),
                        binSize = object@binSize,
                        minX    = minX,
                        maxX    = maxX)
  if (is.null(pPValue)) {
    return(gridExtra::arrangeGrob(pCompartment,
                       pConcordance,
                       ncol    = 1,
                       heights = c( 1, 5),
                       top     = chr,
                       padding = unit(0, "cm")))
  }
  return(gridExtra::arrangeGrob(pPValue,
                     pCompartment,
                     pConcordance,
                     ncol    = 1,
                     heights = c(1, 1, 5),
                     top     = chr,
                     padding = unit(0, "cm")))
}

#' @export
plotCompartmentChanges <- function(object) {

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
  if (is.null(object@compartments)) {
    stop(paste0("Compartments are not computed.  ",
                "Please run 'detectCompartments' first."))
  }

  plots <- lapply(object@chromosomes, plotChr, object = object)
  #plots <- marrangeGrob(lapply(object@chromosomes, plotChr, object = object), ncol = 1, nrow = length(object@chromosomes))
  return(plots)
}
