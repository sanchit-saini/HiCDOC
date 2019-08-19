#' @export
plotCompartmentChanges <- function(object) {

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
  if (is.null(object@compartments)) {
    stop(paste0("Compartments are not computed.  ",
                "Please run 'detectConstrainedKMeans' first."))
  }

  theme_set(theme_minimal())

  plots <- lapply(object@chromosomes,
                  function(chr) {

                    maxX = unlist(object@concordances %>%
                        filter(chromosome == chr) %>%
                        summarize(max = max(position)))[[1]]
                    minX = -maxX / 100.0
                    maxX = maxX + (maxX / 100.0)


                    pConcordance <- object@concordances %>%
                        separate("replicate",
                                 c(NA, "condition", "replicate")) %>%
                        filter(chromosome == chr) %>%
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

                    pCompartment <- object@compartments %>%
                        filter(chromosome == chr) %>%
                        mutate(compartment = factor(value)) %>%
                        ggplot(aes(x = position, fill = compartment)) +
                          geom_histogram(binwidth = object@binSize) +
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

                    pPvalue <- object@DIR %>%
                        filter(chromosome == chr) %>%
                        mutate(padj = -sign(padj) * log10(abs(padj))) %>%
                        ggplot(aes(x = start,
                                   y = 0,
                                   width = object@binSize / 2)) +
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

                    grid.arrange(pPvalue,
                                 pCompartment,
                                 pConcordance,
                                 ncol = 1,
                                 heights = c(1, 1, 5),
                                 top = chr,
                                 padding = unit(0, "cm"))
                    # plot_grid(pPvalue,
                    #           pCompartment,
                    #           pConcordance,
                    #           ncol = 1,
                    #           rel_heights = c(1, 1, 5))
                  })

  names(plots) <- object@chromosomes

  return(plots)
}


