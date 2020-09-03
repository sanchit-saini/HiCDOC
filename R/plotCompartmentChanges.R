plotConcordance <- function(concordance,
                            minX    = minX,
                            maxX    = maxX) {
  concordance %>%
    mutate(condition = paste0("confidence, cond. ", condition)) %>%
    ggplot(aes(x = position, y = value, color = replicate)) +
    geom_line() +
    xlim(minX, maxX) +
    geom_hline(yintercept = 0.0, size = 0.1) +
    facet_grid(rows = vars(condition), margins = FALSE, switch = "y") +
    theme_minimal() + 
    theme(axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          strip.placement = "outside"
    )
}

plotCompartments <- function(compartments,
                             binSize = binSize,
                             minX    = minX,
                             maxX    = maxX) {
  compartments %>%
    mutate(compartment = factor(value)) %>%
    #mutate(condition = paste0("cond.", condition)) %>%
    ggplot(aes(x = position, fill = compartment)) +
    geom_histogram(binwidth = binSize) +
    xlim(minX, maxX) +
    facet_grid(rows = vars(condition),
               margins = FALSE, switch = "y")  +
    theme_minimal() + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom" ,
          panel.spacing = unit(0, "cm"),
          strip.placement = "outside"
    )
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
    mutate(padjname = "signif.") %>%
    ggplot(aes(x = start, y = 0, width = binSize)) +
    geom_tile(aes(fill = padj)) +
    facet_grid(rows = vars(padjname), margins = FALSE, switch = "y") +
    xlim(minX, maxX) +
    theme_minimal() + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "bottom",
          strip.placement = "outside"
    )
}

plotChr <- function(chr, object = object) {
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
    plotsgrobs <- lapply( list(pCompartment + theme(legend.position = "none"), 
                               pConcordance  + theme(legend.position = "none")), ggplot2::ggplotGrob)
    commonwidths <- plotsgrobs[[length(plotsgrobs)]]$widths
    plotsgrobs <- lapply(plotsgrobs, function(x) { # Change widths of plots
      x$widths <- commonwidths
      return(x)
    })
    tempplot <- gridExtra::arrangeGrob(plotsgrobs[[1]],
                            plotsgrobs[[2]],
                            ncol    = 1,
                            heights = c(1, 5),
                            top     = paste0("chromosome ",chr),
                            padding = unit(0, "cm"))
  } else {
    plotsgrobs <- lapply( list(pPValue + theme(legend.position = "none"), 
                               pCompartment + theme(legend.position = "none"), 
                               pConcordance  + theme(legend.position = "none")), ggplot2::ggplotGrob)
    commonwidths <- plotsgrobs[[length(plotsgrobs)]]$widths
    plotsgrobs <- lapply(plotsgrobs, 
                         function(x) { # Change widths of plots to align on position
                           x$widths <- commonwidths
                           return(x)
                         }
    )
    tempplot <- gridExtra::arrangeGrob(plotsgrobs[[1]],
                            plotsgrobs[[2]],
                            plotsgrobs[[3]],
                            ncol    = 1,
                            heights = c(1, 1, 5),
                            top     = paste0("chromosome ",chr),
                            padding = unit(0, "cm"))
    
    legPvalue      <- gtable::gtable_filter(ggplotGrob(pPValue), "guide-box", invert=FALSE)
  }
  
  legCompartment <- gtable::gtable_filter(ggplotGrob(pCompartment), "guide-box", invert=FALSE)
  legConcordance <- gtable::gtable_filter(ggplotGrob(pConcordance), "guide-box", invert=FALSE)
  
  legplot <- gridExtra::arrangeGrob(legPvalue, legCompartment, legConcordance,
                         ncol    = 3,  padding = unit(0, "cm"))
  
  finalplot <- gridExtra::arrangeGrob(tempplot, legplot,
                           ncol    = 1,  padding = unit(0, "cm"), 
                           heights = c(10, 1))
  finalplot <- ggpubr::as_ggplot(finalplot)
  return(finalplot)
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
  return(plots)
}

