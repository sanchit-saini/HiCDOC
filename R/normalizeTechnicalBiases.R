##- normalizeTechnicalBiases -------------------------------------------------#
##----------------------------------------------------------------------------#
#' Normalize the distance effect using a cyclic loess on all the matrices.
#'
#' @rdname normalizeTechnicalBiases
#'
#' @param object A \code{HiCDOCExp} object.
#'
#' @return A \code{HiCDOCExp} object, with the normalized matrices.
#'
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' @export
normalizeTechnicalBiases <- function(object, parallel=FALSE) {

  object@weakBins <- object@weakBins[mixedsort(names(object@weakBins))]
  object@chromosomes <- mixedsort(object@chromosomes)

  # One matrix by condition and replicate
  matrices <- object@interactions %>%
    arrange(order(mixedsort(chromosome))) %>%
    mutate(chromosome = factor(chromosome, levels=object@chromosomes)) %>%
    mutate(chromosome = as.integer(chromosome)) %>%
    group_split(condition, replicate) %>%
    map(function(x) select(x, -c(condition, replicate)))

  # Regions to remove
  remove.regions <- data.frame(cbind(
    unlist(mapply(
      function(bins, name) rep(name, length(bins)),
      object@weakBins,
      as.integer(mixedsort(factor(names(object@weakBins))))
    )),
    (unlist(object@weakBins) - 1) * object@binSize,
    unlist(object@weakBins) * object@binSize - 1
  ))

  # Regions to remove in Granges format
  if (nrow(remove.regions) > 0) {
    colnames(remove.regions) <- c("chromosome", "start", "end")
    remove.regions <- GenomicRanges::makeGRangesFromDataFrame(remove.regions)
  } else {
    remove.regions <- NULL
  }

  hicexp <- make_hicexp(
    data_list = matrices,
    groups = object@conditions,
    remove.regions = remove.regions,
    remove_zeros = FALSE,
    filter = TRUE,
    zero.p = 1,
    A.min = 0
  )

  normalized <- cyclic_loess(hicexp, parallel = parallel)
  output <- hic_table(normalized) %>% as_tibble() %>% select(-D)

  colnames(output) <- c(
    "chromosome", "position.1", "position.2", seq_along(object@replicates)
  )

  object@interactions <- output %>%
    gather(
      as.character(seq_along(object@replicates)),
      key = "i",
      value = "value"
    ) %>%
    mutate(i = factor(as.integer(i))) %>%
    mutate(condition = factor(object@conditions[i])) %>%
    mutate(replicate = factor(object@replicates[i])) %>%
    mutate(chromosome = factor(object@chromosomes[chromosome])) %>%
    select(chromosome, position.1, position.2, condition, replicate, value)

  return(object)
}
