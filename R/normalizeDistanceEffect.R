##- normalizeDistanceEffectChr --------------------------------------------------#
##----------------------------------------------------------------------------#
#' Normalize the distance effect using a loess on the intensity vs distance
#' to diagonal.
#'
#' @rdname normalizeDistanceEffectChr
#'
#' @param object A \code{HiCDOCExp} object.
#' @param chromosomeId The name or number of the chromosome to plot.
#' If number, will be taken in \code{object@chromosomes[chromosomeId]}
#'
#' @return the normalized interaction matrix for this chromosome.
normalizeDistanceEffectChr <- function(object, chromosomeId) {
  testSlotsHiCDOCExp(object,
                     slots = c("interactions"))
  chr <- testchromosome(object, chromosomeId)
  
  message("Chromosome: ", chr)

  interactionsChr <- object@interactions %>%
    filter(chromosome == chr & value>0) %>%
    mutate(distance = position.2 - position.1) %>%
    mutate(
      bin.1 = position.1 / object@binSize + 1,
      bin.2 = position.2 / object@binSize + 1
    )
  
  sample <- interactionsChr %>%
    filter(!(bin.1 %in% object@weakBins[[chromosomeId]])) %>%
    filter(!(bin.2 %in% object@weakBins[[chromosomeId]])) %>%
    sample_n(size = min(object@sampleSize, nrow(interactionsChr))) %>%
    select(distance, value) %>%
    arrange(distance)
  
  if (nrow(sample) == 0) {
    message("The chromosome is empty")
    next
  }
  
  optimizeSpan <- function(
    model,
    criterion = c("aicc", "gcv"),
    spans = c(0.01, 0.9)
  ) {
    criterion <- match.arg(criterion)
    f <- function(span) {
      m      <- update(model, span = span)
      span   <- m$pars$span
      traceL <- m$trace.hat
      sigma2 <- sum(m$residuals^2)/(m$n - 1)
      if (criterion == "aicc") {
        quality <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(m$n - traceL - 2)
      } else if (criterion == "gcv") {
        quality <- m$n * sigma2/(m$n - traceL)^2
      }
      return (quality)
    }
    result <- optimize(f, spans)
    return (result$minimum)
  }
  
  methodtrace <- "approximate"
  if(object@sampleSize <= 1000)
    methodtrace <- "exact"
  
  l <- loess(value ~ distance, data = sample,
             control = loess.control(trace.hat=methodtrace))
  span <- optimizeSpan(l, criterion = "gcv")
  l <- loess(value ~ distance, span = span, data = sample,
             control = loess.control(trace.hat=methodtrace))
  
  sample %<>%
    mutate(loess = predict(l)) %>%
    mutate(loess = pmax(loess, 0)) %>%
    rename(sampleDistance = distance) %>%
    select(-value) %>%
    unique()
  
  sampleDistances <- unique(sort(sample$sampleDistance))
  uniqueDistances <- unique(sort(interactionsChr$distance))
  valueMap <- tibble(
    distance = uniqueDistances,
    sampleDistance = vapply(uniqueDistances, function(x) {
      sampleDistances[which.min(abs(x - sampleDistances))]
    }, FUN.VALUE = c(0))
  ) %>%
    left_join(sample, by = "sampleDistance") %>%
    select(-sampleDistance)
  
  interactionsChr %<>%
    left_join(valueMap, by = "distance") %>%
    mutate(value = value / (loess + 0.00001)) %>%
    select(-distance, -loess) %>%
    select(-c(bin.1, bin.2))
  
  return(interactionsChr)
}



##- normalizeDistanceEffect --------------------------------------------------#
##----------------------------------------------------------------------------#
#' Normalize the distance effect using a loess on the intensity vs distance
#' to diagonal.
#'
#' @rdname normalizeDistanceEffect
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
#' object <- normalizeBiologicalBiases(object)
#' object <- normalizeDistanceEffect(object)
#' @export
normalizeDistanceEffect <- function(object) {
  
  interactionsNorm <- purrr::map_dfr(object@chromosomes,
                                     function(x) normalizeDistanceEffectChr(object, x)) %>%
    mutate(chromosome = factor(chromosome, levels = object@chromosomes))
  object@interactions <- interactionsNorm
  
  return(object)
}
