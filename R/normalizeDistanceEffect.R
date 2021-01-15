##- normalizeDistanceEffectChr --------------------------------------------------#
##----------------------------------------------------------------------------#
#' Normalize the distance effect using a loess on the intensity vs distance
#' to diagonal.
#'
#' @rdname normalizeDistanceEffectChr
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param chromosomeId The name or number of the chromosome to plot.
#' If number, will be taken in \code{object@chromosomes[chromosomeId]}
#'
#' @return the normalized interaction matrix for this chromosome.
normalizeDistanceEffectChr <- function(object, chromosomeId) {
    testSlotsHiCDOC(object,
                       slots = c("interactions",
                                 "binSize",
                                 "weakBins"))
    chr <- testChromosome(object, chromosomeId)

    object@parameters <- checkParameters(object@parameters,
                                         c("sampleSize"))
    message("Chromosome: ", chr)

    interactionsChr <- object@interactions %>%
        dplyr::filter(chromosome == chr & value > 0) %>%
        dplyr::mutate(distance = (bin.2 - bin.1) * object@binSize)

    sample <- interactionsChr %>%
        dplyr::sample_n(size = min(object@parameters$sampleSize,
                                   nrow(interactionsChr))) %>%
        dplyr::select(distance, value) %>%
        dplyr::arrange(distance)

    if (nrow(sample) == 0) {
        message("The chromosome is empty")
        return(NULL)
    }

    optimizeSpan <- function(model,
                             criterion = c("aicc", "gcv"),
                             spans = c(0.01, 0.9)) {
        criterion <- match.arg(criterion)
        f <- function(span) {
            m <- update(model, span = span)
            span <- m$pars$span
            traceL <- m$trace.hat
            sigma2 <- sum(m$residuals ^ 2) / (m$n - 1)
            if (criterion == "aicc") {
                quality <-
                    log(sigma2) + 1 + 2 * (2 * (traceL + 1)) / (m$n - traceL - 2)
            } else if (criterion == "gcv") {
                quality <- m$n * sigma2 / (m$n - traceL) ^ 2
            }
            return (quality)
        }
        result <- optimize(f, spans)
        return (result$minimum)
    }

    methodtrace <- "approximate"
    if (object@parameters$sampleSize <= 1000)
        methodtrace <- "exact"

    l <- loess(value ~ distance,
               data = sample,
               control = loess.control(trace.hat = methodtrace))
    span <- optimizeSpan(l, criterion = "gcv")

    l <- loess(
        value ~ distance,
        span = span,
        data = sample,
        control = loess.control(trace.hat = methodtrace)
    )

    sample %<>%
        dplyr::mutate(loess = predict(l)) %>%
        dplyr::mutate(loess = pmax(loess, 0)) %>%
        dplyr::rename(sampleDistance = distance) %>%
        dplyr::select(-value) %>%
        unique()

    sampleDistances <- unique(sort(sample$sampleDistance))
    uniqueDistances <- unique(sort(interactionsChr$distance))
    valueMap <- dplyr::tibble(
        distance = uniqueDistances,
        sampleDistance = vapply(uniqueDistances, function(x) {
            sampleDistances[which.min(abs(x - sampleDistances))]
        }, FUN.VALUE = c(0))
    ) %>%
        dplyr::left_join(sample, by = "sampleDistance") %>%
        dplyr::select(-sampleDistance)

    interactionsChr %<>%
        dplyr::left_join(valueMap, by = "distance") %>%
        dplyr::mutate(value = value / (loess + 0.00001)) %>%
        dplyr::select(-distance,-loess)

    return(interactionsChr)
}



##- normalizeDistanceEffect --------------------------------------------------#
##----------------------------------------------------------------------------#
#' Normalize the distance effect using a loess on the intensity vs distance
#' to diagonal.
#'
#' @rdname normalizeDistanceEffect
#'
#' @param object A \code{HiCDOCDataSet} object.
#'
#' @return A \code{HiCDOCDataSet} object, with the normalized matrices.
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
    interactionsNorm <-
        purrr::map_dfr(object@chromosomes,
                       function(x) normalizeDistanceEffectChr(object, x)) %>%
        dplyr::mutate(chromosome = factor(chromosome, levels = object@chromosomes))
    object@interactions <- interactionsNorm

    return(object)
}
