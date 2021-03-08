## - .normalizeDistanceEffectOfChromosome ---------------------------------------------#
## --------------------------------------------------------------------------#
#' Normalize the distance effect using a loess on the intensity vs distance
#' to diagonal.
#' @param object A \code{HiCDOCDataSet} object.
#' @param chromosomeName The name or number of the chromosome to plot.
#' If number, will be taken in \code{object@chromosomes[chromosomeName]}
#'
#' @return the normalized interaction matrix for this chromosome.
#' @keywords internal
#' @noRd
.normalizeDistanceEffectOfChromosome <- function(object, chromosomeName) {
    .validateSlots(
        object,
        slots = c(
            "interactions",
            "binSize",
            "weakBins",
            "parameters"
        )
    )

    message("Chromosome ", chromosomeName, ": normalizing distance effect.")

    chromosomeInteractions <-
        object@interactions %>%
        dplyr::filter(chromosome == chromosomeName & interaction > 0) %>%
        dplyr::mutate(distance = (bin.2 - bin.1) * object@binSize)

    sample <-
        chromosomeInteractions %>%
        dplyr::sample_n(
            size = min(
                object@parameters$loessSampleSize,
                nrow(chromosomeInteractions)
            )
        ) %>%
        dplyr::select(distance, interaction) %>%
        dplyr::arrange(distance)

    if (nrow(sample) == 0) {
        message("Chromosome ", chromosomeName, " is empty.")
        return(NULL)
    }

    optimizeSpan <- function(
        model,
        criterion = c("aicc", "gcv"),
        spans = c(0.01, 0.9)
    ) {
        criterion <- match.arg(criterion)
        result <- stats::optimize(
            function(span) {
                model <- stats::update(model, span = span)
                span <- model$pars$span
                trace <- model$trace.hat
                sigma2 <- sum(model$residuals^2) / (model$n - 1)
                if (criterion == "aicc") {
                    quality <-
                        log(sigma2) + 1 + 2 * (2 * (trace + 1)) /
                        (model$n - trace - 2)
                } else if (criterion == "gcv") {
                    quality <- model$n * sigma2 / (model$n - trace)^2
                }
                return(quality)
            },
            spans
        )
        return(result$minimum)
    }

    traceMethod <- "approximate"
    if (object@parameters$loessSampleSize <= 1000) traceMethod <- "exact"

    loess <-
        stats::loess(
            interaction ~ distance,
            data = sample,
            control = stats::loess.control(trace.hat = traceMethod)
        )
    span <- optimizeSpan(loess, criterion = "gcv")

    loess <-
        stats::loess(
            interaction ~ distance,
            span = span,
            data = sample,
            control = stats::loess.control(trace.hat = traceMethod)
        )

    sample %<>%
        dplyr::mutate(loess = stats::predict(loess)) %>%
        dplyr::mutate(loess = pmax(loess, 0)) %>%
        dplyr::rename(sampleDistance = distance) %>%
        dplyr::select(-interaction) %>%
        unique()

    sampleDistances <- unique(sort(sample$sampleDistance))
    uniqueDistances <- unique(sort(chromosomeInteractions$distance))
    valueMap <-
        dplyr::tibble(
            distance = uniqueDistances,
            sampleDistance = vapply(
                uniqueDistances,
                function(distance) {
                    sampleDistances[which.min(abs(distance - sampleDistances))]
                },
                FUN.VALUE = 0
            )
        ) %>%
        dplyr::left_join(sample, by = "sampleDistance") %>%
        dplyr::select(-sampleDistance)

    chromosomeInteractions %<>%
        dplyr::left_join(valueMap, by = "distance") %>%
        dplyr::mutate(interaction = interaction / (loess + 0.00001)) %>%
        dplyr::select(-distance, -loess)

    return(chromosomeInteractions)
}



## - normalizeDistanceEffect ------------------------------------------------#
## --------------------------------------------------------------------------#
#' Normalize the distance effect on a HiCDOCDataSet object using a loess
#' on the intensity vs distance to diagonal.
#'
#' @name normalizeDistanceEffect
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param loessSampleSize A numerical value. The number of bins used when sampling
#' all the bins. If NULL, default to the first not NULL of
#' \code{object$loessSampleSize} and \code{defaultHiCDOCParameters$loessSampleSize}.
#'
#' @return A \code{HiCDOCDataSet} object, with the normalized matrices.
#'
#' @examples
#' object <- HiCDOCDataSetExample()
#' object <- normalizeDistanceEffect(object)
#' @seealso \code{\link[HiCDOC]{normalizeTechnicalBiases}},
#' \code{\link[HiCDOC]{normalizeBiologicalBiases}} and
#' \code{\link[HiCDOC]{HiCDOC}} for the recommended pipeline.
#' @export
normalizeDistanceEffect <- function(object, loessSampleSize = NULL) {
    # Parameters
    if (!is.null(loessSampleSize)) {
        object@parameters$loessSampleSize <- loessSampleSize
    }
    object@parameters <- .validateParameters(object@parameters)
    # Normalization by chromosome
    normalizedInteractions <-
        purrr::map_dfr(
            object@chromosomes,
            function(chromosomeName) {
                .normalizeDistanceEffectOfChromosome(object, chromosomeName)
            }
        ) %>%
        dplyr::mutate(
            chromosome = factor(chromosome, levels = object@chromosomes)
        )
    object@interactions <- normalizedInteractions

    return(object)
}
