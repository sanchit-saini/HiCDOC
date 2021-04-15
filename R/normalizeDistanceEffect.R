#' @description
#' Normalizes the distance effect on the interactions of a given chromosome.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' The name of a chromosome to normalize.
#'
#' @return
#' A tibble of normalized interactions.
#'
#' @keywords internal
#' @noRd
.normalizeDistanceEffectOfChromosome <- function(object, chromosomeName) {

    message("Chromosome ", chromosomeName, ": normalizing distance effect.")

    chromosomeInteractions <-
        object@interactions %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::filter(interaction > 0) %>%
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

#' @title
#' Normalize distance effect.
#'
#' @description
#' Normalizes interactions by their "expected" value relative to the distance
#' that separates their positions. The "expected" values are estimated with a
#' loess regression on the proportion of interactions for each distance.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param loessSampleSize
#' The number of positions used as a sample to estimate the effect of distance
#' on proportion of interactions. Defaults to
#' \code{object$loessSampleSize} which is originally set to
#' \code{defaultHiCDOCParameters$loessSampleSize} = 20000.
#'
#' @return
#' A \code{\link{HiCDOCDataSet}} with normalized interactions.
#'
#' @examples
#' data(exampleHiCDOCDataSet)
#' object <- exampleHiCDOCDataSet
#' object <- normalizeDistanceEffect(object)
#'
#' @seealso
#' \code{\link{normalizeTechnicalBiases}},
#' \code{\link{normalizeBiologicalBiases}},
#' \code{\link{HiCDOC}}
#'
#' @export
normalizeDistanceEffect <- function(object, loessSampleSize = NULL) {

    .validateSlots(
        object,
        slots = c(
            "interactions",
            "chromosomes",
            "binSize",
            "parameters"
        )
    )

    if (!is.null(loessSampleSize)) {
        object@parameters$loessSampleSize <- loessSampleSize
    }
    object@parameters <- .validateParameters(object@parameters)

    normalizedInteractions <-
        purrr::map_dfr(
            object@chromosomes,
            function(chromosomeName) {
                .normalizeDistanceEffectOfChromosome(object, chromosomeName)
            }
        ) %>%
        .sortInteractions(
            object@chromosomes,
            object@conditions,
            object@replicates
        )

    object@interactions <- normalizedInteractions

    return(object)
}
