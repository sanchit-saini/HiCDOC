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
.normalizeDistanceEffectOfChromosome <- function(object) {
    chromosomeName <- as.character(SummarizedExperiment::mcols(object)$Chr[1])
    message("Chromosome ", chromosomeName, ": normalizing distance effect.")
    
    distances <- InteractionSet::pairdist(object, type="mid")
    assay <- SummarizedExperiment::assay(object)
    interaction <- assay
    colnames(interaction) <- paste(object$condition, object$replicate)
    
    # Reordering columns in alphabetic order (useful for tests)
    validAssay <- object@validAssay[[chromosomeName]]
    refOrder <- paste(object$condition, object$replicate)
    interaction <- interaction[, sort(refOrder[validAssay])]
    interaction <- as.vector(interaction)
    chromosomeInteractions <- 
        data.table("distance" = rep(distances, length(validAssay)),
                   "interaction" = interaction)
    chromosomeInteractions <- chromosomeInteractions[!is.na(interaction),]
    idSample <- sample(seq_len(nrow(chromosomeInteractions)),
                   size = min(
                       object@parameters$loessSampleSize,
                       nrow(chromosomeInteractions)
                   ))
    sample <- chromosomeInteractions[idSample]
    setorder(sample, distance)

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
    sample[,loess := stats::predict(loess)]
    sample[,loess := pmax(loess, 0)]
    sample[,interaction := NULL]
    setnames(sample, "distance", "sampleDistance")
    sample <- unique(sample)
    
    uniqueDistances <- unique(sort(chromosomeInteractions$distance))
    sampleDistance <- unique(sort(sample$sampleDistance))
    sampleDistance <- vapply(
        uniqueDistances,
        function(distance) {
            sampleDistance[which.min(abs(distance - sampleDistance))]
        },
        FUN.VALUE = 0
    )
    valueMap <- data.table("distance" = uniqueDistances,
                           "sampleDistance" = sampleDistance)
    valueMap <- merge(valueMap, sample, by="sampleDistance")
    
    loessDistances <- merge(data.table("distance" = distances),
                            valueMap, 
                            by="distance", 
                            sort=FALSE,
                            all.x=T)
    assay <- assay / (loessDistances$loess + 0.00001)
    
    return(assay)
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
normalizeDistanceEffect <- function(object, 
                                    loessSampleSize = NULL, 
                                    parallel=FALSE) {

    .validateSlots(
        object,
        slots = c(
            "chromosomes",
            "parameters"
        )
    )

    if (!is.null(loessSampleSize)) {
        object@parameters$loessSampleSize <- loessSampleSize
    }
    object@parameters <- .validateParameters(object@parameters)
    objectChromosomes <- S4Vectors::split(
        object, 
        SummarizedExperiment::mcols(object)$Chr, drop=FALSE)
    
    normalizedAssays <- .internalApply(parallel,
                                       objectChromosomes,
                                       FUN = .normalizeDistanceEffectOfChromosome) 
    
    normalizedAssays <- do.call("rbind", normalizedAssays)
    SummarizedExperiment::assay(object) <- normalizedAssays
    
    return(object)
}
