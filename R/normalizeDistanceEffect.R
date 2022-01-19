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
    
    currentAssay <- SummarizedExperiment::assay(object)
    
    # Reordering columns in alphabetic order (useful for tests)
    validAssay <- object@validAssay[[chromosomeName]]
    refOrder <- paste(object$condition, object$replicate)
    values <- currentAssay[, validAssay]
    values <- values[,order(refOrder[validAssay])]
    
    distances <- InteractionSet::pairdist(object, type="mid")
    chromosomeValues <- 
        data.table("distance" = rep(distances, length(validAssay)),
                   "value" = as.vector(values))
    chromosomeValues <- chromosomeValues[!is.na(value),]
    idSample <- sample(seq_len(nrow(chromosomeValues)),
                   size = min(
                       object@parameters$loessSampleSize,
                       nrow(chromosomeValues)
                   ))
    sample <- chromosomeValues[idSample]
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
            value ~ distance,
            data = sample,
            control = stats::loess.control(trace.hat = traceMethod)
        )
    span <- optimizeSpan(loess, criterion = "gcv")

    loess <-
        stats::loess(
            value ~ distance,
            span = span,
            data = sample,
            control = stats::loess.control(trace.hat = traceMethod)
        )
    sample[,loess := stats::predict(loess)]
    sample[,loess := pmax(loess, 0)]
    sample[,value := NULL]
    data.table::setnames(sample, "distance", "sampleDistance")
    sample <- unique(sample)
    
    uniqueDistances <- unique(sort(chromosomeValues$distance))
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
    valueMap <- data.table::merge.data.table(valueMap,
                                             sample,
                                             by="sampleDistance")
    
    loessDistances <- data.table::merge.data.table(
        data.table("distance" = distances),
        valueMap, 
        by="distance", 
        sort=FALSE,
        all.x=T)
    currentAssay <- currentAssay / (loessDistances$loess + 0.00001)
    
    return(currentAssay)
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
    if (!is.null(loessSampleSize)) {
        object@parameters$loessSampleSize <- loessSampleSize
    }
    object@parameters <- .validateParameters(object@parameters)
    objectChromosomes <- S4Vectors::split(
        object, 
        SummarizedExperiment::mcols(object)$Chr, drop=FALSE)
    
    normAssay <- .internalApply(parallel,
                                objectChromosomes,
                                FUN = .normalizeDistanceEffectOfChromosome,
                                type="lapply") 
    
    normAssay <- do.call("rbind", normAssay)
    SummarizedExperiment::assay(object) <- normAssay
    
    return(object)
}
