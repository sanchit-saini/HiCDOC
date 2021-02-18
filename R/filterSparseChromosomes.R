#' sparsityChromosome
#' Compute the maximum of sparsity of a given chromosome
#'
#' @param object a HiCDOCDataSet object
#' @param chromosomeId the name of the chromosome
#'
#' @return a one line tibble
#' @keywords internal
sparsityChromosome <-
    function(object,
    chromosomeId) {
        interactionsChr <- object@interactions %>%
            dplyr::filter(chromosome == chromosomeId) %>%
            dplyr::filter(value > 0)
        totalCells <- object@totalBins[chromosomeId]^2

        pctFill <- interactionsChr %>%
            dplyr::mutate(value = ifelse(bin.1 == bin.2, 1, 2)) %>%
            dplyr::group_by(replicate, condition) %>%
            dplyr::summarise(
                pctSparse = 1 - sum(value) / totalCells,
                .groups = "keep"
            ) %>%
            dplyr::ungroup()

        pctFill %<>%
            dplyr::summarise(maxSparsity = max(pctSparse, na.rm = TRUE)) %>%
            dplyr::mutate(
                totalBins = object@totalBins[chromosomeId],
                totalCells = totalCells,
                chromosome = factor(chromosomeId, levels = object@chromosomes)
            ) %>%
            dplyr::mutate(quality = ifelse(
                maxSparsity < 0.001,
                "***",
                ifelse(
                    maxSparsity < 0.01,
                    "**",
                    ifelse(maxSparsity < 0.05, "*",
                        ifelse(maxSparsity < 0.1, ".", "")
                    )
                )
            ))

        # reordering columns
        pctFill %<>% dplyr::select(chromosome, totalCells, maxSparsity, quality)
        return(pctFill)
    }

#' #filterSparseChromosomes
#' The function compute the sparsity of the matrix, by condition and
#' replicate, for each chromosome. Then it remove (if \code{removeChromosomes}
#' is TRUE) the chromosomes with a sparcity >= \code{thresholdSparseMatrix}
#' on at least one condition and one replicate.
#' @param object A HiCDOCDataSet object
#' @param thresholdSparseMatrix A numerical value. The sparsity threshold
#' from which a chromosome will be removed, if \code{removeChromosomes} is TRUE.
#' If NULL, default to the first not NULL of \code{object$sparseThreshold.} and
#' \code{HiCDOCDefaultParameters$sparseThreshold.}.
#' @param removeChromosomes Logical, default to TRUE. Should the function
#' remove sparse chromosomes ?
#'
#' @return A HiCDOCDataSet object
#' @export
#' @seealso \code{\link[HiCDOC]{filterSmallChromosomes}}, 
#' \code{\link[HiCDOC]{filterWeakPositions}} and 
#' \code{\link[HiCDOC]{runHiCDOC}} for the recommended pipeline.
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSparseChromosomes(object)
filterSparseChromosomes <-
    function(object,
    thresholdSparseMatrix = NULL,
    removeChromosomes = TRUE) {
        # Parameters
        if (!is.null(thresholdSparseMatrix)) {
            object@parameters$sparseThreshold <- thresholdSparseMatrix
        }
        object@parameters <- checkParameters(object@parameters)
        thresh <- object@parameters$sparseThreshold
        # Sparsity
        chrQuality <-
            purrr::map_dfr(
                object@chromosomes,
                function(x) {
                      sparsityChromosome(
                          object,
                          x
                      )
                  }
            )
        # Remove chromosomes
        if (removeChromosomes == TRUE) {
            sparseChromosomes <- chrQuality %>%
                dplyr::filter(maxSparsity >= thresh) %>%
                dplyr::pull(chromosome) %>%
                as.character()
            if (length(sparseChromosomes) == 0) {
                message(
                    "No chromosome removed (threshold: ",
                    round(100 * thresh, 2),
                    "%)"
                )
            } else {
                message(
                    length(sparseChromosomes),
                    " chromosome(s) removed: ",
                    sparseChromosomes,
                    " (threshold : ",
                    round(100 * thresh, 2),
                    "%)"
                )
                nonSparseChr <-
                    object@chromosomes[!(object@chromosomes %in%
                        sparseChromosomes)]
                object <-
                    reduceHiCDOCDataSet(object, chromosomes = nonSparseChr)
            }
        }
        chrQuality %<>% dplyr::mutate(
            maxSparsity =
                paste0(round(100 * maxSparsity, 2), "%")
        )
        print(data.frame(chrQuality)) # In data.frame for correct alignement
        cat("\n")
        return(object)
    }
