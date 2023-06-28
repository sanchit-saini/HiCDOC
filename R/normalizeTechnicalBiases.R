#' @title
#' Normalize technical biases.
#'
#' @description
#' Normalizes technical biases such as sequencing depth by using a cyclic loess
#' to recursively normalize each pair of interaction matrices. Depends on
#' \code{csaw}.
#'
#' @details
#' \subsection{Parallel processing}{
#' If \code{parallel = TRUE}, the normalization is launched in parallel mode
#' by chromosome using \code{\link[BiocParallel]{bplapply}} function. 
#' Before to call the function in parallel you should specify
#' the parallel parameters such as:
#'     \itemize{
#'         \item{On Linux:
#'
#'              \code{multiParam <- BiocParallel::MulticoreParam(workers = 10)}
#'          }
#'          \item{On Windows:
#'
#'              \code{multiParam <- BiocParallel::SnowParam(workers = 10)}
#'         }
#'     }
#'     And then you can register the parameters to be used by BiocParallel:
#'
#'     \code{BiocParallel::register(multiParam, default = TRUE)}
#' }
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @param parallel
#' Whether or not to parallelize the processing. Defaults to FALSE
#'
#' @return
#' A \code{\link{HiCDOCDataSet}} with normalized interactions.
#'
#' @examples
#' data(exampleHiCDOCDataSet)
#' object <- filterSmallChromosomes(exampleHiCDOCDataSet)
#' object <- filterSparseReplicates(object)
#' object <- filterWeakPositions(object)
#' # Not printing loess warnings for example purpose. 
#' # Results should be inspected if there is any.
#' suppressWarnings(
#'     object <- normalizeTechnicalBiases(object)
#' )
#' 
#'
#' @seealso
#' \code{\link{filterSparseReplicates}},
#' \code{\link{filterWeakPositions}},
#' \code{\link{normalizeBiologicalBiases}},
#' \code{\link{normalizeDistanceEffect}},
#' \code{\link{HiCDOC}}
#'
#' @export
normalizeTechnicalBiases <- function(object, parallel = FALSE) {
    message("Normalizing technical biases.")

    hic_table <- as.data.table(InteractionSet::interactions(object))
    hic_table <- hic_table[, .(
        chromosome = seqnames1,
        region1 = start1,
        region2 = start2
    )]
    currentAssay <- SummarizedExperiment::assay(object)
    currentAssay[is.na(currentAssay)] <- 0
    
    # Spliting by chromosome
    listAssays <- S4Vectors::split(
        SummarizedExperiment::assay(object),
        SummarizedExperiment::mcols(object)$chromosome,
        drop = FALSE
    )
    # Keeping only valid columns for each chromosome
    for(c in object@chromosomes){
        listAssays[[c]] <- listAssays[[c]][,object@validAssay[[c]]]
    }
    # For each chromosome, apply a normalisation inter-matrix with csaw
    normedAssay <- .internalLapply(
        parallel,
        listAssays, 
        function(x){
            offsets <- csaw::normOffsets(x, se.out = FALSE, span=0.3, iterations=4L)
            offsets <- offsets - mean(log(colSums(x)))
            return(x / exp(offsets))
        }
    )
    
    # Feeling only valid columns for each chromosome
    for(c in object@chromosomes){
        currentAssay[
            as.logical(SummarizedExperiment::mcols(object)$chromosome == c),
            object@validAssay[[c]]
        ] <- normedAssay[[c]]
    }
    currentAssay[currentAssay == 0] <- NA
    SummarizedExperiment::assay(object, withDimnames = FALSE) <- currentAssay
    return(object)
}
