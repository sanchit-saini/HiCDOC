#' @title
#' Normalize technical biases.
#'
#' @description
#' Normalizes technical biases such as sequencing depth by using a cyclic loess
#' to recursively normalize each pair of interaction matrices. Depends on
#' \code{multiHiCcompare}.
#'
#' @details
#' \subsection{Parallel processing}{
#' If \code{parallel=TRUE}, the function
#' \code{\link[multiHiCcompare]{cyclic_loess}}
#' is launched in parallel mode, using \code{\link[BiocParallel]{bplapply}}
#' function. Before to call the function in parallel you should specify
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
#' object <- exampleHiCDOCDataSet
#' object <- filterSparseReplicates(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
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
    hic_table <-
        hic_table[, .(chr = as.numeric(as.factor(as.character(seqnames1))),
                      region1 = start1,
                      region2 = start2)]
    
    currentAssay <- SummarizedExperiment::assay(object)
    currentAssay[is.na(currentAssay)] <- 0
    # Reordering columns in condition order
    refOrder <- paste(object$condition, object$replicate)
    currentAssay <- currentAssay[, order(refOrder)]
    
    table_list <- lapply(seq_len(ncol(currentAssay)), 
                         function(x) cbind(hic_table, currentAssay[,x]))
    
    experiment <-
        multiHiCcompare::make_hicexp(
            data_list = table_list,
            groups = sort(preparedObject$condition),
            remove_zeros = FALSE,
            filter = TRUE,
            zero.p = 1,
            A.min = 0,
            remove.regions = NULL
        )
    normalized <-
        multiHiCcompare::cyclic_loess(experiment, parallel = parallel)
    normalized <- multiHiCcompare::hic_table(normalized)
    
    # Re-sorting the rows in the same order as original
    data.table::setindexv(normalized, c("chr", "region1", "region2"))
    data.table::setindexv(hic_table, c("chr", "region1", "region2"))
    hic_table <- data.table::merge.data.table(hic_table,
                                              normalized,
                                              sort = FALSE)
    
    currentAssay <- as.matrix(hic_table[, 5:ncol(hic_table)])
    # Reordering columns in original order
    currentAssay <- currentAssay[,match(refOrder, sort(refOrder))]
    colnames(currentAssay) <- NULL
    currentAssay[currentAssay == 0] <- NA
    SummarizedExperiment::assay(object) <- currentAssay
    return(object)
}
