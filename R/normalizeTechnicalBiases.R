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
    .validateSlots(object,
                   slots = c("chromosomes"))
    
    message("Normalizing technical biases.")
    
    hic_table <- as.data.table(InteractionSet::interactions(object))
    hic_table <-
        hic_table[, .(chr = as.numeric(as.factor(as.character(seqnames1))),
                      region1 = start1,
                      region2 = start2)]
    hic_table[, D := InteractionSet::pairdist(object, type = "diag")]
    matAssay <- SummarizedExperiment::assay(object)
    matAssay[is.na(matAssay)] <- 0
    # Reordering columns in alphabetic order (useful for tests)
    refOrder <- paste(object$condition, object$replicate)
    matAssay <- matAssay[, order(refOrder)]
    matAssay <- as.data.table(matAssay)
    setnames(matAssay, paste0("IF", seq_len(ncol(matAssay))))
    
    hic_table <- cbind(hic_table, matAssay)
    
    table_list <- split(hic_table, hic_table$chr)
    
    # plug into parallelized loess function
    normalized <- .internalApply(
        parallel = parallel,
        table_list,
        FUN = .cloess,
        iterations = 3,
        verbose = FALSE,
        span = NA
    )
    normalized <- data.table::rbindlist(normalized)
    
    ifcolumns <- colnames(hic_table)[4:ncol(hic_table)]
    hic_table[, (ifcolumns) := NULL]
    data.table::setindexv(normalized, c("chr", "region1", "region2"))
    data.table::setindexv(hic_table, c("chr", "region1", "region2"))
    hic_table <- merge(hic_table, normalized, sort = FALSE)
    
    matAssay <- hic_table[, 5:ncol(hic_table)]
    matAssay <- as.matrix(matAssay)
    # Reordering columns in original order
    matAssay <- matAssay[,match(refOrder, sort(refOrder))]
    colnames(matAssay) <- NULL
    matAssay[matAssay == 0] <- NA
    SummarizedExperiment::assay(object) <- matAssay
    return(object)
}
