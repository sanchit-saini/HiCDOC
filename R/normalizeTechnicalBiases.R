#' @title
#' Normalize technical biases.
#'
#' @description
#' Normalizes technical biases such as sequencing depth by using a cyclic loess
#' to recursively normalize each pair of interaction matrices. Depends on
#' \code{multiHiCcompare}.
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
#'              \code{BiocParallel::register(multiParam, default = TRUE)}
#'          }
#'          \item{On Windows:
#'
#'              \code{multiParam <- BiocParallel::SnowParam(workers = 10)}
#'              \code{BiocParallel::register(multiParam, default = TRUE)}
#'         }
#'     }
#' }
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @param parallel
#' Whether or not to parallelize the processing. Defaults to TRUE.
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
normalizeTechnicalBiases <- function(object, parallel = TRUE) {

    .validateSlots(
        object,
        slots = c(
            "interactions",
            "chromosomes",
            "binSize",
            "weakBins",
            "validReplicates",
            "validConditions"
        )
    )

    message("Normalizing technical biases.")

    # One matrix per condition and replicate
    matrices <-
        object@interactions %>%
        dplyr::mutate(
            chromosome = as.integer(factor(
                chromosome,
                levels = object@chromosomes
            )),
            bin.1 = (bin.1 - 1) * object@binSize,
            bin.2 = (bin.2 - 1) * object@binSize
        ) %>%
        dplyr::group_split(condition, replicate)

    groups <-
        matrices %>%
        purrr::map(
            function(group) {
                dplyr::select(group, c(condition, replicate)) %>%
                dplyr::slice(1)
            }
        ) %>%
        purrr::reduce(rbind)

    matrices %<>%
        purrr::map(
            function(group) dplyr::select(group, -condition, -replicate)
        )

    # Regions to ignore during normalization
    weakRegions <-
        data.frame(
            "chromosome" = unlist(mapply(
                function(bins, chromosomeName) {
                    rep(chromosomeName, length(bins))
                },
                object@weakBins,
                names(object@weakBins)
            )),
            "bin" = unlist(object@weakBins)
        )

    if (nrow(weakRegions) > 0) {

        weakRegions %<>%
            dplyr::mutate(
                start = (bin - 1) * object@binSize,
                end = (bin * object@binSize) - 1
            ) %>%
            dplyr::select(-bin) %>%
            dplyr::mutate(
                chromosome = as.integer(factor(
                    chromosome,
                    levels = object@chromosomes
                ))
            ) %>%
            GenomicRanges::makeGRangesFromDataFrame()

    } else {
        weakRegions <- NULL
    }

    # Cyclic loess normalization
    experiment <-
        multiHiCcompare::make_hicexp(
            data_list = matrices,
            groups = groups$condition,
            remove.regions = weakRegions,
            remove_zeros = FALSE,
            filter = TRUE,
            zero.p = 1,
            A.min = 0
        )
    normalized <- multiHiCcompare::cyclic_loess(experiment, parallel = parallel)

    result <-
        multiHiCcompare::hic_table(normalized) %>%
        dplyr::as_tibble() %>%
        dplyr::select(-D)
    colnames(result) <-
        c("chromosome", "bin.1", "bin.2", seq_along(groups$replicate))
    result %<>%
        dplyr::mutate(
            bin.1 = as.integer(bin.1 / object@binSize + 1),
            bin.2 = as.integer(bin.2 / object@binSize + 1)
        )

    object@interactions <-
        result %>%
        tidyr::gather(
            as.character(seq_along(groups$replicate)),
            key = "index",
            value = "interaction"
        ) %>%
        dplyr::mutate(index = factor(as.integer(index))) %>%
        dplyr::mutate(condition = groups$condition[index]) %>%
        dplyr::mutate(replicate = groups$replicate[index]) %>%
        dplyr::mutate(chromosome = object@chromosomes[chromosome]) %>%
        dplyr::select(
            chromosome,
            condition,
            replicate,
            bin.1,
            bin.2,
            interaction
        ) %>%
        .sortInteractions(
            object@chromosomes,
            object@conditions,
            object@replicates
        )

    return(object)
}
