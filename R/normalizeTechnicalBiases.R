## - normalizeTechnicalBiases -----------------------------------------------#
## --------------------------------------------------------------------------#
#' Normalize the distance effect using a cyclic loess on all the matrices.
#'
#' @rdname normalizeTechnicalBiases
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param parallel Logical, defautl to FALSE. Should parallel computing be
#' used in \code{\link[multiHiCcompare]{cyclic_loess}} ?
#'
#' @return A \code{HiCDOCDataSet} object, with the normalized matrices.
#'
#' @examples
#' object <- HiCDOCDataSetExample()
#' object <- normalizeTechnicalBiases(object)
#' @seealso \code{\link[HiCDOC]{normalizeBiologicalBiases}},
#' \code{\link[HiCDOC]{normalizeDistanceEffect}} and
#' \code{\link[HiCDOC]{HiCDOC}} for the recommended pipeline.
#' @export
normalizeTechnicalBiases <- function(object, parallel = FALSE) {

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
        dplyr::group_split(condition, replicate) %>%
        purrr::map(function(x) dplyr::select(x, -c(condition, replicate)))

    # Regions to remove
    remove.regions <-
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

    # Regions to remove in Granges format
    if (nrow(remove.regions) > 0) {
        remove.regions %<>%
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
        remove.regions <- NULL
    }

    hicexp <-
        multiHiCcompare::make_hicexp(
            data_list = matrices,
            groups = object@conditions,
            remove.regions = remove.regions,
            remove_zeros = FALSE,
            filter = TRUE,
            zero.p = 1,
            A.min = 0
        )

    normalized <- multiHiCcompare::cyclic_loess(hicexp, parallel = parallel)
    result <-
        multiHiCcompare::hic_table(normalized) %>%
        dplyr::as_tibble() %>%
        dplyr::select(-D)

    colnames(result) <-
        c("chromosome", "bin.1", "bin.2", seq_along(object@replicates))
    result %<>%
        dplyr::mutate(
            bin.1 = as.integer(bin.1 / object@binSize + 1),
            bin.2 = as.integer(bin.2 / object@binSize + 1)
        )

    object@interactions <-
        result %>%
        tidyr::gather(
            as.character(seq_along(object@replicates)),
            key = "index",
            value = "interaction"
        ) %>%
        dplyr::mutate(index = factor(as.integer(index))) %>%
        dplyr::mutate(condition = factor(object@conditions[index])) %>%
        dplyr::mutate(replicate = factor(object@replicates[index])) %>%
        dplyr::mutate(chromosome = factor(object@chromosomes[chromosome])) %>%
        dplyr::select(
            chromosome,
            condition,
            replicate,
            bin.1,
            bin.2,
            interaction
        )

    return(object)
}
