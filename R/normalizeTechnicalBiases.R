##- normalizeTechnicalBiases -------------------------------------------------#
##----------------------------------------------------------------------------#
#' Normalize the distance effect using a cyclic loess on all the matrices.
#'
#' @rdname normalizeTechnicalBiases
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param parallel Logical, defautl to FALSE. Should parallel computing be used
#' in \code{\link[multiHiCcompare]{cyclic_loess}} ?
#'
#' @return A \code{HiCDOCDataSet} object, with the normalized matrices.
#'
#' @examples
#' object <- HiCDOCExample()
#' object <- filterSmallChromosomes(object)
#' object <- filterWeakPositions(object)
#' object <- normalizeTechnicalBiases(object)
#' @export
normalizeTechnicalBiases <- function(object, parallel=FALSE) {

    object@weakBins <- object@weakBins[mixedsort(names(object@weakBins))]

    # One matrix by condition and replicate
    matrices <- object@interactions %>%
        arrange(order(mixedsort(chromosome))) %>%
        dplyr::mutate(chromosome = factor(chromosome, levels=object@chromosomes)) %>%
        dplyr::mutate(chromosome = as.integer(chromosome)) %>%
        dplyr::group_split(condition, replicate) %>%
        purrr::map(function(x) dplyr::select(x, -c(condition, replicate)))

    # Regions to remove
    remove.regions <- data.frame(cbind(
        unlist(mapply(
            function(bins, name) rep(name, length(bins)),
            object@weakBins,
            as.integer(mixedsort(factor(names(object@weakBins))))
        )),
        (unlist(object@weakBins) - 1) * object@binSize,
        unlist(object@weakBins) * object@binSize - 1
    ))

    # Regions to remove in Granges format
    if (nrow(remove.regions) > 0) {
        colnames(remove.regions) <- c("chromosome", "start", "end")
        remove.regions <- GenomicRanges::makeGRangesFromDataFrame(remove.regions)
    } else {
        remove.regions <- NULL
    }
    
    hicexp <- multiHiCcompare::make_hicexp(
        data_list = matrices,
        groups = object@conditions,
        remove.regions = remove.regions,
        remove_zeros = FALSE,
        filter = TRUE,
        zero.p = 1,
        A.min = 0
    )
    
    normalized <- multiHiCcompare::cyclic_loess(hicexp, parallel = parallel)
    output <- multiHiCcompare::hic_table(normalized) %>% as_tibble() %>% dplyr::select(-D)
    
    colnames(output) <- c(
        "chromosome", "position.1", "position.2", seq_along(object@replicates)
    )

    object@interactions <- output %>%
        gather(
            as.character(seq_along(object@replicates)),
            key = "i",
            value = "value"
        ) %>%
        dplyr::mutate(i = factor(as.integer(i))) %>%
        dplyr::mutate(condition = factor(object@conditions[i])) %>%
        dplyr::mutate(replicate = factor(object@replicates[i])) %>%
        dplyr::mutate(chromosome = factor(object@chromosomes[chromosome])) %>%
        dplyr::select(chromosome, position.1, position.2, condition, replicate, value)

    return(object)
}
