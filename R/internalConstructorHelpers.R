#' #' @description
#' #' Determines the number of bases per bin by finding the minimum distance
#' #' between \code{position.1} and \code{position.2} that is not zero.
#' #'
#' #' @param object
#' #' A \code{\link{HiCDOCDataSet}}.
#' #'
#' #' @return
#' #' The number of bases per bin.
#' #'
#' #' @keywords internal
#' #' @noRd
#' .determineBinSize <- function(object) {
#'     object@interactions %>%
#'         dplyr::filter(position.1 != position.2) %>%
#'         dplyr::mutate(size = abs(position.1 - position.2)) %>%
#'         dplyr::pull(size) %>%
#'         min()
#' }

#' #' @description
#' #' Determines the position corresponding to each bin of each chromosome.
#' #'
#' #' @param object
#' #' A \code{\link{HiCDOCDataSet}}.
#' #'
#' #' @return
#' #' A tibble of 4 columns: chromosome, bin, start, end.
#' #'
#' #' @keywords internal
#' #' @noRd
#' .determinePositions <- function(object) {
#' 
#'     maxPositions <-
#'         object@interactions %>%
#'         dplyr::group_by(chromosome) %>%
#'         dplyr::summarize(
#'             maxPosition.1 = max(position.1),
#'             maxPosition.2 = max(position.2)
#'         ) %>%
#'         dplyr::ungroup() %>%
#'         dplyr::rowwise() %>%
#'         dplyr::mutate(
#'             maxPosition = max(maxPosition.1, maxPosition.2),
#'             .keep = "unused"
#'         )
#' 
#'     positions <-
#'         lapply(
#'             object@chromosomes,
#'             function(chromosomeName) {
#'                 return(data.frame(
#'                     "chromosome" = chromosomeName,
#'                     "start" = seq(
#'                         0,
#'                         maxPositions[
#'                             maxPositions$chromosome == chromosomeName,
#'                         ]$maxPosition,
#'                         by = object@binSize
#'                     )
#'                 ))
#'             }
#'         ) %>%
#'         dplyr::bind_rows() %>%
#'         dplyr::mutate(
#'             chromosome = factor(chromosome, levels = object@chromosomes)
#'         ) %>%
#'         dplyr::group_by(chromosome) %>%
#'         dplyr::arrange(chromosome, start, .by_group = TRUE) %>%
#'         dplyr::mutate(end = dplyr::lead(start) - 1) %>%
#'         dplyr::mutate(bin = dplyr::row_number()) %>%
#'         dplyr::ungroup()
#' 
#'     positions %<>%
#'         dplyr::mutate(
#'             end = ifelse(is.na(end), start + object@binSize - 1, end)
#'         ) %>%
#'         dplyr::select(chromosome, bin, start, end)
#' 
#'     return(positions)
#' }

#' #' @description
#' #' Replaces positions by their corresponding bins.
#' #'
#' #' @param object
#' #' A \code{\link{HiCDOCDataSet}}.
#' #'
#' #' @return
#' #' A tibble of interactions with bin columns replacing position columns.
#' #'
#' #' @keywords internal
#' #' @noRd
#' .replacePositionsByBins <- function(object) {
#' 
#'     interactions <-
#'         object@interactions %>%
#'         dplyr::left_join(
#'             object@positions %>% dplyr::select(
#'                 chromosome,
#'                 position.1 = start,
#'                 bin.1 = bin
#'             ),
#'             by = c("chromosome", "position.1")
#'         ) %>%
#'         dplyr::left_join(
#'             object@positions %>% dplyr::select(
#'                 chromosome,
#'                 position.2 = start,
#'                 bin.2 = bin
#'             ),
#'             by = c("chromosome", "position.2")
#'         ) %>%
#'         dplyr::select(
#'             chromosome,
#'             condition,
#'             replicate,
#'             bin.1,
#'             bin.2,
#'             interaction
#'         )
#' 
#'     return(interactions)
#' }

#' @description
#' Determines the number of bins per chromosome.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A list of the number of bins per chromosome.
#'
#' @keywords internal
#' @noRd
.determineChromosomeSizes <- function(object) {
    totalBins <- S4Vectors::runLength(
        GenomeInfoDb::seqnames(InteractionSet::regions(object@interactions)))
    names(totalBins) <- object@chromosomes
    return(totalBins)
}

#' #' @description
#' #' Reformats the interactions to fit them into a numeric upper triangular
#' #' sparse matrix.
#' #'
#' #' @param interactions
#' #' A tibble of interactions.
#' #'
#' #' @return
#' #' A reformated tibble of interactions.
#' #'
#' #' @keywords internal
#' #' @noRd
#' .reformatInteractions <- function(interactions) {
#' 
#'     interactions %>%
#'     dplyr::mutate(minBin = pmin(bin.1, bin.2)) %>%
#'     dplyr::mutate(maxBin = pmax(bin.1, bin.2)) %>%
#'     dplyr::select(chromosome, condition, replicate, interaction, minBin, maxBin) %>%
#'     dplyr::rename(bin.1 = minBin) %>%
#'     dplyr::rename(bin.2 = maxBin) %>%
#'     dplyr::mutate(interaction = as.numeric(interaction)) %>%
#'     dplyr::filter(interaction > 0) %>%
#'     dplyr::select(chromosome, bin.1, bin.2, condition, replicate, interaction)
#' }

#' #' @description
#' #' Sorts interactions and factorizes chromosomes, conditions, and replicates.
#' #'
#' #' @param interactions
#' #' A tibble of interactions.
#' #' @param chromosomeNames
#' #' The names of chromosomes to sort with.
#' #' @param conditionNames
#' #' The names of conditions to sort with.
#' #' @param replicateNames
#' #' The names of replicates to sort with.
#' #'
#' #' @return
#' #' A sorted tibble of interactions.
#' #'
#' #' @keywords internal
#' #' @noRd
#' .sortInteractions <- function(
#'     interactions,
#'     chromosomeNames,
#'     conditionNames,
#'     replicateNames
#' ) {
#' 
#'     interactions %<>%
#'         dplyr::mutate(
#'             chromosome = factor(
#'                 chromosome,
#'                 levels = gtools::mixedsort(unique(chromosomeNames))
#'             ),
#'             condition = factor(
#'                 condition,
#'                 levels = gtools::mixedsort(unique(conditionNames))
#'             ),
#'             replicate = factor(
#'                 replicate,
#'                 levels = gtools::mixedsort(unique(replicateNames))
#'             ),
#'             interaction = as.numeric(interaction)
#'         ) %>%
#'         dplyr::filter(interaction > 0)
#' 
#'     if ("bin.1" %in% colnames(interactions)) {
#'         interactions %<>%
#'             dplyr::arrange(
#'                 chromosome,
#'                 condition,
#'                 replicate,
#'                 bin.1,
#'                 bin.2
#'             )
#'     } else if ("position.1" %in% colnames(interactions)) {
#'         interactions %<>%
#'             dplyr::arrange(
#'                 chromosome,
#'                 condition,
#'                 replicate,
#'                 position.1,
#'                 position.2
#'             )
#'     } else {
#'         interactions %<>%
#'             dplyr::arrange(
#'                 chromosome,
#'                 condition,
#'                 replicate
#'             )
#'     }
#' 
#'     return(interactions)
#' }

#' @description
#' Fills parameters and slots describing the data. Called by a
#' \code{\link{HiCDOCDataSet}} constructor.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}} ready for analysis.
#'
#' @keywords internal
#' @noRd
.fillHiCDOCDataSet <- function(object) {

    # Fill all other slots than interactions
    .validateSlots(object, "interactions")

    # Chromosomes and their size (max bin)
    object@chromosomes <- GenomeInfoDb::seqlevels(object@interactions)
    object@totalBins <- .determineChromosomeSizes(object)
    object@parameters <- defaultHiCDOCParameters
    
    # Valid conditions and replicats by chromosome (==not empty)
    # maybe do a function for valid conditions and replicats ?
    valids <- S4Vectors::split(SummarizedExperiment::assay(object@interactions), 
                          S4Vectors::mcols(object@interactions)$Chr, drop=FALSE)
    valids <- lapply(valids, colSums, na.rm=TRUE)
    valids <- lapply(valids, function(x) (x>0 & !is.na(x)))

    object@validConditions <-
        lapply(valids, function(x) object@interactions$condition[x])
    object@validReplicates <-
        lapply(valids, function(x) object@interactions$replicat[x])

    # Weakbins
    object@weakBins <- vector("list", length(object@chromosomes))
    names(object@weakBins) <- object@chromosomes

    return(object)
}
