#' @description
#' Returns the provided limit if valid, or the minimium and maximium of the
#' chromosome's positions.
#'
#' @param xlim
#' A numeric vector of a minimum and a maximum limit for the x axis.
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' The name of a chromosome.
#'
#' @return
#' A numeric vector of a minimum and a maximum limit for the x axis.
#'
#' @keywords internal
#' @noRd
.validateXlim <- function(xlim, object, chromosomeName) {
    if (!is.null(xlim)) {
        if (length(xlim) != 2) {
            message(
                "Expected a vector of two numbers but received '",
                paste(xlim, collapse = " "),
                "'. Setting xlim to NULL."
            )
            xlim <- NULL
        } else {
            xlim <- sort(xlim)
        }
    }
    if (is.null(xlim)) {
        regions <- InteractionSet::regions(object)
        regions <- regions[GenomeInfoDb::seqnames(regions) == chromosomeName]
        xlim <- c(
            min(GenomicRanges::start(regions), na.rm = TRUE),
            max(GenomicRanges::start(regions), na.rm = TRUE)
        )
    }
    return(xlim)
}

#' @description
#' Extracts legends from ggplot2 objects. Based on \code{gtable::gtable_filter}.
#'
#' @param grob
#' A grid graphical object.
#'
#' @return
#' The legend as a grob.
#'
#' @keywords internal
#' @noRd
.extractLegends <- function(grob) {
    matches <- grepl("guide-box", .subset2(grob$layout, "name"), fixed = FALSE)
    grob$layout <- grob$layout[matches, , drop = FALSE]
    grob$grobs <- grob$grobs[matches]
    return(grob)
}


#' @description
#' Complete the levels of replicates to get balanced condition x replicate
#' @param replicates
#' Vector of replicates for one conditino
#' @param expectedLength
#' Expected length of replicates levels
#' @param condition
#' Name of the condition
#'
#' @return
#' A vector with fictif levels if some missing
#'
#' @keywords internal
#' @noRd
.completeLevels <- function(replicates, expectedLength, condition) {
    if (length(replicates) < expectedLength) {
        complete <- paste0("R.", seq(expectedLength))
        complete[seq(length(replicates))] <- replicates
    } else {
        complete <- replicates
    }
    return(complete)
}

#' @title
#' Compute PCA
#'
#' @description
#' Helper function that computes Principal Components of centroids.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosomeName
#' A chromosome name or index in \code{chromosomes(object)}.
#'
#' @return
#' A list, with a \code{data.table}, which contains the PCA,
#' and the variability explained by the 2 first axes.
#'
#' @keywords internal
#' @noRd
.computePca <- function(object, chromosomeName) {
    df <- object@centroids[
        chromosome == chromosomeName,
        .(condition, compartment, centroid)
    ]
    if (nrow(df) == 0) {
        message("No centroids for chromosome ", chromosomeName, ".")
        return(NULL)
    }
    conditions <- df$condition
    compartments <- df$compartment
    
    df <- lapply(df$centroid, unlist)
    df <- do.call("rbind", df)
    
    pca <- stats::prcomp(df)
    varpca <- pca$sdev ^ 2
    propvar <- varpca / sum(varpca)
    
    pca <- as.data.table(pca$x)
    pca[, condition := conditions]
    pca[, compartment := compartments]
    
    return(list(PCA = pca, propvar = propvar))
}
