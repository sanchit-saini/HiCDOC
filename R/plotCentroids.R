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

#' @title
#' Plot centroids.
#'
#' @description
#' Plots the result of the PCA on the compartments' centroids.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosome
#' A chromosome name or index in \code{chromosomes(object)}.
#' @param size
#' Size of each point. Defaults to 2.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' data(exampleHiCDOCDataSetProcessed)
#' plotCentroids(exampleHiCDOCDataSetProcessed, chromosome = 1)
#'
#' @export
plotCentroids <- function(object, chromosome, size = 2) {
    .validateSlots(object, slots = "centroids")
    if (length(chromosome) > 1) {
        warning(
            "`chromosome` should be of length 1, ",
            "taking the first one"
        )
        chromosome < chromosome[1]
    }
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")

    pcaData <- .computePca(object, chromosomeName)
    pca     <- pcaData$PCA
    propvar <- pcaData$propvar
    propvar <- paste(round(100 * propvar, 2), "%")

    plot <- ggplot(
        pca,
        aes(
            x = PC1,
            y = PC2,
            color = compartment,
            shape = condition
        )
    ) + geom_point(size = size) + labs(
        title = paste0("PCA on centroids of chromosome ", chromosomeName),
        x = paste("PC1 ", propvar[1]),
        y = paste("PC2 ", propvar[2])
    )
    return(plot)
}
