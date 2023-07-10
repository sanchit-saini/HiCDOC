#' Interactions plot using facet_grid
#'
#' @param interactions
#' a data.table
#' @param xylim
#' a length 2 vector : limits of the matrix
#' @param transform
#' character : transformation to be passed to \code{scale_fill_gradient2}
#' @param colours
#' vector of colors of length 3 to be passed to \code{scale_fill_gradient2}
#' @param midpoint
#' midpoint value to be passed to \code{scale_fill_gradient2}
#' @param chromosomeName
#' Character, name of the chromosome
#' @return
#' a ggplot object
#'
#' @keywords internal
#' @noRd
.plotInteractionsGrid <- function(
    dataplot,
    xylim,
    transform,
    colours,
    midpoint,
    chromosomeName
) {
    plot <- ggplot(
        data = dataplot,
        aes(x = start1, y = start2, z = interaction)
    ) + geom_tile(
        aes(fill = interaction),
        na.rm = TRUE
    ) + geom_tile(
        data = dataplot[start1 != start2, ],
        aes(x = start2, y = start1, fill = interaction),
        na.rm = TRUE
    ) + coord_fixed(ratio = 1) + theme_bw() + xlim(xylim) + scale_y_reverse(
        limits = rev(xylim)
    ) + facet_grid(
        condition ~ replicate,
        drop = FALSE
    ) + labs(
        title = paste0("Chromosome ", chromosomeName),
        x = "",
        y = ""
    ) + scale_fill_gradient2(
        low=colours[1],
        mid=colours[2],
        high=colours[3],
        midpoint = midpoint,
        trans = transform,
        name = "Intensity",
        na.value = "transparent"
    )
    return(plot)
}


#' Interactions plot using facet_wrap
#'
#' @param interactions
#' a data.table
#' @param xylim
#' a length 2 vector : limits of the matrix
#' @param transform
#' character : transformation to be passed to \code{scale_fill_gradient2}
#' @param colours
#' vector of colors of length 3 to be passed to \code{scale_fill_gradient2}
#' @param midpoint
#' midpoint value to be passed to \code{scale_fill_gradient2}
#' @param chromosomeName
#' Character, name of the chromosome
#' @param totalRows
#' Integer, number of rows in facet_wrap
#' @param totalCols
#' Integer, number of colums in facet_wrap
#'
#' @return
#' a ggplot object
#'
#' @keywords internal
#' @noRd
.plotInteractionsWrap <- function(
    dataplot,
    xylim,
    transform,
    colours,
    midpoint,
    chromosomeName,
    totalRows,
    totalCols
) {
    plot <- ggplot(
        data = dataplot,
        aes(x = start1, y = start2, z = interaction)
    ) + geom_tile(
        aes(fill = interaction),
        na.rm = TRUE
    ) + geom_tile(
        data = dataplot[start1 != start2, ],
        aes(x = start2, y = start1, fill = interaction),
        na.rm = TRUE
    ) + coord_fixed(ratio = 1) + theme_bw() + xlim(xylim) + scale_y_reverse(
        limits = rev(xylim)
    ) + facet_wrap(
        . ~ variable,
        ncol = totalCols,
        nrow = totalRows,
        drop = FALSE
    ) + labs(
        title = paste("Chromosome ", chromosomeName),
        x = "",
        y = ""
    ) + scale_fill_gradient2(
        low=colours[1],
        mid=colours[2],
        high=colours[3],
        midpoint = midpoint,
        trans = transform,
        name = "Intensity",
        na.value = "transparent"
    )
    return(plot)
}

#' @title
#' Plot interaction matrices.
#'
#' @description
#' Plots the interaction matrices as heatmaps.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param chromosome
#' A chromosome name or index in \code{chromosomes(object)}.
#' @param transform
#' Transformation of the color scale. Default to NULL (no transformation). See
#' \code{\link[ggplot2]{scale_fill_gradient2}} for other accepted values.
#' @param colours
#' A character vector colours of length 3 to use for the gradient. See
#' \code{\link[ggplot2]{scale_fill_gradient2}} for more info. Defaults to
#' \code{c("low"=#2c7bb6", "mid"=#ffffbf", "high"="#d7191c")}.
#' @param midpoint
#' midpoint value to be passed to \code{scale_fill_gradient2}.
#' Default to 0.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' data(exampleHiCDOCDataSet)
#' plotInteractions(exampleHiCDOCDataSet, chromosome = 1)
#'
#' @export
plotInteractions <- function(
    object,
    chromosome,
    transform = NULL,
    colours = c("low" = "#2c7bb6", "mid" = "#ffffbf", "high" = "#d7191c"),
    midpoint = 0
) {
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")
    if (is.null(transform)) {
        transform <- "identity"
    }

    rowsChromosome <- (S4Vectors::mcols(object)$chromosome == chromosomeName)
    assayChromosome <- SummarizedExperiment::assay(object[rowsChromosome, ])
    assayChromosome <- data.table::as.data.table(assayChromosome)
    data.table::setnames(
        assayChromosome,
        paste(object$condition, object$replicate, sep = "_")
    )

    interactionsChromosome <- InteractionSet::interactions(
        object[rowsChromosome, ]
    )
    interactionsChromosome <- as.data.table(interactionsChromosome)
    dataplot <- base::cbind(
        interactionsChromosome[, .(
            seqnames = seqnames1,
            start1,
            start2
        )],
        assayChromosome
    )
    dataplot <- data.table::melt.data.table(
        dataplot,
        id.vars = c("seqnames", "start1", "start2"),
        value.name = "interaction",
        variable.factor = FALSE
    )
    dataplot <- dataplot[!is.na(interaction)]
    if (nrow(dataplot) == 0) {
        message("No interactions for chromosome ", chromosomeName, ".")
        return(NULL)
    }
    dataplot[, c("condition", "replicate") := data.table::tstrsplit(
        variable,
        "_",
        fixed = TRUE
    )]
    dataplot[, condition := factor(
        condition,
        levels = sort(unique(object$condition))
    )]
    dataplot[, replicate := factor(
        replicate,
        levels = sort(unique(object$replicate))
    )]
    
    regionsChromosome <- as.data.table(InteractionSet::regions(object))
    regionsChromosome <- regionsChromosome[seqnames == chromosomeName]
    xylim <- c(min(regionsChromosome$start), max(regionsChromosome$start))
    
    if(length(colours) != 3){
        stop("`colours` must be a vector of length 3.")
    } else {
        if(is.null(names(colours))) names(colours) <- c("low", "mid", "high")
        if(!identical(sort(names(colours)), c("high", "low", "mid"))){
            stop("`colours` are supposed to be named 'low', 'mid' and 'high")
        }
    }
    
    if (length(unique(object$replicate)) <= max(table(object$condition))) {
        plot <- .plotInteractionsGrid(
            dataplot,
            xylim,
            transform,
            colours,
            midpoint,
            chromosomeName
        )
    } else {
        totalLevels <- table(object$condition)
        totalCols <- max(totalLevels)
        totalRows <- length(unique(object$condition))

        existing <- by(object$replicate, object$condition, unique)
        existing <- lapply(existing, as.character)
        existing <- lapply(existing, .completeLevels, totalCols)
        existing <- mapply(
            paste,
            names(existing),
            existing,
            sep = "_",
            SIMPLIFY = FALSE
        )
        allLevels <- unlist(existing, use.names = FALSE)
        dataplot[, variable := factor(variable, levels = allLevels)]

        plot <- .plotInteractionsWrap(
            dataplot,
            xylim,
            transform,
            colours,
            midpoint,
            chromosomeName,
            totalRows,
            totalCols
        )
    }
    return(plot)
}
