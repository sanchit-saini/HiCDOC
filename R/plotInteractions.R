#' Interactions plot using facet_grid
#'
#' @param interactions
#' a tibble
#' @param xylim
#' a length 2 vector : limits of the matrix
#' @param transform
#' character : transformation to be passed to \code{scale_fill_gradientn}
#' @param colours
#' vector of colors to be passed to \code{scale_fill_gradientn}
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
    chromosomeName
) {
    plot <-
        ggplot(
            data = dataplot,
            aes(x = start1, y = start2, z = interaction)
        ) +
        geom_tile(aes(fill = interaction), na.rm = TRUE) +
        geom_tile(
            data = dataplot[start1 != start2, ],
            aes(x = start2, y = start1, fill = interaction),
            na.rm = TRUE
        ) +
        coord_fixed(ratio = 1) +
        theme_bw() +
        xlim(xylim) +
        scale_y_reverse(limits = rev(xylim)) +
        facet_grid(
            condition ~ replicate,
            drop = FALSE
        ) +
        labs(title = paste0("Chromosome ", chromosomeName), x = "", y = "") +
        scale_fill_gradientn(
            colours = colours,
            trans = transform,
            name = "Intensity",
            na.value = "transparent"
        )
    return(plot)
}


#' Interactions plot using facet_wrap
#'
#' @param interactions
#' a tibble
#' @param xylim
#' a length 2 vector : limits of the matrix
#' @param transform
#' character : transformation to be passed to \code{scale_fill_gradientn}
#' @param colours
#' vector of colors to be passed to \code{scale_fill_gradientn}
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
    chromosomeName,
    totalRows,
    totalCols
) {
    plot <-
        ggplot(
            data = dataplot,
            aes(x = start1, y = start2, z = interaction)
        ) +
        geom_tile(aes(fill = interaction), na.rm = TRUE) +
        geom_tile(
            data = dataplot[start1 != start2, ],
            aes(x = start2, y = start1, fill = interaction),
            na.rm = TRUE
        ) +
        coord_fixed(ratio = 1) +
        theme_bw() +
        xlim(xylim) +
        scale_y_reverse(limits = rev(xylim)) +
        facet_wrap(
            . ~ variable,
            ncol = totalCols,
            nrow = totalRows,
            drop = FALSE
        ) +
        labs(title = paste("Chromosome ", chromosomeName), x = "", y = "") +
        scale_fill_gradientn(
            colours = colours,
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
#' \code{\link[ggplot2]{scale_fill_gradientn}} for other accepted values.
#' @param colours
#' A character vector of colours to use for the gradient. See
#' \code{\link[ggplot2]{scale_fill_gradientn}} for more info. Defaults to
#' \code{c("#F6FFB8", "#FF00CC", "#310038")}.
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
    colours = c("#F6FFB8", "#FF00CC", "#310038")
) {

    chromosomeName <- .validateNames(object, chromosome, "chromosomes")
    if (is.null(transform)) transform <- "identity"
    
    rowsChromosome <- (S4Vectors::mcols(object)$Chr == chromosomeName)
    assayChromosome <- SummarizedExperiment::assay(object[rowsChromosome,])
    assayChromosome <- as.data.table(assayChromosome)
    setnames(assayChromosome, paste(object$condition, object$replicate, sep="_"))
    
    interactionsChromosome <- InteractionSet::interactions(object[rowsChromosome,]) 
    interactionsChromosome <- as.data.table(interactionsChromosome)
    dataplot <- cbind(interactionsChromosome[,.(seqnames = seqnames1, start1, start2)], 
                      assayChromosome) 
    dataplot <- data.table::melt.data.table(dataplot, 
                    id.vars = c("seqnames", "start1", "start2"),
                    value.name = "interaction",
                    variable.factor = FALSE)
    dataplot <- dataplot[!is.na(interaction)]
    dataplot[, c("condition", "replicate") := 
                 data.table::tstrsplit(variable, "_", fixed=TRUE)]
    dataplot[,condition := factor(condition, levels=sort(unique(object$condition)))]
    dataplot[,replicate := factor(replicate, levels=sort(unique(object$replicate)))]

    if (nrow(dataplot) == 0) {
        message("No interactions for chromosome ", chromosomeName, ".")
        return(NULL)
    }
    
    regionsChromosome <- as.data.table(InteractionSet::regions(object))
    regionsChromosome <- regionsChromosome[seqnames == chromosomeName]
    xylim <- c(min(regionsChromosome$start), max(regionsChromosome$start))
    
    if (length(unique(object$replicate)) <= max(table(object$condition))) {
        plot <- .plotInteractionsGrid(dataplot, xylim, transform, colours,
                    chromosomeName)
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
            SIMPLIFY = FALSE)
        allLevels <- unlist(existing, use.names = FALSE)
        dataplot[,variable := factor(variable, levels = allLevels)]

        plot <- .plotInteractionsWrap(dataplot, xylim, transform, colours,
                    chromosomeName, totalRows, totalCols)
    }
    return(plot)
}
