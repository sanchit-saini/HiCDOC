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
    interactions, 
    xylim, 
    transform, 
    colours, 
    chromosomeName
){
    plot <- 
        ggplot(
            data = interactions,
            aes(x = position.1, y = position.2, z = interaction)
        ) +
        geom_raster(aes(fill = interaction), na.rm = TRUE) +
        geom_raster(
            data = interactions[interactions$bin.1 != interactions$bin.2, ],
            aes(x = position.2, y = position.1, fill = interaction),
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
        labs(title = paste("Chromosome ", chromosomeName), x = "", y = "") +
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
    interactions, 
    xylim, 
    transform, 
    colours, 
    chromosomeName,
    totalRows,
    totalCols
){
    plot <- 
        ggplot(
            data = interactions,
            aes(x = position.1, y = position.2, z = interaction)
        ) +
        geom_raster(aes(fill = interaction), na.rm = TRUE) +
        geom_raster(
            data = interactions[interactions$bin.1 != interactions$bin.2, ],
            aes(x = position.2, y = position.1, fill = interaction),
            na.rm = TRUE
        ) +
        coord_fixed(ratio = 1) +
        theme_bw() +
        xlim(xylim) +
        scale_y_reverse(limits = rev(xylim)) +
        facet_wrap(
            . ~ facet,
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

#' Complete the levels of replicates to get balanced condition x replicate
#'
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
.completeLevels <- function(replicates, expectedLength, condition){
    if(length(replicates) < expectedLength){
        complete <- paste0("R.", seq_len(expectedLength))
        complete[1:length(replicates)] <- replicates
    } else {
        complete <- replicates
    }
    return(complete)
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
#' Transformation of the color scale. Set to NULL for no transformation. See
#' \code{\link[ggplot2]{scale_fill_gradientn}} for other accepted values.
#' Defaults to "log2".
#' @param colours
#' A character vector of colours to use for the gradient. See
#' \code{\link[ggplot2]{scale_fill_gradientn}} for more info. Defaults to
#' \code{c("#000066", "#ffffbf", "#990000")}.
#'
#' @return
#' A \code{ggplot}.
#'
#' @examples
#' object <- HiCDOCDataSetExample()
#' p <- plotInteractions(object, chromosome = 1)
#'
#' @export
plotInteractions <- function(
    object,
    chromosome,
    transform = "log2",
    colours = c("#000066", "#ffffbf", "#990000")
) {

    .validateSlots(
        object,
        slots = c(
            "interactions",
            "conditions",
            "totalBins",
            "binSize",
            "positions"
        )
    )
    chromosomeName <- .validateNames(object, chromosome, "chromosomes")
    if (is.null(transform)) transform <- "Identity"

    positions <-
        object@positions %>%
        dplyr::filter(chromosome == chromosomeName)
    interactions <-
        object@interactions %>%
        dplyr::filter(chromosome == chromosomeName) %>%
        dplyr::filter(interaction > 0) %>%
        dplyr::left_join(
            positions %>% dplyr::select(
                bin.1 = bin,
                position.1 = start
            ),
            by = "bin.1"
        ) %>%
        dplyr::left_join(
            positions %>% dplyr::select(
                bin.2 = bin,
                position.2 = start
            ),
            by = "bin.2"
        )

    if (nrow(interactions) == 0) {
        message("No interaction data to plot.")
        return(NULL)
    }
    
    xylim <- c(min(positions$start), max(positions$start))
    
    if(prod(dim(table(interactions$condition, interactions$replicate))) <=
       sum(table(object@conditions))){
        plot <- .plotInteractionsGrid(interactions, xylim, transform, colours,
                                      chromosomeName)
    } else {
        totalLevels <- table(object@conditions)
        totalCols <- max(totalLevels)
        totalRows <- length(unique(object@conditions))
        
        existing <- by(interactions$replicate, interactions$condition, unique)
        existing <- lapply(existing, as.character)
        existing <- lapply(existing, .completeLevels, totalCols)
        existing <- mapply(paste, unique(object@conditions), existing, sep="_", SIMPLIFY = FALSE)
        allLevels <- unlist(existing, use.names = FALSE)
        
        # Construct facet variable
        interactions %<>% 
            dplyr::mutate(facet = paste(condition ,replicate, sep="_")) %>%
            dplyr::mutate(facet = factor(facet, levels = allLevels))
        
        plot <- .plotInteractionsWrap(interactions, xylim, transform, colours,
                                      chromosomeName, totalRows, totalCols)
    }
    return(plot)
}
