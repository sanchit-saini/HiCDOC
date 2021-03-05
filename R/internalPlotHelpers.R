#' Return valid xlim
#'
#' test if the xlim is valid (length(2)) or \code{NULL} and return
#' the xlim for the plot.
#' If \code{NULL} or invalid return the min and max of \code{positions}.
#'
#' @param xlim The entred xlim
#' @param positions The positions from the data
#'
#' @return A length 2 numerical vector
#' @keywords internal
#' @noRd
.validateXlim <- function(xlim, object, chromosomeName) {
    if (!is.null(xlim)) {
        if (length(xlim) != 2) {
            message(
                "Incorrect values for xlim (numerical of length 2 expected). ",
                "Setting xlim to NULL."
            )
            xlim <- NULL
        } else {
            xlim <- sort(xlim)
        }
    }
    if (is.null(xlim)) {
        positions <-
            object@positions %>%
            dplyr::filter(chromosome == chromosomeName) %>%
            dplyr::select(start) %>%
            dplyr::pull()
        xlim <- c(min(positions, na.rm = TRUE), max(positions, na.rm = TRUE))
    }
    return(xlim)
}


#' Function to extract legends from ggplot2 objects
#'
#' Inspired from gtable::gtable_filter (not used to avoid dependency)
#'
#' @param x a grob
#'
#' @return the legend as a grob
#' @keywords internal
#' @noRd
.extractLegends <- function(x) {
    matches <- grepl("guide-box", .subset2(x$layout, "name"), fixed = FALSE)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    return(x)
}
