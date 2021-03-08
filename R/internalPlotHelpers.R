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
                paste0(c(
                    "Expected a vector of two numbers but received '",
                    paste(xlim, collapse = " "),
                    "'. Setting xlim to NULL."
                ))
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
