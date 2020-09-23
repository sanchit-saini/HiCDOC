##- interactions -------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Accessors for the 'interactions' slot of an HiCDOCExp object
#'
#' The \code{interactions} slot contains the (transformed) interaction profiles.
#'
#' @docType methods
#' @name interactions
#' @rdname interactions
#' @aliases interactions interactions,HiCDOCExp-method
#' @param object An \code{HiCDOCExp} object.
#' @return A tibble
#' @examples
#' exp <- HiCDOCExample()
#' interactions(exp)
#'
#' @export
setMethod(
  f = "interactions",
  signature = "HiCDOCExp",
  definition = function(object) {
  object@interactions
  }
)


##- differences --------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Extracts differentially interacting regions of an HiCDOCExp object
#'
#' This function extracts the differentially interacting regions from
#' \code{\link{HiCDOCExp}}.
#'
#' @docType methods
#' @name differences
#' @rdname differences
#'
#' @aliases differences differences,HiCDOCExp-method
#'
#' @param object An \code{HiCDOCExp} object.
#' @param pvalue Numeric cutoff value for adjusted p-values. Only regions with
#'         adjusted p-values equal or lower than specified are returned.
#'         Default to 1, all regions are returned.
#'
#' @return A \code{GenomicRanges} object of the selected differentially
#' interacting regions.
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp)
#' differences(exp, pvalue = 0.05)
#'
#' @export
setMethod(
  f = "differences",
  signature = "HiCDOCExp",
  definition = function(object, pvalue=1) {
    if (is.null(object@differences)) {
      message(
        "No 'differences' slot found in the ",
        "HiCDOCExp object. Run HiCDOC first."
      )
    }
    else if (length(object@differences) == 0) {
      message("No 'differences' found.")
      return(GRanges())
    }
    else {
      ##- checking input value ---------------------------------#
      ##--------------------------------------------------------#
      if (length(pvalue) != 1) {
        stop("'pvalue' must be a single value.", call. = FALSE)
      }

      if (is.null(pvalue) || !is.numeric(pvalue) ||
        !is.finite(pvalue)) {
        stop("'pvalue' value must be numeric.", call. = FALSE)
      }

      if ((pvalue > 1) || (pvalue < 0)) {
        stop(
          "'pvalue' value ", pvalue, ", outside the interval [0,1].",
          call. = FALSE
        )
      }
      ##- end check -------------------------------------------#

      gr <- object@differences %>%
        filter(abs(padj) <= pvalue) %>%
        mutate(start = start + 1)
      if (nrow(gr) == 0) {
        message(paste0(
          "No 'differences' found at p-value ",
          pvalue,
          ": best is: ",
          min(abs(object@differences$padj)),
          "."
        ))
        return(GRanges())
      }
      return (makeGRangesFromDataFrame(
        gr,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE
      ))
    }
  }
)


##- concordances -------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Extracts concordances
#'
#' This function extracts the concordances from \code{\link{HiCDOCExp}}.
#'
#' @docType methods
#' @name concordances
#' @rdname concordances
#'
#' @aliases concordances concordances,HiCDOCExp-method
#'
#' @param object An \code{HiCDOCExp} object.
#'
#' @return A \code{tibble} object of the concordance
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp)
#' concordances(exp)
#'
#' @export
setMethod(
  f = "concordances",
  signature = "HiCDOCExp",
  definition = function(object) {
    if (is.null(object@concordances)) {
      message(
        "No 'concordances' slot found in the ",
        "HiCDOCExp object. Run HiCDOC first."
      )
    } else {
      return(object@concordances)
    }
  }
)


##- compartments -------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Extracts compartments
#'
#' This function extracts the compartments from \code{\link{HiCDOCExp}}.
#'
#' @docType methods
#' @name compartments
#' @rdname compartments
#'
#' @aliases compartments compartments,HiCDOCExp-method
#'
#' @param object An \code{HiCDOCExp} object.
#'
#' @return A \code{tibble} object of the compartments
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp)
#' compartments(exp)
#'
#' @export
setMethod(
  f = "compartments",
  signature = "HiCDOCExp",
  definition = function(object) {
    if (is.null(object@compartments)) {
      message(
        "No 'compartments' slot found in the ",
        "HiCDOCExp object. Run HiCDOC first."
      )
    } else {
      grl <- object@compartments %>%
          mutate(start = position + 1) %>%
          mutate(end = start + object@binSize - 1) %>%
          select(-position)
      return(grl)
      # grl <- object@compartments %>%
      #   mutate(start = position + 1) %>%
      #   mutate(end = start + object@binSize - 1) %>%
      #   select(-position) %>%
      #   mutate(condition = factor(condition)) %>%
      #   rename(compartment = value) %>%
      #   GenomicRanges::makeGRangesListFromDataFrame(
      #     keep.extra.columns = TRUE,
      #     ignore.strand = TRUE,
      #     split.field = "condition"
      #   )
      # grl2 <- lapply(grl, function(x) { split(x, ~ compartment) })
      # grl3 <- as(lapply(grl2, function(x) {
      #   unlist(as(lapply(names(x), function(y) {
      #     z <- GenomicRanges::reduce(x[[y]])
      #     z$compartment <- factor(y)
      #     return(z)
      #   }), "GRangesList"))
      # }), "GRangesList")
      # return(grl3)
    }
  }
)


##- parameters ---------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Accessors for the 'parameters' slot of an HiCDOCExp object
#'
#' The \code{parameters} slot holds the parameter values
#' used in an experiment as a named \code{list}. Default values
#' exist for parameters, but these can also be supplied as input
#' values in the \code{useParameters} argument of the \code{\link{HiCDOC}}
#' function or using the assignment function \code{\link{parameters<-}}.
#'
#' Parameters in a HiCDOC experiment.
#'
#' \subsection{Global parameters}{
#'  \describe{
#'    \item{\code{minDepth}}{The cutoff to filter the base-level
#'        coverage. Bases where at least one sample has (normalized)
#'        coverage greater than \code{minDepth} be been retained.
#'        Default to \code{10}.}
#'  }
#' }
#'
#' @docType methods
#' @name parameters
#' @rdname parameters
#' @aliases parameters parameters,HiCDOCExp-method
#' parameters<- parameters<-,HiCDOCExp-method
#' @param object An \code{HiCDOCExp} object.
#' @param value  A named \code{list} containing valid parameters. See details.
#' @return The named list of the parameters used in the analysis.
#' @seealso
#' \code{useParameters} argument in \code{\link{HiCDOC}} function.
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp)
#' print(parameters(exp))
#'
#' parameters(srnaExp) <- list("minSize" = 1, "maxSize" = 1500)
#'
#' @export
setMethod(
  f = "parameters",
  signature = "HiCDOCExp",
  definition = function(object) {
    if (is.null(object@parameters)) {
      message(
        "No 'parameters' slot found in the HiCDOCExp ",
        "object. Run HiCDOC first or assign a named ",
        "list of valid parameters. See help(parameters) ",
        "for details."
      )
    } else {
      object@parameters
      class(object@parameters) <- "HiCDOCParameters"
      return(invisible(object@parameters))
    }
  }
)

#' @name parameters
#' @rdname parameters
#' @exportMethod "parameters<-"
setReplaceMethod(
  "parameters",
  signature(object = "HiCDOCExp", value = "ANY"),
  function(object, value) {

    ##- checking input value ---------------------------------#
    ##--------------------------------------------------------#
    defaultParNames <- names(HiCDOCDefaultParameters)

    if (!is.null(object@parameters)) {
      HiCDOCDefaultParameters <- object@parameters
    }

    if (!is(value, "list")) {
      print(value)
      print(typeof(value))
      print(class(value))
      stop(
        "'value' must be a named list. See ",
        "help(parameters) for details.",
        call. = FALSE
      )
    }

    valueNames <- names(value)

    if (any(duplicated(valueNames))) {
      stop(
        "duplicate name parameters in 'value'. See ",
        "help(parameters) for details.",
        call. = FALSE
      )
    }

    if (!all(valueNames %in% defaultParNames)) {
      stop(
        "'value' must be a named list of valid ",
        "parameters. See help(parameters) for details.",
        call. = FALSE
      )
    }

    ##- individual parameters
    HiCDOCDefaultParameters[valueNames] <- value
    checkParameters(HiCDOCDefaultParameters)

    ##- end check -------------------------------------------#

    object@parameters <- HiCDOCDefaultParameters
    object
  }
)

##- show ---------------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname HiCDOCExp
#' @param object An \code{HiCDOCExp} object.
#' @return The \code{show} method informatively display object contents.
#' @export
setMethod(
  f = "show",
  signature = "HiCDOCExp",
  definition = function(object) {
    cat("Object of class HiCDOCExp.\n", "Sample information\n")
    print("TODO")
  }
)


##- print method for parameters ----------------------------------------------#
##----------------------------------------------------------------------------#
#' Dispatch print method for the parameters used by an \code{HiCDOC} object.
#'
#' @docType methods
#' @name parameters
#' @rdname parameters
#' @aliases parameters parameters,HiCDOCExp-method
#' @param x The first element of the parameters used by an \code{HiCDOC} object
#' @param ... The other elements of the parameters
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(srnaExp)
#' print(parameters(exp))
#'
#' @export
printHiCDOCParameters <- function(x, ...) {

  cat("\n Global parameters: \n", "------------------ \n")
  df <- data.frame(value = unlist(x[seq_len(6)]))
  print(df)

  cat("\n Constrained K-means parameters: \n", "---------------------- \n")
  df <- data.frame(value = unlist(x[7:10]))
  print(df)
}
