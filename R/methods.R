#### chromosomes ####
#' Retrieves the vector of chromosome names.
#' @rdname HiCDOCDataSet-methods
#' @usage
#' NULL
#' @export
setMethod("chromosomes", "HiCDOCDataSet", function(object) {
    object@chromosomes
})

#### conditions ####
#' Retrieves the vector of condition names.
#' @rdname HiCDOCDataSet-methods
#' @usage
#' NULL
#' @export
setMethod("conditions", "HiCDOCDataSet", function(object) {
    object$condition
})

#### replicates ####
#' Retrieves the vector of replicate names.
#' @rdname HiCDOCDataSet-methods
#' @usage
#' NULL
#' @export
setMethod("replicates", "HiCDOCDataSet", function(object) {
    object$replicate
})

# TODO : j'aimerai bien doubler les fonctions assay et regions.
# #### assay ####
# #' Retrieves the assay matrix of interactions
# #' @rdname HiCDOCDataSet-methods
# #' @usage
# #' NULL
# #' @export
# setMethod("assay",
#           signature(object="HiCDOCDataSet"),
#           function(object) {
#               SummarizedExperiment::assay(object@interactions)
#           })
#
# #### regions ####
# #' Retrieves the regions of interactions
# #' @rdname HiCDOCDataSet-methods
# #' @usage
# #' NULL
# #' @export
# setMethod("regions",
#           signature(object="HiCDOCDataSet"),
#           function(object) InteractionSet::regions(object@interactions))
#
# #### interactions ####
# #' Retrieves a data.table of the interactions.
# #' @rdname HiCDOCDataSet-methods
# #' @usage
# #' NULL
# #' @export
# setMethod("interactions",
#           signature(object="HiCDOCDataSet"),
#           function(object) {
#     if (is.null(object@interactions)) return(NULL)
#     interactions <- InteractionSet::interactions(object@interactions)
#     return(interactions)
# })

#### compartments ####
#' Retrieves a \code{GenomicRange} of the compartment of every position
#' @rdname HiCDOCDataSet-methods
#' @usage
#' NULL
#' @export
setMethod("compartments", "HiCDOCDataSet", function(object) {
    object@compartments
})

#### differences ####
#' Retrieves a \code{GenomicRange} of the significant compartment differences
#' @rdname HiCDOCDataSet-methods
#' @usage
#' NULL
#' @export
setMethod("differences", "HiCDOCDataSet", function(object, threshold = NULL) {
    if (is.null(object@differences)) {
        return(NULL)
    }

    if (
        !is.null(threshold) && (
            !is.numeric(threshold) ||
            length(threshold) > 1
        )
    ) {
        stop("'threshold' must be a number.", call. = FALSE)
    }
    differences <- object@differences
    if (!is.null(threshold)) {
        differences <- differences[differences$pvalue.adjusted <= threshold]
    }

    if (length(differences) == 0) {
        if (is.null(threshold)) {
            message("No differences found.")
        } else {
            message(
                "No differences found with adjusted p-value <= ",
                threshold,
                "."
            )
        }
        return(NULL)
    }

    return(differences)
})

#### concordances ####
#' Retrieves a \code{GenomicRange} of the concordance (confidence in assigned
#' compartment) of every position in every replicate.
#' @rdname HiCDOCDataSet-methods
#' @usage
#' NULL
#' @export
setMethod("concordances", "HiCDOCDataSet", function(object) {
    object@concordances
})

#### show ####
#' @describeIn HiCDOCDataSet-methods
#' Describes the object and its methods.
#' @usage
#' NULL
#' @export
setMethod("show", "HiCDOCDataSet", function(object) {
    cat(
        "Object of class 'HiCDOCDataSet'\n\n",
        "- Inputs:\n",
        paste0(
            "  ",
            vapply(
                object@input,
                function(x) paste0(x),
                character(ifelse(
                    length(object@input) > 0,
                    length(object@input[[1]]),
                    1
                ))
            ),
            "\n"
        ),
        "\n",
        "- Chromosomes:\n  ",
        if (
            is.null(object@chromosomes) ||
            length(object@chromosomes) == 0
        ) {
            "None"
        } else {
            paste(object@chromosomes, collapse = ", ")
        },
        "\n\n",
        "- Replicates:\n",
        if (
            is.null(replicates(object)) ||
            length(replicates(object)) == 0
        ) {
            "  None\n"
        } else {
            paste0(
                "  condition ",
                conditions(object),
                ", replicate ",
                replicates(object),
                "\n"
            )
        },
        "\n",
        "- Parameters:\n",
        paste0(
            "  ",
            vapply(
                seq_along(parameters(object)),
                function(x) {
                    paste(
                        names(parameters(object))[x],
                        '=',
                        parameters(object)[x]
                    )
                },
                character(1)
            ),
            "\n"
        ),
        "\n",
        "- Methods:\n",
        "  chromosomes(object)\n",
        "  conditions(object)\n",
        "  replicates(object)\n",
        "  compartments(object)\n",
        "  differences(object)\n",
        "  concordances(object)\n",
        "  parameters(object)\n",
        "  parameters(object) <- list()\n\n",
        sep = ""
    )
})

#### parameters ####
#' Access the parameters of a \code{\link{HiCDOCDataSet}}.
#' @rdname HiCDOCDataSet-parameters
#' @usage
#' NULL
#' @export
setMethod("parameters", "HiCDOCDataSet", function(object) {
    object@parameters
})

#### parameters<- ####
#' Change the parameters of a \code{\link{HiCDOCDataSet}}.
#' @rdname HiCDOCDataSet-parameters
#' @usage
#' NULL
#' @export
setReplaceMethod("parameters", "HiCDOCDataSet", function(object, value) {
    defaultParameterNames <- names(defaultHiCDOCParameters)

    if (!is(value, "list")) {
        stop(
            "'parameters' must be a named list.\n",
            "No parameters were updated. ",
            "See 'help(parameters)' for details.",
            call. = FALSE
        )
    }
    parameterNames <- names(value)

    duplicatedParameterNames <- unique(
        parameterNames[duplicated(parameterNames)]
    )

    if (length(duplicatedParameterNames) > 0) {
        stop(
            "Duplicate parameter",
            if (length(duplicatedParameterNames) != 1) "s",
            " provided: ",
            paste(duplicatedParameterNames, collapse = ", "),
            "\nNo parameters were updated. ",
            "See 'help(parameters)' for details.",
            call. = FALSE
        )
    }

    invalidParameterNames <- parameterNames[
        !(parameterNames %in% defaultParameterNames)
    ]

    if (length(invalidParameterNames) > 0) {
        stop(
            "Invalid parameter",
            if (length(invalidParameterNames) != 1) "s",
            " provided: ",
            paste(invalidParameterNames, collapse = ", "),
            "\nNo parameters were updated. ",
            "See 'help(parameters)' for details.",
            call. = FALSE
        )
    }

    object@parameters[parameterNames] <- value

    return(object)
})
