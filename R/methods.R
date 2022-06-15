#' @title
#' Methods to access a \code{\link{HiCDOCDataSet}} components.
#'
#' @name
#' HiCDOCDataSet-methods
#'
#' @description
#' Retrieve information and results from a \code{\link{HiCDOCDataSet}}.
#' 
#' @examples
#' # Load an example dataset already processed 
#' # (i.e. after the detection of compartments)
#' data(exampleHiCDOCDataSetProcessed)
#' 
#' exampleHiCDOCDataSetProcessed
#' chromosomes(exampleHiCDOCDataSetProcessed)
#' sampleConditions(exampleHiCDOCDataSetProcessed)
#' sampleReplicates(exampleHiCDOCDataSetProcessed)
#' compartments(exampleHiCDOCDataSetProcessed)
#' differences(exampleHiCDOCDataSetProcessed)
#' concordances(exampleHiCDOCDataSetProcessed)
#' 
#' @return
#' A character vector (for \code{chromosomes}, \code{sampleConditions},
#' \code{sampleReplicates}), 
#' or a GRanges object
#' (for \code{compartments}, \code{concordances}, \code{differences}).
NULL

#### chromosomes ####
#' Retrieves the vector of chromosome names.
#' @rdname HiCDOCDataSet-methods
#' @export
setMethod("chromosomes", "HiCDOCDataSet", function(object) {
    object@chromosomes
})

#### sampleConditions ####
#' Retrieves the vector of condition names, one for each sample.
#' @rdname HiCDOCDataSet-methods
#' @export
setMethod("sampleConditions", "HiCDOCDataSet", function(object) {
    object$condition
})

#### sampleReplicates ####
#' Retrieves the vector of replicate names, one for each sample.
#' @rdname HiCDOCDataSet-methods
#' @export
setMethod("sampleReplicates", "HiCDOCDataSet", function(object) {
    object$replicate
})

#### compartments ####
#' Retrieves a \code{GenomicRange} of the compartment of every position
#' @rdname HiCDOCDataSet-methods
#' @export
setMethod("compartments", "HiCDOCDataSet", function(object, passChecks = TRUE) {
    if(passChecks == TRUE){
        compartments <- object@compartments
        compartments[compartments$assignment.check == TRUE & 
                         compartments$centroid.check == TRUE &
                         compartments$assignment.check == TRUE]
    } else {
        object@compartments
    }
})

#### differences ####
#' Retrieves a \code{GenomicRange} of the significant compartment differences
#' @rdname HiCDOCDataSet-methods
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
#' @export
setMethod("concordances", "HiCDOCDataSet", function(object, passChecks = TRUE) {
    if(passChecks == TRUE){
        concordances <- object@concordances
        passingChecks <- object@checks[centroid.check == TRUE &
                                           PC1.check == TRUE &
                                           assignment.check == TRUE, chromosome]
        concordances[concordances@seqnames %in% passingChecks]
    } else {
        object@concordances
    }
})

#### show ####
#' @rdname HiCDOCDataSet-methods
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
            is.null(sampleReplicates(object)) ||
            length(sampleReplicates(object)) == 0
        ) {
            "  None\n"
        } else {
            paste0(
                "  condition ",
                sampleConditions(object),
                ", replicate ",
                sampleReplicates(object),
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
        "  sampleConditions(object)\n",
        "  sampleReplicates(object)\n",
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
#' @export
setMethod("parameters", "HiCDOCDataSet", function(object) {
    object@parameters
})

#### parameters<- ####
#' Change the parameters of a \code{\link{HiCDOCDataSet}}.
#' @rdname HiCDOCDataSet-parameters
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
