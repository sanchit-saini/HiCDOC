#' @description
#' Returns parameters, updated with their default values if invalid.
#'
#' @param parameters
#' A list of parameters.
#'
#' @return
#' A list of valid parameters.
#'
#' @keywords internal
#' @noRd
.validateParameters <- function(parameters) {

    defaultParameterNames <- names(defaultHiCDOCParameters)
    inputParameterNames <- names(parameters)

    numericParameters <-
        vapply(
            parameters,
            function(parameter) {
                is.numeric(parameter) && length(parameter) == 1
            },
            FUN.VALUE = TRUE
        )

    if (!all(numericParameters)) {
        notNumericParameters <- inputParameterNames[!numericParameters]
        warning(
            "Non-numeric parameter",
            if (length(notNumericParameters) != 1) "s were" else "was",
            " replaced with ",
            if (length(notNumericParameters) != 1) "their" else "its",
            " default",
            if (length(notNumericParameters) != 1) "s",
            ":\n",
            paste0(
                notNumericParameters,
                ":",
                " ",
                parameters[notNumericParameters],
                " -> ",
                defaultHiCDOCParameters[notNumericParameters],
                collapse = "\n"
            ),
            call. = FALSE
        )
        parameters[notNumericParameters] <-
            defaultHiCDOCParameters[notNumericParameters]
    }

    known <- inputParameterNames %in% defaultParameterNames

    if (!all(known)) {
        unknownParameterNames <- inputParameterNames[!known]
        warning(
            "Unknown parameter",
            if (length(unknownParameterNames) != 1) "s were" else "was",
            " removed:\n",
            paste(unknownParameterNames, collapse = "\n"),
            call. = FALSE
        )
        parameters[unknownParameterNames] <- NULL
    }

    present <- defaultParameterNames %in% inputParameterNames

    if (!all(present)) {
        missingParameterNames <- defaultParameterNames[!present]
        warning(
            "Missing parameter",
            if (length(missingParameterNames) != 1) "s were" else "was",
            " replaced with ",
            if (length(missingParameterNames) != 1) "their" else "its",
            " default",
            if (length(missingParameterNames) != 1) "s",
            ":\n",
            paste0(
                missingParameterNames,
                ":",
                defaultHiCDOCParameters[missingParameterNames],
                collapse = "\n"
            ),
            call. = FALSE
        )
        parameters[missingParameterNames] <-
            defaultHiCDOCParameters[missingParameterNames]
    }

    return(parameters)
}

#' @description
#' Returns valid chromosome, condition, or replicate names, from given names
#' or id.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @param names
#' One or several names or ids to look up.
#'
#' @param category
#' The category in which to look up the names or ids. One of "chromosomes",
#' "conditions", "replicates". Defaults to "chromosomes".
#'
#' @return
#' Valid names.
#'
#' @keywords internal
#' @noRd
.validateNames <- function(object, names, category = "chromosomes") {

    validNames <- unique(slot(object, category))

    if (all(names %in% validNames)) return(names)

    if (is.numeric(names) && all(names %in% seq_len(length(validNames)))) {
        return(validNames[names])
    }

    unknown <- names[which(!(names %in% validNames))]

    stop(
        "Unknown ",
        substr(category, 1, nchar(category) - 1),
        if (length(unknown) != 1) "s",
        ": ",
        paste(unknown, collapse = ", "),
        call. = FALSE
    )
}

#' @description
#' Checks that the provided object is a \code{\link{HiCDOCDataSet}}, and that
#' the provided slots are in that object.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @param slots
#' The names of slots to verify. Defaults to NULL.
#'
#' @return
#' Raises an error if the object is not a \code{\link{HiCDOCDataSet}}, or if one
#' of the slots is missing.
#'
#' @keywords internal
#' @noRd
.validateSlots <- function(object, slots = NULL) {

    if (!is(object, "HiCDOCDataSet")) {
        stop(
            "The provided object is not a 'HiCDOCDataSet'.",
            call. = FALSE
        )
    }

    if (
        "interactions" %in% slots && (
            !.hasSlot(object, "interactions") ||
            is.null(slot(object, "interactions")) ||
            nrow(object@interactions) == 0
        )
    ) {
        stop(
            "No interactions found in the 'HiCDOCDataSet'.",
            call. = FALSE
        )
    }

    if (!is.null(slots)) {

        allSlots <- slotNames("HiCDOCDataSet")

        presentSlots <-
            allSlots[
                vapply(
                    allSlots,
                    function(x) {
                        .hasSlot(object, x) && !is.null(slot(object, x))
                    },
                    FUN.VALUE = TRUE
                )
            ]

        missingSlots <- slots[!(slots %in% presentSlots)]

        if ("resolution" %in% missingSlots) {
            stop(
                "Resolution is unknown.\n",
                "This 'HiCDOCDataSet' wasn't built properly.",
                call. = FALSE
            )
        }

        if ("totalBins" %in% missingSlots) {
            stop(
                "Chromosome lengths are unknown.\n",
                "This 'HiCDOCDataSet' wasn't built properly.",
                call. = FALSE
            )
        }

        if ("positions" %in% missingSlots) {
            stop(
                "Positions are unknown.\n",
                "This 'HiCDOCDataSet' wasn't built properly.",
                call. = FALSE
            )
        }

        if (any(c("validReplicates", "validConditions") %in% missingSlots)) {
            stop(
                "Cannot process potentially sparse replicates.\n",
                "First, run 'filterSparseReplicates()' on the object.",
                call. = FALSE
            )
        }

        if ("weakBins" %in% missingSlots) {
            stop(
                "Cannot process potentially weak positions.\n",
                "First, run 'filterWeakPositions()' on the object.",
                call. = FALSE
            )
        }

        compartmentSlots <- c(
            "compartments",
            "concordances",
            "differences",
            "distances",
            "centroids",
            "selfInteractionRatios"
        )

        if (any(compartmentSlots %in% missingSlots)) {
            stop(
                "No compartments found.\n",
                "First, run 'detectCompartments()' on the object.",
                call. = FALSE
            )
        }

        if (length(missingSlots) > 0) {
            stop(
                "Missing slots: ",
                paste(missingSlots, collapse = ", "),
                call. = FALSE
            )
        }
    }
}
