## - .validateParameters ----------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Check object@parameters, return default if NULL
#'
#' @param parameters list of parameters
#'
#' @return list of updated parameters, default from HiCDOCDefaultParameters()
#' if null.
#' @keywords internal
#' @noRd
.validateParameters <- function(parameters) {
    defaultParameterNames <- names(HiCDOCDefaultParameters)
    inputParameterNames <- names(parameters)

    numeric <-
        vapply(
            parameters,
            function(x) is.numeric(x) && length(x) == 1,
            FUN.VALUE = TRUE
        )
    if (!all(numeric)) {
        notNumericNames <- inputParameterNames[!numeric]
        warning(
            "Non-numeric value",
            if (length(notNumericNames) != 1) "s",
            " will be replaced with ",
            if (length(notNumericNames) != 1) "their" else "its",
            " default",
            if (length(notNumericNames) != 1) "s",
            " from 'HiCDOCDefaultParameters':\n",
            paste0(
                notNumericNames,
                ":",
                " ",
                parameters[notNumericNames],
                " -> ",
                HiCDOCDefaultParameters[notNumericNames],
                collapse = "\n"
            ),
            call. = FALSE
        )
        parameters[notNumericNames] <- HiCDOCDefaultParameters[notNumericNames]
    }

    known <- inputParameterNames %in% defaultParameterNames
    if (!all(known)) {
        unknownParameterNames <- inputParameterNames[!known]
        warning(
            "Unknown parameter",
            if (length(unknownParameterNames) != 1) "s",
            " will be removed: ",
            paste(unknownParameterNames, collapse = ", "),
            call. = FALSE
        )
        parameters[unknownParameterNames] <- NULL
    }

    present <- defaultParameterNames %in% inputParameterNames
    if (!all(present)) {
        missingParameterNames <- defaultParameterNames[!present]
        warning(
            "Missing parameter",
            if (length(missingParameterNames) != 1) "s",
            " will be filled with ",
            if (length(missingParameterNames) != 1) "their" else "its",
            " default",
            if (length(missingParameterNames) != 1) "s",
            " from 'HiCDOCDefaultParameters':\n",
            paste0(
                missingParameterNames,
                ":",
                HiCDOCDefaultParameters[missingParameterNames],
                collapse = "\n"
            ),
            call. = FALSE
        )
        parameters[missingParameterNames] <-
            HiCDOCDefaultParameters[missingParameterNames]
    }

    return(parameters)
}

## - .validateNameOrId -----------------------------------------------------------#
## -------------------------------------------------------------------------#
#' Test the existence of a set of identifier
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param id The identifiers to test (string or integer)
#' @param category In which slot should we test the presence of the identifiers ?
#' One of "chromosomes", "conditions", "replicates"
#'
#' @return The chromosome name or an error.
#' @keywords internal
#' @noRd
.validateNameOrId <- function(object, id, category = "chromosomes") {
    names <- unique(slot(object, category))
    if (all(id %in% names)) return(id)
    if (is.numeric(id) && all(id %in% seq_len(length(names)))) return(names[id])
    unknown <- id[which(!(id %in% names))]
    stop(
        "Unknown ",
        substr(category, 1, nchar(category) - 1),
        if (length(unknown) != 1) "s",
        ": ",
        paste(unknown, collapse = ", "),
        call. = FALSE
    )
}

## - .validateSlots --------------------------------------------------------#
## --------------------------------------------------------------------------#
#' Test the existence of slots
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param slots Character vector, names of slots to verify. Default to NULL.
#' If NULL, check only for the class of \code{object}
#'
#' @return An error if the object is not a HiCDOCDataSet object or a slot is
#' missing.
#' @keywords internal
#' @noRd
.validateSlots <- function(object, slots = NULL) {
    if (!is(object, "HiCDOCDataSet")) {
        stop(
            "The provided object is not an HiCDOC object.",
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
        if ("interactions" %in% missingSlots) {
            stop(
                "No interactions found. ",
                "Provide an HiCDOC object with interactions.",
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
        if (any(missingSlots %in% compartmentSlots)) {
            stop(
                "No compartments found. ",
                "Call 'detectCompartments()' first.",
                call. = FALSE
            )
        }
    }
}
