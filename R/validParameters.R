## - checkParameters ----------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Check object@parameters, return default if NULL
#'
#' @param parameters list of parameters
#'
#' @return list of updated parameters, default from HiCDOCDefaultParameters()
#' if null.
#' @keywords internal
#' @noRd
checkParameters <- function(parameters) {
    defaultParNames <- names(HiCDOCDefaultParameters)
    currentParNames <- names(parameters)
    # Test for numerical values only
    testNum <- vapply(
        parameters,
        function(x) is.numeric(x) && length(x) == 1, FALSE
    )
    if (!all(testNum)) {
        warning(paste0(
            "non numeric values found in parameters,",
            "replaced by HiCDOCDefaultParameters corresponding values"
        ),
        call. = FALSE
        )
        notNum <- currentParNames[which(!testNum)]
        parameters[notNum] <- HiCDOCDefaultParameters[notNum]
    }

    # Test for correct names
    includePar <- currentParNames %in% defaultParNames
    if (!(all(includePar))) {
        unknownPAr <- currentParNames[!includePar]
        warning(paste0(
            "unknown parameters: ",
            paste(unknownPAr, collapse = ", "),
            " they will be removed."
        ),
        call. = FALSE
        )
        parameters[unknownPAr] <- NULL
    }

    includePar <- defaultParNames %in% currentParNames
    if (!(all(includePar))) {
        knownPAr <- defaultParNames[includePar]
        warning(paste0(
            "Unfixed parameters: ",
            paste(unknownPAr, collapse = ", "),
            " they will be removed."
        ),
        call. = FALSE
        )
        tempPar <- HiCDOCDefaultParameters
        tempPar[knownPar] <- parameters[knownPAr]
        parameters <- tempPar
    }

    return(parameters)
}

## - testValidId -----------------------------------------------------------#
## -------------------------------------------------------------------------#
#' Test the existence of a set of identifier
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param id The identifiers to test (string or integer)
#' @param what In what slot should we test the presence of the identifiers ?
#' One of "chromosome", "condition", "replicate"
#'
#' @return The chromosome name or an error.
#' @keywords internal
#' @noRd
testValidId <- function(object, id, what = "chromosomes") {
    valid <- unique(slot(object, what))
    if (all(id %in% valid)) {
        return(id)
    }
    if (is.numeric(id) == TRUE && all(id %in% seq_len(length(valid)))) {
        return(valid[id])
    }
    if (length(id) == 1) {
        what <- substr(what, 1, nchar(what) - 1)
        stop(paste0("Unknown ", what, ": ", id), call. = FALSE)
    } else {
        stop(paste0("Unknown ", what, ": ", paste(id, collapse = ", ")),
            call. = FALSE
        )
    }
}

## - testSlotsHiCDOC --------------------------------------------------------#
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
testSlotsHiCDOC <- function(object, slots = NULL) {
    if (!is(object, "HiCDOCDataSet")) {
        stop("The object provided is not from class HiCDOCDataSet", 
             call. = FALSE)
    }

    if (!is.null(slots)) {
        existing <- slotNames("HiCDOCDataSet")
        existing <- existing[vapply(
            existing,
            function(x) .hasSlot(object, x) && !is.null(slot(object, x)),
            TRUE
        )]
        missing <- slots[!(slots %in% existing)]
        if (length(missing) > 0) {
            missingInteractions <- FALSE
            missingCompartments <- FALSE
            compartmentSlots <- c(
                "compartments",
                "concordances",
                "differences",
                "distances",
                "centroids",
                "diagonalRatios"
            )
            for (slot in missing) {
                if (slot == "interactions") missingInteractions <- TRUE
                if (slot %in% compartmentSlots) missingCompartments <- TRUE
            }
            stop(
                paste0(
                    "Missing slot(s): ",
                    paste(missing, collapse = ", "),
                    if (missingInteractions) {
                        "\nPlease provide an HiCDOC object, with interactions."
                    },
                    if (missingCompartments) {
                        "\nPlease run 'detectCompartments' first."
                    }
                ),
                call. = FALSE
            )
        }
    }
}
