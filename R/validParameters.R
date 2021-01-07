## - testChromosome -----------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Test the existence of a given chromosome.
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param chromosomeId The condition name, or an error.
#'
#' @return The chromosome name or an error.
testChromosome <- function(object, chromosomeId) {
  if (chromosomeId %in% object@chromosomes) {
    return(chromosomeId)
  }
  if (is.numeric(chromosomeId) == TRUE &&
    chromosomeId %in% seq_len(length(object@chromosomes))) {
    return(object@chromosomes[chromosomeId])
  }
  stop(paste("Unknown chromosome:", chromosomeId), call. = FALSE)
}

## - testCondition -----------------------------------------------------------#
## ---------------------------------------------------------------------------#
#' Test the existence of a given condition in a HiCDOC object
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param conditionId    A character or numeric value
#'
#' @return The condition name, or an error.
testCondition <- function(object, conditionId) {
  if (conditionId %in% object@conditions) {
    return(conditionId)
  }
  uniquecond <- unique(object@conditions)
  if (is.numeric(conditionId) && conditionId %in% seq_len(length(uniquecond))) {
    return(uniquecond[conditionId])
  }
  stop(paste("Unknown condition:", conditionId), call. = FALSE)
}

## - testSlotsHiCDOC -------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Test the existence of slots in HiCDOCDataSet object
#'
#' @param object A \code{HiCDOCDataSet} object.
#' @param slots Character vector, names of slots to verify. Default to NULL.
#' If NULL, check only for the class of \code{object}
#'
#' @return An error if the object is not a HiCDOCDataSet object or a slot is
#' missing.
testSlotsHiCDOC <- function(object, slots = NULL) {
  if (!is(object, "HiCDOCDataSet")) {
    stop("The object provided is not from class HiCDOCDataSet", call. = FALSE)
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
            "\nPlease provide an HiCDOC object, with filled interactions."
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

#' testNbCores
#' Test if the nbCores parameter is correct, for parallel computations.
#'
#' @param nbCores integer. The number of cores to use in parallel::mcmapply.
#'
#' @return integer
testNbCores <- function(nbCores) {
  nbAvail <- parallel::detectCores()
  if (is.null(nbCores)) {
    stop(paste0(
      "Please give a nbCores parameters, suggested: ",
      nbAvail - 1, "."
    ),
    call. = FALSE
    )
  } else {
    nbAvail <- parallel::detectCores()
    if (nbCores > nbAvail) {
      message(paste0("nbCores is too big, used = ", nbAvail - 1, "."))
      nbCores <- nbAvail - 1
    }
  }
  return(nbCores)
}
