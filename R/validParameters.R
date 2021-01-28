## - checkParameters ----------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Check object@parameters, return default if NULL
#'
#' @param parameters list of parameters
#'
#' @return list of updated parameters, default from HiCDOCDefaultParameters()
#' if null.
checkParameters <- function(parameters) {
    defaultParNames <- names(HiCDOCDefaultParameters)
    currentParNames <- names(parameters)
    # Test for numerical values only
    testNum <- vapply(parameters, 
                      function(x) is.numeric(x) && length(x) == 1, FALSE)
    if(!all(testNum)){
        warning(paste0("non numeric values found in parameters,",
                "replaced by HiCDOCDefaultParameters corresponding values"),
                call. = FALSE)
        notNum <- currentParNames[which(!testNum)]
        parameters[notNum] <- HiCDOCDefaultParameters[notNum]
    }
    
    # Test for correct names
    includePar <- currentParNames %in% defaultParNames
    if(!(all(includePar))){
        unknownPAr <- currentParNames[!includePar]
        warning(paste0("unknown parameters: ", paste(unknownPAr, collapse=", "),
                       " they will be removed."),
                       call. = FALSE)
        parameters[unknownPAr] <- NULL
    }
    
    includePar <- defaultParNames %in% currentParNames
    if(!(all(includePar))){
        knownPAr <- defaultParNames[includePar]
        warning(paste0("Unfixed parameters: ", paste(unknownPAr, collapse=", "),
                       " they will be removed."),
                       call. = FALSE)
        tempPar <- HiCDOCDefaultParameters
        tempPar[knownPar] <- parameters[knownPAr]
        parameters <- tempPar
    }
    
    return(parameters)
}

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

