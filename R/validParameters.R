##- testChromosome -----------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Test the existence of a given chromosome.
#'
#' @param object        A \code{HiCDOCExp} object.
#' @param chromosomeId  A chromosome.
#'
#' @return A chromosome or an error.
testChromosome <- function(object, chromosomeId) {
  if (chromosomeId %in% object@chromosomes) return (chromosomeId)
  if (chromosomeId %in% seq_len(length(object@chromosomes)))
    return (object@chromosomes[as.numeric(chromosomeId)])
  stop(paste("Unknown chromosome:", chromosomeId), call. = FALSE)
}

##- testSlotsHiCDOCExp -------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Test the existence of slots in HiCDOCExp object
#'
#' @param object  A \code{HiCDOCExp} object.
#' @param slots   Names of slots to verify.
#'
#' @return An error if the object is not a HiCDOCExp object or a slot is
#' missing.
testSlotsHiCDOCExp <- function(object, slots = NULL) {
  if (class(object) != "HiCDOCExp")
    stop("The object provided is not from class HiCDOCExp", call. = FALSE)

  if (!is.null(slots)) {
    existing <- slotNames(object)
    existing <- existing[vapply(
      existing,
      function(x) !is.null(slot(object, x)),
      TRUE
    )]
    missing <- slots[!(slots %in% existing)]
    if (length(missing) > 0) {
      missingInteractions = FALSE
      missingCompartments = FALSE
      compartmentSlots = c(
        "compartments",
        "concordances",
        "differences",
        "distances",
        "centroids",
        "diagonalRatios"
      )
      for (slot in missing) {
        if (slot == "interactions") missingInteractions = TRUE
        if (slot %in% compartmentSlots) missingCompartments = TRUE
      }
      stop(
        paste0(
          "Missing slot(s): ",
          paste(missing, collapse = ", "),
          if (missingInteractions)
            "\nPlease provide an HiCDOC object, with filled interactions.",
          if (missingCompartments)
            "\nPlease run 'detectCompartments' first."
        ),
        call. = FALSE
      )
    }
  }
}
