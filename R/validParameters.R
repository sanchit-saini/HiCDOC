#' Return valid chromosome as string
#' 
#' Test if the chromosome entred exists. 
#' If a numerical value is verified, and exists, return the 
#' corresponding character in object@chromosome
#'
#' @param object An HiCDOC object
#' @param chrId Numerical or character value
#'
#' @return
testchromosome <- function(object, chrId){
  if(chrId %in% object@chromosomes) {
    chr <- chrId      
  } else {
    if(is.numeric(chrId) && chrId %in% seq_len(length(object@chromosomes))){
      chr <- object@chromosomes[chrId]
    } else {
      stop("Unknown chromosomeId") 
    }
  }
  return(chr)
}


#' Test the existence of slots in HiCDOCExp object
#'
#' @param object a R object
#' @param slots the names of the slots to test the existence
#'
#' @return An error if the object is not a HiCDOCExp object or
#' a slot is missing, NULL otherwise
testSlotsHiCDOCExp <- function(object, slots = NULL){
  if(class(object) != "HiCDOCExp"){
    stop("The object provided is not from class HiCDOCExp", call.=FALSE)
  }
  if(is.null(slots) == F){
    existing <- slotNames(object)
    existing <- existing[vapply(existing, function(x) !is.null(slot(object, x)), TRUE)]
    missing <- slots[!(slots %in% existing)]
    if(length(missing) > 0){
      for(i in missing){
        if(i == "interactions"){
          message(paste0(
            "Interaction matrix is not loaded yet. ",
            "Please provide an HiCDOC object, with a non empty slot interactions."
          ))
        }
        if(i %in% c("concordances", "compartments", "differences")){
          message(paste0(
            i, " are not computed. Please run 'detectCompartments' first."
          ))
        }
      }
      stop(paste("missing slot(s):", paste(missing, collapse = ", ")), call. = FALSE)
    }
  } else {
    return(NULL)
  }
}
