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
