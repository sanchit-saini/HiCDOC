#' Look for the weak positions for a given chromosome
#'
#' The function identifies the weak positions from the interaction matrix
#' for a given chromosome., in a recursive way.
#'
#' @param inter_chr Interaction matrix for the chromosome
#' @param totalBins_chr Number of total bins for this chromosome
#' @param binsize Size of the bins in the HiCDOC object
#' @param replicates Levels of replicates. Given from object@conditions
#' @param conditions Levels of conditions (repeated by replicates).
#' Given from object@conditions
#' @param threshold Threshold for the selection of "empty". The function will
#' return the bin if the mean of value if inferior to `threshold`,
#' for at least one condition and a replicate (after we reconstruct the 0 values).
#' Default to 0.
#'
#' @return
#'
#' @examples
#' object <- HiCDOCExample()
#' chr <- object@chromosomes[[1]]
#' filterWeakChr(object@interactions[object@interactions$chromosome == chr,],
#'              object@totalBins[[chr]],
#'              object@binSize,
#'              object@replicates,
#'              object@conditions,
#'              object@filterThreshold)
filterWeakChr <- function(interChr, totalBinsChr, binsize, replicates, conditions, threshold=0){

  # Initialization
  chr <- interChr[1,"chromosome"] %>% pull() %>% as.character()
  nbNewEmptyBins <- 1
  nbRemovedBins <- 0

  fullbindata <- tibble("replicate"= factor(replicates), "condition" = factor(conditions)) %>% 
    tidyr::expand(tidyr::nesting(replicate, condition), 
                  "position" = as.integer((seq_len(totalBinsChr)-1) * binsize))
  interChr <- interChr %>% filter(value>0)

  # Recursive removal of bins - deleting a bin can create a new weak bin. 
  while(nbNewEmptyBins>0 & nbRemovedBins<=totalBinsChr){
    existing <- interChr %>% 
      tidyr::pivot_longer(cols = c(`position.1`, `position.2`), names_to="pos_1or2", 
                          names_prefix="position.", values_to="position") %>%
      select(-pos_1or2) 
    
    fullbindatatest <- fullbindata %>% 
      dplyr::left_join(existing, by = c("replicate", "condition", "position")) %>%  
      mutate(value = ifelse(is.na(value)==T,0,value)) %>%
      group_by(replicate, condition, position) %>%
      mutate(mean = mean(value)) %>%
      filter(mean <= threshold)

    weakposChr <- fullbindatatest %>%
      pull(position) %>%
      unique() %>%
      sort()

    nbNewEmptyBins <- length(weakposChr) - nbRemovedBins

    # Remove the positions in rows and columns in empty bins
    if(nbNewEmptyBins>0) {
      interChr[interChr$position.1 %in% weakposChr,]$value <- 0
      interChr[interChr$position.2 %in% weakposChr,]$value <- 0
      nbRemovedBins <- nbRemovedBins + nbNewEmptyBins
      interChr <- interChr %>% filter(value>0)
    }
  }
  message(paste0("Chromosome ", chr, ": ",
                     length(weakposChr), " position(s) removed"))
  return(list("pos"= weakposChr, "interactions" = interChr ))
}


#' Remove the weak bins of an HiCDOC object
#' 
#' The function indentifies and remove the weak bins of an HiCDOC object, 
#' chromosome by chromosome. 
#' To be kept, the bins (positions) must have a mean value greater than 
#' object@filterThreshold, on all replicates and all conditions. 
#' The mean is computed on the row of the reconstructed full interaction 
#' matrix for 1 chromosome, 1 condition and 1 replicate.
#' 
#' @param object A HiCDOC object
#'
#' @return A \code{HicDOCExp} object with a reduced \code{interactions} slot, 
#' and the weak bins identified by chromosomes in \code{object@weakBins}.
#' @export
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- filterWeakPositions(exp)
filterWeakPositions <- function(object) {
  weakPositions <- lapply(object@chromosomes,
                          function(chr) filterWeakChr(object@interactions[object@interactions$chromosome == chr,],
                                                     object@totalBins[[chr]],
                                                     object@binSize,
                                                     object@replicates,
                                                     object@conditions,
                                                     object@filterThreshold)
  )
  names(weakPositions) <- object@chromosomes

  weakBins <- weakPositions %>% purrr:::map("pos")
  weakBins <- lapply(weakBins, function(x) x / object@binSize + 1)
  nullweak <- vapply(weakBins, function(x) length(x) == 0, FUN.VALUE = TRUE)
  weakBins[nullweak] <- list(NULL)
  
  interactions <- weakPositions %>% purrr::map_dfr("interactions")

  # Affecting new values to the object
  object@weakBins <- weakBins
  object@interactions <- interactions

  nbweak <- sum(vapply(weakBins, length, c(0)))

  message(cat(paste0("Chromosome ", names(weakBins), ": ",
                       lapply(weakBins, length), " position(s) removed"), sep="\n"))
  message("Removed ",
          nbweak,
          " position",
          if (nbweak > 1)
            "s")

  return (object)
}
