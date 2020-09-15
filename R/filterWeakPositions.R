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
#' exp <- HiCDOCExample()
#' chr <- exp@chromosomes[[1]]
#' fun_weak_chr(object@interactions[object@interactions$chromosome == chr,],
#'              object@totalBins[[chr]],
#'              object@binSize,
#'              object@replicates,
#'              object@conditions,
#'              object@filterThreshold)
#' @export
fun_weak_chr <- function(inter_chr, totalBins_chr, binsize, replicates, conditions, threshold=0){

  # Initialization
  chr <- inter_chr[1,"chromosome"]
  weakbins_chr <- NULL
  nbNewEmptyBins <- 1
  nbRemovedBins <- 0
  fullbindata <- tibble("replicate"= factor(replicates), "condition" = factor(conditions)) %>%
    tidyr::expand(tidyr::nesting(replicate, condition), "position" = (seq_len(totalBins_chr)-1) * binsize)
  inter_chr <- inter_chr %>% filter(value>0)

  # Recursive removal of bins - deleting a bin can create a new weak bin.
  while(nbNewEmptyBins>0 & nbRemovedBins<=totalBins_chr){
    existing <- inter_chr %>%
      tidyr::pivot_longer(cols = c(`position.1`, `position.2`), names_to="pos_1or2", names_prefix="position.", values_to="position") %>%
      select(-pos_1or2)

    fullbindatatest <- fullbindata %>%
      dplyr::left_join(existing, by = c("replicate", "condition", "position")) %>%
      mutate(value = ifelse(is.na(value)==T,0,value)) %>%
      group_by(replicate, condition, position) %>%
      mutate(mean = mean(value)) %>%
      filter(mean <= threshold)

    weakpos_chr <- fullbindatatest %>%
      pull(position) %>%
      unique() %>%
      sort()

    nbNewEmptyBins <- length(weakpos_chr) - nbRemovedBins

    # Remove the positions in rows and columns in empty bins
    if(nbNewEmptyBins>0) {
      inter_chr[inter_chr$position.1 %in% weakpos_chr,]$value <- 0
      inter_chr[inter_chr$position.2 %in% weakpos_chr,]$value <- 0
      nbRemovedBins <- nbRemovedBins + nbNewEmptyBins
      inter_chr <- inter_chr %>% filter(value>0)
    }
  }
  return(list("pos"= weakpos_chr, "interactions" = inter_chr ))
}


#' Indentifies and remove the weak bins of an HiCDOC object, by chromosome.
#'
#' @param object A \code{HicDOCExp} object.
#'
#' @return A \code{HicDOCExp} object, with a mark on weak bins.
#' @export
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- filterWeakPositions(exp)
filterWeakPositions <- function(object) {
  weakPositions <- lapply(object@chromosomes,
                          function(chr) fun_weak_chr(object@interactions[object@interactions$chromosome == chr,],
                                                     object@totalBins[[chr]],
                                                     object@binSize,
                                                     object@replicates,
                                                     object@conditions,
                                                     object@filterThreshold)
  )
  names(weakPositions) <- object@chromosomes

  weakBins <- weakPositions %>% map("pos")
  weakBins <- lapply(weakBins, function(x) x / object@binSize + 1)
  nullweak <- vapply(weakBins, function(x) length(x) == 0, FUN.VALUE = TRUE)
  weakBins[nullweak] <- list(NULL)

  interactions <- weakPositions %>% purrr::map("interactions") %>% bind_rows()

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
