##- Helper function for check input parameters -------------------------------#
##----------------------------------------------------------------------------#
checkParameters <- function(value) {
    # TODO
}

#' @export
removeSmallChromosomes <- function(object) {
    bigChromosomes <- object@interactionMatrix %>%
        group_by(chromosome, replicate) %>%
        summarise(nCounts = n()) %>%
        summarise(meanNCounts = mean(nCounts)) %>%
        filter(meanNCounts >= object@minNCounts) %>%
        pull(chromosome)
    message("Keeping ", length(bigChromosomes), " chromosomes.")
    bigChromosomes <- sort(bigChromosomes)
    object@chromosomes <- droplevels(bigChromosomes)
    object@interactionMatrix %<>%
        filter(chromosome %in% bigChromosomes) %>%
        mutate(chromosome = droplevels(chromosome))
    return(object)
}
