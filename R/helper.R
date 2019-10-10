##- Helper function for check input parameters -------------------------------#
##----------------------------------------------------------------------------#
checkParameters <- function(value) {
    # TODO
}

makeFullMatrix <- function(data) {
    data %>% filter(`position 1` != `position 2`) %>%
        rename(tmp = `position 1`, `position 1` = `position 2`) %>%
        rename(`position 2` = tmp) %>%
        bind_rows(data)
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
