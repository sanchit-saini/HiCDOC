buildABComparison <- function(object) {
    diagonal <- object@interactionMatrix %>%
        filter(`position 1` == `position 2`) %>%
        rename(position = `position 1`) %>%
        select(-`position 2`) %>%
        group_by(chromosome, position, condition) %>%
        summarise(diagValue = median(value)) %>%
        ungroup()
    offDiagonal <- object@interactionMatrix %>%
        makeFullMatrix() %>%
        filter(`position 1` != `position 2`) %>%
        group_by(chromosome, `position 1`, condition, replicate) %>%
        summarise(value = mean(value)) %>%
        summarise(offDiagValue = median(value)) %>%
        ungroup() %>%
        rename(position = `position 1`)
    full_join(diagonal, offDiagonal,
              by = c("chromosome", "position", "condition")) %>%
        right_join(object@compartments,
                   by = c("chromosome", "position", "condition")) %>%
        rename(compartment = value) %>%
        replace_na(list(diagValue = 0, offDiagValue = 0)) %>%
        mutate(diffValue = diagValue - offDiagValue) %>%
        select(-c(position, condition, diagValue, offDiagValue)) %>%
        mutate(chromosome = factor(chromosome, levels = object@chromosomes))
}

#' @export
predictAB <- function(object) {
    if (is.null(object@interactionMatrix)) {
        stop(paste0("Interaction matrix is not loaded yet.  ",
                    "Please provide a matrix first."))
    }
    if (is.null(object@compartments)) {
        stop(paste0("Compartments are not computed.  ",
                    "Please run 'detectConstrainedKMeans' first."))
    }

    compartments <- buildABComparison(object) %>%
        group_by(chromosome, compartment) %>%
        summarise(value = median(diffValue)) %>%
        ungroup() %>%
        spread(key = compartment, value = value, fill = 0) %>%
        mutate(A = if_else(`1` >= `2`, 1, 2)) %>%
        select(-c(`1`, `2`))
    object@compartments %<>%
        left_join(compartments, by = "chromosome") %>%
        mutate(value = factor(if_else(value == A, "A", "B"))) %>%
        select(-c(A))
    object@concordances %<>%
        left_join(compartments, by = "chromosome") %>%
        mutate(change = if_else(A == 1, 1, -1)) %>%
        mutate(value = change * value) %>%
        select(-c(A, change))
    object@distances %<>%
        left_join(compartments, by = "chromosome") %>%
        mutate(cluster = factor(if_else(cluster == A, "A", "B"))) %>%
        select(-c(A))
    return(object)
}
