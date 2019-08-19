##- Parse input data ---------------------------------------------------------#
##----------------------------------------------------------------------------#

#' @export
parseInteractionMatrix3Columns <- function(object) {
    object@interactionMatrix <- read_tsv(
        file = object@inputMatrixPath,
        comment = '#',
        col_names = TRUE
    )
    if (colnames(object@interactionMatrix)[[1]] != "chromosome") {
        stop("First column of the input matrix should be named 'chromosome'.")
    }
    if (colnames(object@interactionMatrix)[[2]] != "position 1") {
        stop("First column of the input matrix should be named 'position 1'.")
    }
    if (colnames(object@interactionMatrix)[[3]] != "position 2") {
        stop("First column of the input matrix should be named 'position 2'.")
    }
    object@replicates <- colnames(object@interactionMatrix)[4:ncol(object@interactionMatrix)]
    object@conditions <- rep(NA, length(object@replicates))
    object@conditions[grep("replicate 1", object@replicates)] <- 1
    object@conditions[grep("replicate 2", object@replicates)] <- 2
    if (any(is.na(object@conditions))) {
        stop("Column names of interaction matrix is not well formed.\n",
             "It should be 'chromosome', 'position 1', 'position 2',
             'replicate 1.1', 'replicate 1.2', etc.", call. = FALSE)
    }
    object@interactionMatrix %<>%
        gather(object@replicates, key = "replicate", value = "value") %>%
        separate(replicate, c(NA, "condition", "replicate")) %>%
        mutate(chromosome = factor(chromosome),
               condition = factor(condition),
               replicate = factor(replicate))
    return(object)
}


##- Parse a single cool Matrix -----------------------------------------------#
##----------------------------------------------------------------------------#

h5readCatch <- function(file = file, name = name) {
    return(tryCatch(h5read(file = file, name = name),
                    error = function(e) {
                        stop(paste0("The file '",
                                    file,
                                    "' does not seem to have the .cool format:",
                                    "path '",
                                    name ,
                                    "' is missing."),
                             call. = FALSE)
                    }))
}

parseCoolMatrix <- function(fileName) {
    bins <- tibble(chromosome = factor(h5readCatch(file = fileName,
                                            name = "bins/chrom")),
                   start      = h5readCatch(file = fileName,
                                            name = "bins/start"),
                   end        = h5readCatch(file = fileName,
                                            name = "bins/end"))
    step <- bins$end - bins$start
    pc <- sum(step == max(step)) / length(step)
    if (length(step) < 0.9) {
        stop(paste0("Cannot parse cool file ",
                      fileName,
                      "fixed width only."))
    }
    step <- max(step)

    bins %<>%
        select(-(end)) %>%
        rename(position = start) %>%
        rowid_to_column("id") %>%
        mutate(id = id - 1)

    data <- tibble(id1 = h5readCatch(file = fileName,
                                     name = "pixels/bin1_id"),
                   id2 = h5readCatch(file = fileName,
                                     name = "pixels/bin2_id"),
                   value = h5readCatch(file = fileName,
                                       name = "pixels/count")) %>%
        left_join(bins, by = c("id1" = "id")) %>%
        rename(`chromosome 1` = chromosome,
               `position 1` = position) %>%
        select(-id1) %>%
        left_join(bins, by = c("id2" = "id")) %>%
        rename(`chromosome 2` = chromosome,
               `position 2` = position) %>%
        select(-id2) %>%
        filter(`chromosome 1` == `chromosome 2`) %>%
        select(-`chromosome 2`) %>%
        rename(chromosome = `chromosome 1`) %>%
        select(chromosome, `position 1`, `position 2`, value)

    return(data)
}

#' @export
parseInteractionMatrixCool <- function(object) {
    matrices <- lapply(object@inputMatrixPath, parseCoolMatrix)
    replicates <- lapply(seq_along(object@conditions),
                         function(x) {
                             length(which(object@conditions[seq(x)] ==
                                              object@conditions[[x]]))
                         })
    replicates <- paste0("replicate ", object@conditions, ".", replicates)
    for (i in seq_along(matrices)) {
        matrices[[i]]$replicate <- factor(replicates[[i]])
    }
    object@interactionMatrix <- bind_rows(matrices) %>%
        mutate(chromosome = factor(chromosome)) %>%
        mutate(replicate = factor(replicate))

    return(object)
}
