##- Parse input data ---------------------------------------------------------#
##----------------------------------------------------------------------------#

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
    object@chromosomes <- unique(object@interactionMatrix$chromosome)
    object@replicates <- colnames(object@interactionMatrix)[4:ncol(object@interactionMatrix)]
    object@nReplicates <- length(object@replicates)
    object@conditions <- rep(NA, length(object@replicates))
    object@conditions[grep("replicate 1", object@replicates)] <- 1
    object@conditions[grep("replicate 2", object@replicates)] <- 2
    if (any(is.na(object@conditions))) {
        stop("Column names of interaction matrix is not well formed.\n",
             "It should be 'chromosome', 'position 1', 'position 2',
             'replicate 1.1', 'replicate 1.2', etc.", call. = FALSE)
    }
    object@binSize <- min(object@interactionMatrix$`position 1`
                          [object@interactionMatrix$`position 1` > 0])
    object@interactionMatrix %<>%
        gather(object@replicates, key = "replicate", value = "value")
    return(object)
}
