#' @export
normalizeCyclicLoess <- function(object) {
  matrix <- object@interactionMatrix %>%
    spread(replicate, value)
  matrices <- lapply(object@replicates,
                     function(x) {
                       matrix[, c("chromosome", "position 1", "position 2", x)]
                     })
  matrices <- object@interactionMatrix %>%
    split(.$replicate) %>%
    map(~ select(.x, -replicate))
  hicexp <- make_hicexp(
    data_list = matrices,
    groups = object@conditions,
    remove.regions = NULL
  )
  normalized <- cyclic_loess(hicexp, parallel = TRUE)
  output <- as_tibble(hic_table(normalized)) %>%
    select(-D)
  colnames(output) <- c("chromosome",
                        "position 1",
                        "position 2",
                        object@replicates)
  object@interactionMatrix <- output %>%
    gather(object@replicates, key = "replicate", value = "value")
  return(object)
}

# Split joint matrix into one matrix per replicate
# matrices = vector(mode = 'list', length = length(replicates))
# for (r in 1:length(replicates)) {
#   matrices[[r]] = input[,c(1,2,3,r+3)]
# }

# Create hicexp object

# Normalize

# Format normalized matrix
# normalized_df = as.data.frame(hic_table(normalized))
# normalized_df = normalized_df[,!(names(normalized_df) %in% 'D')]
# names(normalized_df) = c(
#   'chromosome', 'position 1', 'position 2',
#   paste(rep('replicate', length(replicates)), replicates)
# )
#
# message("Normalize cyclic loess writing to output")
#
# # Write to output
# #writeLines(comments, args$o)
#
# write.table(
#   normalized_df,
#   file = args$o,
#   quote = FALSE,
#   sep = '\t',
#   row.names = FALSE
# )
#
# message("Done cyclic loess")
