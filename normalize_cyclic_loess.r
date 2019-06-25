#!/usr/bin/env Rscript

library('argparse')
library('multiHiCcompare')

parser = ArgumentParser(
  description = 'Normalize technical biases with cyclic loess'
)
parser$add_argument('-i', required=TRUE, help='Input matrix')
parser$add_argument('-o', required=TRUE, help='Output matrix')
args = parser$parse_args()

# Read input
input = read.table(
  file = args$i,
  sep = '\t',
  comment.char = '#',
  header = TRUE
)

# Save comment lines before header
comments = c()
lines = file(args$i, 'r')
while (TRUE) {
  line = readLines(lines, 1)
  if (grepl('^#', line)) {
    comments = c(comments, line)
  }
  else {
    break
  }
}

# Determine replicates
replicates = as.vector(sapply(names(input[,4:ncol(input)]), function(x) {
  gsub('^\\w*\\.', '', x)
}))

# Determine conditions
groups = as.numeric(sapply(replicates, function(x) {
  gsub('\\..*', '', x)
}))

# Split joint matrix into one matrix per replicate
matrices = vector(mode = 'list', length = length(replicates))
for (r in 1:length(replicates)) {
  matrices[[r]] = input[,c(1,2,3,r+3)]
}

# Create hicexp object
hicexp = make_hicexp(
  data_list = matrices,
  groups = groups,
  remove.regions = NULL
)

# Normalize
normalized = cyclic_loess(hicexp, parallel = TRUE)

# Format normalized matrix
normalized_df = as.data.frame(hic_table(normalized))
normalized_df = normalized_df[,!(names(normalized_df) %in% 'D')]
names(normalized_df) = c(
  'chromosome', 'position 1', 'position 2',
  paste(rep('replicate', length(replicates)), replicates)
)

# Write to output
writeLines(comments, args$o)

write.table(
  normalized_df,
  file = args$o,
  quote = FALSE,
  sep = '\t',
  row.names = FALSE,
  append = TRUE
)
