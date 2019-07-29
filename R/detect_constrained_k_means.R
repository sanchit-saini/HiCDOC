#!/usr/bin/env Rscript

library('argparser')
library('tibble')
library('dplyr')
library('tidyr')
library('magrittr')
library('readr')
library('ggplot2')
library('ggExtra')
library('purrr')
library('reshape2')
library('Rcpp')

euclideanDistance <- function(p1, p2) {
   sqrt(sum((p1 - p2)^2))
}

concordanceFunction <- function(x, p1, p2) {
  eps <- 1e-10
  log2((euclideanDistance(x, p1) + eps) / (euclideanDistance(x, p2) + eps))
}


parser <- arg_parser("Detect compartments using constrainted k-means")

parser <- add_argument(parser = parser, arg = "--input",
                       help = "Input matrix")
parser <- add_argument(parser = parser, arg = "--compartments",
                       help = "Compartments")
parser <- add_argument(parser = parser, arg = "--distance1",
                       help = "Distances to the centroid 2")
parser <- add_argument(parser = parser, arg = "--distance2",
                       help = "Distances to the centroid 2")
parser <- add_argument(parser = parser, arg = "--concordance",
                       help = "Concordance")
argv <- parse_args(parser)

sourceCpp("/home/mazytnicki/Desktop/Projects/Kurylo/HiCDOC/constrainedClustering.cpp")

distance_files = c(argv$distance1, argv$distance2)

# Read input
input <- read_tsv(
  file = argv$input,
  comment = '#',
  col_names = TRUE
)

# Save comment lines before header
# comments = c()
# lines = file(argv$input, 'r')
# while (TRUE) {
#   line = readLines(lines, 1)
#   if (grepl('^#', line)) {
#     comments = c(comments, line)
#   }
#   else {
#     break
#   }
# }

if (any(input$`position 1` > input$`position 2`)) {
  stop("Input is not well formed")
}

step <- min(input$`position 1`[input$`position 1` > 0])
input %<>%
  mutate(bin1 = `position 1` / step + 1) %>%
  mutate(bin2 = `position 2` / step + 1)

head(input)

replicate_names <- colnames(input)[4:(ncol(input) - 2)]
n_replicates <- length(replicate_names)

# mustLinkSeed <- asplit(rbind(t(combn(grep("replicate 1", colnames(input)), 2)),
#                       t(combn(grep("replicate 2", colnames(input)), 2))) - 4, 1)
mustLinkSeed <- list(rep1 = grep("replicate 1", colnames(input)) - 4,
                     rep2 = grep("replicate 2", colnames(input)) - 4)


distances <- tibble()
compartments <- tibble()
concordances <- tibble()

for (chr in unique(input$chromosome)) {
#chr <- 17
  inputChromosome <- filter(input, chromosome == chr)
  inputReplicate  <- tibble(bin1 = inputChromosome$bin1,
                            bin2 = inputChromosome$bin2,
                            data = 0)
  positions <- unique(inputChromosome$`position 1`)
  n <- max(inputChromosome$bin1, inputChromosome$bin2)
  bigmat <- matrix(nrow = 0, ncol = n)
  message(paste0("Chromosome ", chr, ", of dim. ", n))
  for (col in 4:(ncol(input) - 2)) {
    message(paste0("  Replicate ", colnames(input)[[col]]))
    inputReplicate$data <- inputChromosome[[col]]
    mat <- acast(inputReplicate, bin1 ~ bin2, fill = 0, value.var = "data")
    mat <- mat + t(mat) - diag(diag(mat))
    if (!isSymmetric(mat)) {
      stop(paste0("Matrix ", chr, "/", colnames(input)[[col]], " is not symmetric"))
    }
    bigmat <- rbind(bigmat, mat)
  }
  mustLink <- do.call("rbind", lapply(mustLinkSeed, function(x) { matrix(rep(x, n)*n + rep(0:(n - 1), each = length(x)), nrow = n, byrow = TRUE) } ))
  #mustLink <- do.call("rbind", lapply(mustLinkSeed, function(x) { matrix(rep(x, n)*n + rep(1:n, each=length(x)), nrow=n, byrow=TRUE) } ))
  #clusters <- ckmeans(bigmat, 2, mustLink, c())
  clusters <- constrainedClustering(bigmat, mustLink, 0.001, 10) + 1
  centroids <- lapply(1:2,
                      function(i) {
                        colMeans(bigmat[ which(clusters == i),  ])
                      }
                     )
  compartments %<>% bind_rows(tibble(
                                chromosome = chr,
                                position   = positions,
                                replicate  = "replicate 1",
                                compartment = head(clusters, n))) %>%
                    bind_rows(tibble(
                                chromosome = chr,
                                position   = positions,
                                replicate  = "replicate 2",
                                compartment = tail(clusters, n)))
                              
  minMaxDistances <- lapply(centroids,
                            function(x) {
                              concordanceFunction(x, centroids[[1]], centroids[[2]])
                            }
                           )
  
  concordances <- bind_rows(concordances,
                            tibble(chromosome = chr,
                                   position = rep(positions,
                                                  length(replicate_names)),
                                   replicate = rep(replicate_names,
                                                   each = ncol(bigmat)),
                                   concordance = 
                                     apply(bigmat, 1, function(x) {
                                         (2 * concordanceFunction(x,
                                                                  centroids[[1]],
                                                                  centroids[[2]]) -
                                           minMaxDistances[[1]]) / 
                                           (minMaxDistances[[2]] - minMaxDistances[[1]])
                                     })
                           ))
    
    distances %<>%
                 bind_rows(tibble(chromosome = chr,
                      position = rep(positions, length(replicate_names)),
                      cluster = 1,
                      replicate = rep(replicate_names, each = ncol(bigmat)),
                      distance = apply(bigmat, 1, function(x) {
                                      euclideanDistance(x, centroids[[1]])}))) %>%
                 bind_rows(tibble(chromosome = chr,
                      position = rep(positions, length(replicate_names)),
                      cluster = 2,
                      replicate = rep(replicate_names, each = ncol(bigmat)),
                      distance = apply(bigmat, 1, function(x) {
                                      euclideanDistance(x, centroids[[2]])})))
}

for (i in c(1, 2)) {
  output <- distances %>%
    filter(cluster == i) %>%
    select(-c(cluster)) %>%
    spread(key = replicate, value = distance)
  write_tsv(format(as.data.frame(output), scientific = FALSE),
            path = distance_files[i],
            quote_escape = FALSE,
            col_names = TRUE)
}

output <- compartments %>%
  spread(key = replicate, value = compartment)
write_tsv(format(as.data.frame(output), scientific = FALSE),
          path = argv$compartments,
          quote_escape = FALSE,
          col_names = TRUE)

output <- concordances %>%
  spread(key = replicate, value = concordance)
write_tsv(format(as.data.frame(output), scientific = FALSE),
          path = argv$concordance,
          quote_escape = FALSE,
          col_names = TRUE)