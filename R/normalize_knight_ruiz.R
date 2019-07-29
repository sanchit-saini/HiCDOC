#!/usr/bin/env Rscript

library('argparser')
library('HiCcompare')
library('magrittr')
library('tidyr')
library('readr')

parser <- arg_parser("Normalize with the Knight-Ruiz method")

parser <- add_argument(parser = parser, arg = "--input",
                       help = "Input matrix")
parser <- add_argument(parser = parser, arg = "--output",
                       help = "Output matrix")
argv <- parse_args(parser)

# argv <- parse_args(parser, c("--input",  "/home/mazytnicki/Desktop/Projects/Kurylo/Test/normalized.tsv",
#                              "--output", "/home/mazytnicki/Desktop/Projects/Kurylo/Test/test.tsv"))

# Read input
input <- read_tsv(
  file = argv$input,
  comment = '#',
  col_names = TRUE
)

head(input)

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

step <- min(input$`position 1`[input$`position 1` > 0])
input %<>%
  mutate(bin1 = `position 1` / step + 1) %>%
  mutate(bin2 = `position 2` / step + 1)

outputTidy <- tibble(chromosome = double(),
                     replicate = character(),
                     values = double(),
                     bin1 = integer(),
                     bin2 = integer())
for (chr in unique(input$chromosome)) {
  inputChromosome <- filter(input, chromosome == chr)
  inputReplicate  <- tibble(bin1 = inputChromosome$bin1,
                            bin2 = inputChromosome$bin2,
                            data = 0)
  n <- max(inputChromosome$bin1, inputChromosome$bin2)
  message(paste0("Chromosome ", chr, ", of dim. ", n))
  for (col in 4:(ncol(input) - 2)) {
    message(paste0("  Replicate ", colnames(input)[[col]]))
    inputReplicate$data <- inputChromosome[[col]]
    mat <- matrix(0, nrow = n, ncol = n)
    tmp <- as.matrix(inputReplicate)
    mat[ tmp[, 1:2] ] <- tmp[, 3]
    mat[ rev(tmp[, 1:2]) ] <- tmp[, 3]
    nullRows <- which(colSums(mat) == 0)
    mat[nullRows, nullRows] <- 1
    matKR <- KRnorm(mat)
    matKR[nullRows, ] <- 0
    matKR[, nullRows] <- 0
    vecKR <- as.vector(t(matKR))
    vecKR[is.na(vecKR)] <- 0
    tmpOutput <- tibble(bin1 = rep(seq(n), each = n),
                        bin2 = rep(seq(n), times = n),
                        values = vecKR,
                        chromosome = chr,
                        replicate = colnames(input)[[col]])
    outputTidy %<>% bind_rows(tmpOutput)
  }
}
message("Here 1")
outputTidy %<>%
  filter(bin1 <= bin2) %>%
  filter(rowSums(.[4:(ncol(.)-2)]) != 0.0) %>%
  mutate(`position 1` = (`bin1` - 1) * step,
                       `position 2` = (`bin2` - 1) * step) %>%
  select(-c(`bin1`, `bin2`))

message("Here 2")
message(nrow(outputTidy))
message(ncol(outputTidy))

write_tsv(outputTidy, path="/home/mazytnicki/Desktop/Projects/Kurylo/Test/tmp.tsv")

output <- outputTidy %>%
  spread(replicate, values)

message("Colnames of output")
message(colnames(output))

# Write to output
# writeLines(comments, argv$output)

write_tsv(format(as.data.frame(output), scientific = FALSE),
          path = argv$output,
          quote_escape = FALSE,
          col_names = TRUE)
