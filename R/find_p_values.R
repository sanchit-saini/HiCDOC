#!/usr/bin/env Rscript

library('argparser')
library('ggplot2')
library('tidyr')
library('dplyr')
library('readr')


parser <- arg_parser("Find suitable p-values")

parser <- add_argument(parser = parser, arg = "--compartments",
                       help = "Input compartments file")
parser <- add_argument(parser = parser, arg = "--concordance",
                       help = "Input concordance file")
parser <- add_argument(parser = parser, arg = "--pvalues",
                       help = "Output pvalues in BED file")
parser <- add_argument(parser = parser, arg = "--distribution",
                       help = "Output plot")
argv <- parse_args(parser)

#argv <- parse_args(parser, c("--concordance",  "/home/mazytnicki/Desktop/Projects/Kurylo/Test/concordance.tsv",
#                             "--compartments", "/home/mazytnicki/Desktop/Projects/Kurylo/Test/compartments.tsv"))

# Read concordance
concordance <- read_tsv(
  file = argv$concordance,
  comment = '#',
  col_names = TRUE
)

# Compute median of differences
median_conc <- concordance %>%
  gather("replicate", "values", 3:ncol(concordance)) %>%
  separate("replicate", c(NA, "condition", "replicate")) %>%
  group_by(.dots = c("chromosome", "position", "condition")) %>%
  summarize(values = median(values))
group1 <- median_conc %>% filter(condition == "1")
group2 <- median_conc %>% filter(condition == "2")
diff <- tibble(chromosome = group1$chromosome,
                   position   = group1$position,
                   values     = group2$values - group1$values)
sortedConcordance <- sort(diff$values)

# Plot the distribution
p <- ggplot(diff, aes(x = values)) + geom_histogram() +
     ggtitle("Distribution of the difference of concordances") +
     xlab("concordane")
ggsave(argv$distribution, p)

# Read compartments
compartments <- read_tsv(file = argv$compartments,
                         comment = '#',
                         col_names = TRUE)

# Compute number of compartments
step <- compartments$position[2] - compartments$position[1]
pvalues <- compartments %>%
  gather("replicate", "values", 3:ncol(compartments)) %>%
  group_by(.dots = c("chromosome", "position")) %>%
  separate("replicate", c(NA, "condition", "replicate")) %>%
  summarize(nCompartments = n_distinct(values),
            compartment = first(values)) %>%
  ungroup() %>%
  left_join(diff, by = c("chromosome", "position")) %>%
  rename(start = position) %>%
  mutate(end = start + step) %>%
  mutate(name = paste0("bin_", row_number())) %>%
  mutate(percentRank = percent_rank(values)) %>%
  filter(nCompartments == 2) %>%
  mutate(pvalue = if_else(percentRank < 0.5,
                         2 * percentRank,
                         2 * (1 - percentRank))) %>%
  mutate(pvalue = if_else(pvalue < 0, 0, pvalue)) %>%
  mutate(pvalue = if_else(pvalue > 1, 1, pvalue)) %>%
  mutate(strand = ".") %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  mutate(padj = if_else(compartment == 1, -padj, padj)) %>%
  select(-c(nCompartments, compartment, values, percentRank, pvalue))

# Write bed file
write.table(format(as.data.frame(pvalues), scientific = FALSE),
            file = argv$pvalues,
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F)

