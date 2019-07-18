#!/usr/bin/env Rscript

library('argparser')
library('ggplot2')
library('tidyr')
library('dplyr')


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

# Read concordance
concordance <- as_tibble(lapply(read.table(
  file = argv$concordance,
  sep = '\t',
  comment.char = '#',
  header = TRUE
), as.double))

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
compartments <- as_tibble(lapply(read.table(file = argv$compartments,
                                            sep = '\t',
                                            comment.char = '#',
                                            header = TRUE),
                                 as.integer))

# Compute number of compartments
step <- compartments$position[2] - compartments$position[1]
pvalues <- compartments %>%
  gather("replicate", "values", 3:ncol(compartments)) %>%
  group_by(.dots = c("chromosome", "position")) %>%
  summarize(nCompartments = n_distinct(values)) %>%
  ungroup() %>%
  left_join(diff, by = c("chromosome", "position")) %>%
  rename(start = position) %>%
  mutate(end = start + step) %>%
  mutate(name = paste0("bin_", row_number())) %>%
  mutate(percentRank = percent_rank(values)) %>%
  filter(nCompartments == 2) %>%
  mutate(pvalue = ifelse(percentRank < 0.5,
                         2 * percentRank,
                         2 * (1 - percentRank))) %>%
  mutate(pvalue = ifelse(pvalue < 0, 0, pvalue)) %>%
  mutate(pvalue = ifelse(pvalue > 1, 1, pvalue)) %>%
  mutate(strand = ".") %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  select(-c(nCompartments, values, percentRank, pvalue))

# Write bed file
write.table(format(as.data.frame(pvalues), scientific = FALSE),
            file = argv$pvalues,
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F)

