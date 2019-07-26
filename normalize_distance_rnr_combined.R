#!/usr/bin/env Rscript

library('argparser')
library('dplyr')
library('tidyr')
library('magrittr')
library('readr')
library('ggplot2')
library('ggExtra')
library('purrr')

parser <- arg_parser("Normalize using loess")

parser <- add_argument(parser = parser, arg = "--input",
                       help = "Input matrix")
parser <- add_argument(parser = parser, arg = "--output",
                       help = "Output matrix")
parser <- add_argument(parser = parser, arg = "--span",
                       help = "Loess span", default = 0.2)
parser <- add_argument(parser = parser, arg = "--before",
                       help = "Output figure before normalization")
parser <- add_argument(parser = parser, arg = "--after",
                       help = "Output figure after normalization")
argv <- parse_args(parser)

argv <- parse_args(parser, c("--input",  "/home/mazytnicki/Desktop/Projects/Kurylo/data/500k/matrices/500k_raw.tsv",
                             "--output", "/home/mazytnicki/Desktop/Projects/Kurylo/Test/test.tsv",
                             "--before", "/home/mazytnicki/Desktop/Projects/Kurylo/Test/test_before.png",
                             "--after",  "/home/mazytnicki/Desktop/Projects/Kurylo/Test/test_after.png"))


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

head(input)

step <- min(input$`position 1`[input$`position 1` > 0])
input_tidy <- input %>%
  gather(4:ncol(.), key = "replicate", value = "values") %>%
  mutate(distance = `position 2` - `position 1`)

head(input_tidy)

normalize_mean <- function(input_tidy) {
  input_mean <- input_tidy %>%
    select(c(distance, values)) %>%
    group_by(distance) %>%
    summarize(mean_value = mean(values))
  
  head(input_mean)

  l <- loess(mean_value ~ distance,
             data = input_mean,
             span = argv$span)
  input_mean %<>%
    mutate(loess = predict(l))
  
  p <- ggplot() + 
    stat_bin_hex(aes(x = input_tidy$distance, y = input_tidy$values)) +
    geom_line(aes(x = input_mean$distance,
                  y = input_mean$mean_value,
                  color = "mean")) +
    geom_line(aes(x = input_mean$distance,
                  y = input_mean$loess,
                  color = "loess")) +
    scale_colour_manual("", 
                        breaks = c("mean", "loess"),
                        values = c("black", "red")) +
    scale_y_log10()
  
  p <- ggMarginal(p, margins = "x", type = "histogram", fill = "transparent")
  ggsave(argv$before, p)
  
  head(input_mean)
  input_mean %<>% select(-mean_value)
  
  output <- 
    input_tidy %>%
    left_join(input_mean, by = "distance")
  
  return(output)
}

normalize_sample <- function(input_tidy) {

  input_sampled <- input_tidy %>%
    select(c(distance, values)) %>%
    rename(sampled_distance = distance) %>%
    sample_n(size = 100000) %>%
    arrange(sampled_distance)
  
  l <- loess(values ~ sampled_distance,
             data = input_sampled,
             span = argv$span)
  input_sampled %<>%
    mutate(loess = predict(l))
  
  p <- ggplot(input_sampled,
              aes(x = sampled_distance)) +
    stat_bin_hex(aes(y = values)) +
    geom_line(aes(y = loess)) +
    scale_x_log10()
    scale_y_log10()
  p <- ggMarginal(p, margins = "x", type = "histogram", fill = "transparent")
  p
  ggsave(argv$before, p)

  input_sampled %<>%
    select(-values) %>%
    unique()
  
  head(input_sampled)
  tail(input_sampled)
  
  sampled_distances <- unique(sort(input_sampled$sampled_distance))
  unique_distances <- unique(sort(input_tidy$distance))
  distance_map <- tibble(distance = unique_distances,
                         sampled_distance =
                           sapply(unique_distances, function(x) {
                             sampled_distances[which.min(abs(x - sampled_distances))]
                           }))
  head(distance_map)
  tail(distance_map)
  
  output <- input_tidy %>%
    left_join(distance_map, by = "distance") %>%
    left_join(input_sampled, by = "sampled_distance") %>%
    select(-sampled_distance)
  
  return(output)
}

output <- normalize_sample(input_tidy)
#output <- normalize_mean(input_tidy)

output %<>%
  mutate(ratio = log2(values / loess))

p <- ggplot(output,
            aes(x = distance, y = ratio)) +
  stat_bin_hex() +
  geom_smooth() +
  scale_x_log10()
  scale_y_log10()
p <- ggMarginal(p, margins = "x", type = "histogram", fill = "transparent")
p

output %<>%
  select(-c(values, distance, loess)) %>%
  spread(replicate, ratio)
  
# Write to output
#writeLines(comments, argv$output)

write_tsv(format(as.data.frame(output), scientific = FALSE),
          path = argv$output,
          quote_escape = FALSE,
          col_names = TRUE)
