#!/usr/bin/env Rscript

library('argparser')
library('plotly')
library('dplyr')
library('tidyr')
library('magrittr')
library('readr')


parser <- arg_parser("Find suitable p-values")

parser <- add_argument(parser = parser, arg = "--compartments",
                       help = "Input compartments file")
parser <- add_argument(parser = parser, arg = "--concordance",
                       help = "Input concordance file")
parser <- add_argument(parser = parser, arg = "--pvalues",
                       help = "Output pvalues in BED file")
parser <- add_argument(parser = parser, arg = "--output",
                       help = "Output plot")
argv <- parse_args(parser)

argv <- parse_args(parser, c("--concordance",  "/home/mazytnicki/Desktop/Projects/Kurylo/Test/concordance.tsv",
                             "--compartments", "/home/mazytnicki/Desktop/Projects/Kurylo/Test/compartments.tsv",
                             "--pvalues",      "/home/mazytnicki/Desktop/Projects/Kurylo/Test/pvalues.bed"))

# Read concordance
concordance <- read_tsv(
  file = argv$concordance,
  comment = '#',
  col_names = TRUE
)

# Compute median of differences
concordance %<>%
  gather("replicate", "values", 3:ncol(concordance)) %>%
  separate("replicate", c(NA, "condition", "replicate"))
    

# compartments <- as_tibble(lapply(read.table(file = argv$compartments,
#                                             sep = '\t',
#                                             comment.char = '#',
#                                             header = TRUE),
#                                  as.integer))
compartments <- read_tsv(file = argv$compartments,
                                comment = '#',
                                col_names = TRUE)

compartments %<>%
  gather("replicate", "values", 3:ncol(compartments)) %>%
  filter(grepl(".*\\.1", replicate)) %>%
  separate(replicate, c(NA, "condition", NA))

# Read p-values
pvalues <- read_tsv(
    file = argv$pvalues,
    comment = '#',
    col_names = FALSE
  )

if (nrow(pvalues) > 0) {
  pvalues %<>%
    rename(chromosome = X1,
           position = X2, 
           pvalue = X6) %>%
    select(-c(X3, X4, X5))
}


# as_tibble(read.table(
#   file = argv$pvalues,
#   sep = '\t',
#   comment.char = '#',
# )) %>%
#   rename(chromosome = V1,
#          position = V2, 
#          pvalue = V6) %>%
#   select(-c(V3, V4, V5)) -> pvalues


# compartments %>% group_by(.dots = c("chromosome", "condition")) %>%
#   mutate(change = c(1, diff(values) != 0)) %>%
#   mutate(change = c(change[-n()], 1)) %>%
#   filter(change == 1) %>%
#   mutate(values_prec = 3 - values) %>%
#   gather(`values_prec`, `values`, key = "tmp", value = "y") %>%
#   filter(position != max(position) | tmp != "values_prec") %>%
#   arrange(position, .by_group = TRUE) %>%
#   ungroup() %>% 
#   mutate(y = y - 1) %>%
#   select(-c(change, tmp)) %>%
#   transform(id = as.integer(factor(condition))) %>%
#   plot_ly(x = ~position,
#           y = ~y,
#           color = ~condition,
#           type = 'scatter',
#           mode = 'lines',
#           fill = 'tozeroy',
#           width = 0,
#           yaxis = ~paste0("y", id)) %>%
#   add_trace(showlegend = FALSE) %>%
#   subplot(nrows = 2, shareX = TRUE) -> pCompartments

concordance %<>% 
  unite("cond.rep", condition, replicate, remove = FALSE) %>% 
  transform(id = as.integer(factor(condition))) %>%
  transform(cond.rep = factor(cond.rep))

tmps <- lapply(c("1", "2"),
             function(xx) {
               concordance %>%
                 filter(condition == xx) %>%
                 plot_ly(x = ~position,
                         y = ~values,
                         color = ~cond.rep,
                         type = 'scatter',
                         mode = 'lines',
                         legendgroup = ~cond.rep,
                         transforms = list(list( 
                              type = 'filter',
                              target = ~chromosome,
                              operation = '=',
                              value = unique(concordance$chromosome)[1]))
                         ) %>% 
                 add_trace(showlegend = FALSE)
             }
)

subplot(tmps, nrows = 2, shareX = TRUE, titleY = TRUE)  %>%
layout(title = 'Hi-C results',
       xaxis = list(title = "position"),
       yaxis = list(title = c("Conc 1", "Conc 2")),
       updatemenus = list(
         list(
           type = 'dropdown',
           active = 0,
           buttons = apply(as.data.frame(unique(compartments$chromosome)), 1, 
                           function(x) list(method = 'restyle', args = list('transforms[0].value', x), label = x)))
         )
       ) -> pConcordance
  
compartmentsNames <- compartments %>%
  mutate("1" = if_else(values == 1, 1, 0), "2" = if_else(values == 2, 1, 0)) %>%
  gather("compartment", "value", c("1", "2")) %>%
  select(-values)

empty_ax <- list(
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)

# compartmentsChanges <- compartments %>%
#   spread(condition, values) %>%
#   mutate("1_to_2" = if_else(`1` < `2`, 1, 0)) %>%
#   mutate("2_to_1" = if_else(`1` > `2`, 1, 0)) %>%
#   select(-c(`1`, `2`)) %>%
#   gather("change", "value", `1_to_2`, `2_to_1`) %>%
#   filter(value > 0)
# compartmentsChanges %>%
#   plot_ly(x = ~position,
#           y = ~value,
#           type = 'bar',
#           color = ~change,
#           legendgroup = ~change,
#           width = 1,
#           transforms = list(list( 
#                type = 'filter',
#                target = ~chromosome,
#                operation = '=',
#                value = unique(compartments$chromosome)[1]))
#           ) %>% 
#   layout(barmode = 'stack',
#          bargap = 0,
#          yaxis = c(empty_ax, list(title="Changes"))) -> pCompartmentsChanges

compartmentsChanges <- pvalues %>%
  mutate("logPV" = sign(pvalue) * (-log10(abs(pvalue) + 0.00000001)))

compartmentsChanges %>%
  plot_ly(x = ~position,
          y = 1,
          type = 'bar',
          color = ~logPV,
          width = 1,
          transforms = list(list( 
               type = 'filter',
               target = ~chromosome,
               operation = '=',
               value = unique(compartments$chromosome)[1]))
          ) %>% 
  layout(barmode = 'stack',
         bargap = 0,
         yaxis = c(empty_ax, list(title="Changes"))) -> pCompartmentsChanges
pCompartmentsChanges

ps <- lapply(c("1", "2"),
             function(xx) {
               compartmentsNames %>%
                 filter(condition == xx) %>%
                 plot_ly(x = ~position,
                         y = ~value,
                         type = 'bar',
                         color = ~compartment,
                         legendgroup = ~compartment, 
                         width = 1,
                         transforms = list(list( 
                              type = 'filter',
                              target = ~chromosome,
                              operation = '=',
                              value = unique(compartments$chromosome)[1]))) %>%
                 layout(
                   barmode = 'stack',
                   bargap = 0,
                   yaxis = c(empty_ax, list(title = paste0("C ", xx))))
               }
             )
pCompartments <- subplot(ps, nrows = 2, shareX = TRUE, titleY = TRUE) %>%
  layout(showlegend = FALSE)

# ps <- lapply(c("1", "2"),
#              function(xx) {
#                 compartments %>%
#                   filter(condition == xx) %>%
#                   group_by(.dots = c("chromosome", "condition")) %>%
#                   mutate(change = c(1, diff(values) != 0)) %>%
#                   mutate(change = c(change[-n()], 1)) %>%
#                   filter(change == 1) %>%
#                   mutate(values_prec = 3 - values) %>%
#                   gather(`values_prec`, `values`, key = "tmp", value = "y") %>%
#                   filter(position != 0 | tmp != "values_prec") %>%
#                   filter(position != max(position) | tmp != "values_prec") %>%
#                   arrange(position, .by_group = TRUE) %>%
#                   ungroup() %>% 
#                   mutate(y = y - 1) %>%
#                   select(-c(change, tmp)) %>%
#                   transform(id = as.integer(factor(condition))) %>%
#                   plot_ly(x = ~position,
#                           y = ~y,
#                           color = ~condition,
#                           type = 'scatter',
#                           mode = 'lines',
#                           fill = 'tozeroy',
#                           width = 0,
#                           yaxis = ~paste0("y", id),
#                           transforms = list(list( 
#                                type = 'filter',
#                                target = ~chromosome,
#                                operation = '=',
#                                value = unique(compartments$chromosome)[1]))) %>%
#                   add_trace(showlegend = FALSE) %>%
#                   layout(yaxis = c(empty_ax, list(title = c("C 1", "C 2"))))
#  }) 
# subplot(ps, nrows = 2, shareX = TRUE, titleY = TRUE) -> pCompartments

subplot(pCompartmentsChanges, pCompartments, nrows = 2, shareX = TRUE, heights = c(0.3, 0.6), titleY = TRUE) -> pHeader

#subplot(pHeader, pConcordance, nrows = 2, shareX = TRUE, titleX = TRUE, titleY = TRUE, heights = c(0.3, 0.6)) -> p
subplot(pHeader, pConcordance, nrows = 2, shareX = TRUE, titleX = TRUE, titleY = TRUE, heights = c(0.3, 0.6)) %>% layout(autosize = F, width = 1000, height = 1000) -> p

htmlwidgets::saveWidget(p, argv$output, knitrOptions = list(out.width = 20, out.height = 20))

