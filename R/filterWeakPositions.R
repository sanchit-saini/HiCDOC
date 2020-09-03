#' @export
filterWeakPositions <- function(object) {

  weakBins <- fullInteractions(object) %>%
    mutate(bin = position.1 / object@binSize + 1) %>%
    group_by(chromosome, condition, replicate, bin) %>%
    summarize(mean = mean(value)) %>%
    filter(mean <= object@filterThreshold)

  if (nrow(weakBins) > 0) {

    weakBins %<>% group_by(chromosome, bin) %>% summarize()

    object@weakBins %<>% modifyList(weakBins %>%
      group_split() %>%
      setNames(
        group_keys(weakBins) %>% pull(chromosome)
      ) %>%
      map(function(x) unique(sort(pull(x, bin))))
    )

    notempty <- vapply(object@weakBins, function(x) length(x)>0, FUN.VALUE=TRUE)
    weakBins <- as_tibble(do.call(rbind, mapply(
      function(chromosome, bin) cbind(chromosome, bin),
      names(object@weakBins[notempty]),
      object@weakBins[notempty],
      SIMPLIFY = FALSE
    ))) %>% mutate(
      chromosome = factor(chromosome, levels = object@chromosomes),
      bin = as.double(bin)
    )

    object@interactions %<>%
      mutate(bin = position.1 / object@binSize + 1) %>%
      anti_join(
        weakBins,
        by = c("chromosome", "bin")
      ) %>%
      select(-bin)

    object@interactions %<>%
      mutate(bin = position.2 / object@binSize + 1) %>%
      anti_join(
        weakBins,
        by = c("chromosome", "bin")
      ) %>%
      select(-bin)
  }

  message(
    "Removed ",
    nrow(weakBins),
    " position",
    if (nrow(weakBins) != 1) "s"
  )

  return (object)
}
