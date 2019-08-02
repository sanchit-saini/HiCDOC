normalizeMean <- function(object) {

  input_mean <- object@interactionMatrix %>%
    select(c(distance, values)) %>%
    group_by(distance) %>%
    summarize(mean_value = mean(values))

  head(input_mean)

  l <- loess(mean_value ~ distance,
             data = input_mean,
             span = object@loessSpan)
  input_mean %<>%
    mutate(loess = predict(l)) %>%
    mutate(loess = pmax(loess, 0))

  object@distancePlotRaw <- ggMarginal(p, margins = "x", type = "histogram", fill = "transparent")

  head(input_mean)
  input_mean %<>% select(-mean_value)

  output <-
    object@interactionMatrix %>%
    left_join(input_mean, by = "distance")

  return(output)
}

normalizeSample <- function(object) {

  inputSampled <- object@interactionMatrix %>%
    select(-c(chromosome, replicate, `position 1`, `position 2`)) %>%
    sample_n(size = min(object@sampleSize, nrow(object@interactionMatrix))) %>%
    rename(sampledDistance = distance) %>%
    select(c(sampledDistance, value)) %>%
    arrange(sampledDistance)

  optimizeSpan <- function(model,
                          criterion = c("aicc", "gcv"),
                          spans = c(0.01, 0.9)) {
      getLoessCriterion <- function(x) {
          span <- x$pars$span
          traceL <- x$trace.hat
          sigma2 <- sum(x$residuals^2)/(x$n - 1)
          aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL - 2)
          gcv <- x$n * sigma2/(x$n - traceL)^2
          result <- list(span = span, aicc = aicc, gcv = gcv)
          return(result)
      }
      criterion <- match.arg(criterion)
      f <- function(span) {
          mod <- update(model, span = span)
          getLoessCriterion(mod)[[criterion]]
      }
      result <- optimize(f, spans)
      return(list(span = result$minimum, criterion = result$objective))
  }
  l <- loess(value ~ sampledDistance, data = inputSampled)
  span <- optimizeSpan(l, criterion = "gcv")$span
  l <- loess(value ~ sampledDistance, span = span, data = inputSampled)

  inputSampled %<>%
    mutate(loess = predict(l)) %>%
    mutate(loess = pmax(loess, 0)) %>%
    select(-value) %>%
    unique()

  sampledDistances <- unique(sort(inputSampled$sampledDistance))
  uniqueDistances <- unique(sort(object@interactionMatrix$distance))
  valueMap <- tibble(distance = uniqueDistances,
                         sampledDistance =
                           sapply(uniqueDistances, function(x) {
                             sampledDistances[which.min(abs(x - sampledDistances))]
                           })) %>%
    left_join(inputSampled, by = "sampledDistance") %>%
    select(-sampledDistance)

  object@interactionMatrix %<>%
      left_join(valueMap, by = "distance")

  return(object)
}

#output <- normalizeMean(input_tidy)

#' @export
normalizeDistanceCombined <- function(object) {

  object@interactionMatrix %<>%
        mutate(distance = `position 2` - `position 1`)

  object <- normalizeSample(object)

  object@interactionMatrix %<>%
    mutate(value = log2((value + 0.0001) / (loess + 0.0001))) %>%
    select(-c(distance, loess))

  return(object)
}
