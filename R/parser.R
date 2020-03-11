##- Parse input data ---------------------------------------------------------#
##----------------------------------------------------------------------------#

#' @export
parseInteractionMatrix3Columns <- function(object) {
  object@interactions <- read.table(
    file = object@inputPath,
    sep = '\t',
    header = TRUE,
    comment.char = '#',
    check.names = FALSE
  )
  if (colnames(object@interactions)[[1]] != "chromosome") {
    stop("First column of the input matrix should be named 'chromosome'.")
  }
  if (colnames(object@interactions)[[2]] != "position 1") {
    stop("Second column of the input matrix should be named 'position 1'.")
  }
  if (colnames(object@interactions)[[3]] != "position 2") {
    stop("Third column of the input matrix should be named 'position 2'.")
  }
  conditions.replicates <- colnames(object@interactions)[4:ncol(object@interactions)]
  if (!all(grepl("^replicate (.+?)\\..+$", conditions.replicates))) {
    stop(paste(
      "Fourth to last columns of the input matrix",
      "should be named 'replicate x.y'",
      "with x the condition number",
      "and y the replicate number."
    ))
  }
  object@conditions <- gsub("^replicate (.+?)\\..+$", '\\1', conditions.replicates)
  object@replicates <- gsub("^replicate .+?\\.(.+)$", '\\1', conditions.replicates)
  object@interactions %<>%
    gather(
      conditions.replicates,
      key = condition.replicate,
      value = value
    ) %>%
    separate(
      condition.replicate,
      c(NA, "condition", "replicate"),
      remove = FALSE
    ) %>%
    select(-condition.replicate) %>%
    mutate(
      chromosome = factor(chromosome),
      condition = factor(condition),
      replicate = factor(replicate)
    )
  return(object)
}


##- Parse a single cool Matrix -----------------------------------------------#
##----------------------------------------------------------------------------#

h5readCatch <- function(file, name) {
  return(tryCatch(
    h5read(file = file, name = name),
    error = function(e) {
      stop(paste0(
        "The file '",
        file,
        "' does not seem to have the .cool format: ",
        "path '",
        name,
        "' is missing."
      ), call. = FALSE)
    }
  ))
}

parseCoolMatrix <- function(fileName) {
  bins <- tibble(
    chromosome = factor(h5readCatch(
      file = fileName,
      name = "bins/chrom"
    )),
    start = h5readCatch(
      file = fileName,
      name = "bins/start"
    ),
    end = h5readCatch(
      file = fileName,
      name = "bins/end"
    )
  )
  step <- bins$end - bins$start
  pc <- sum(step == max(step)) / length(step)
  if (length(step) < 0.9) {
    stop(paste0(
      "Cannot parse cool file ",
      fileName,
      ": fixed width only."
    ))
  }
  step <- max(step)

  bins %<>%
    select(-(end)) %>%
    rename(position = start) %>%
    rowid_to_column("id") %>%
    mutate(id = id - 1)

  data <- tibble(
    id1 = h5readCatch(
      file = fileName,
      name = "pixels/bin1_id"
    ),
    id2 = h5readCatch(
      file = fileName,
      name = "pixels/bin2_id"
    ),
    value = h5readCatch(
      file = fileName,
      name = "pixels/count"
    )
  ) %>%
    left_join(bins, by = c("id1" = "id")) %>%
    rename(
      `chromosome 1` = chromosome,
      `position 1` = position
    ) %>%
    select(-id1) %>%
    left_join(bins, by = c("id2" = "id")) %>%
    rename(
      `chromosome 2` = chromosome,
      `position 2` = position
    ) %>%
    select(-id2) %>%
    filter(`chromosome 1` == `chromosome 2`) %>%
    select(-`chromosome 2`) %>%
    rename(chromosome = `chromosome 1`) %>%
    select(chromosome, `position 1`, `position 2`, value)

  return(data)
}

mergeMatrices <- function(object, matrices) {
  for (i in seq_along(matrices)) {
    matrices[[i]]$replicate <- object@replicates[[i]]
    matrices[[i]]$condition <- object@conditions[[i]]
  }
  object@interactions <- bind_rows(matrices) %>%
    mutate(chromosome = factor(chromosome)) %>%
    mutate(replicate = factor(replicate)) %>%
    mutate(condition = factor(condition))

  return(object)
}


##- parseInteractionMatrixCool -----------------------------------------------#
##----------------------------------------------------------------------------#
#' Read interaction matrices in .cool format.
#'
#' This function parses the .cool files provided in
#' \code{object@inputMatrixPath}.
#'
#' @name parseInteractionMatrixCool
#' @rdname parseInteractionMatrixCool
#'
#' @aliases parseInteractionMatrixCool
#'
#' @param object An \code{HiCDOCExp} object.
#'
#' @return object An \code{HiCDOCExp} object.
#'
#' @export
parseInteractionMatrixCool <- function(object) {
  matrices <- lapply(object@inputPath, parseCoolMatrix)
  return(mergeMatrices(object, matrices))
}

parseHicMatrix <- function(fileName, resolution = resolution) {
  message(paste0("Parsing .hic matrix '", fileName, "'."))
  return(parseHic(fileName, resolution))
}

##- parseInteractionMatrixHic ------------------------------------------------#
##----------------------------------------------------------------------------#
#' Read interaction matrices in .hic format.
#'
#' This function parses the .hic files provided in
#' \code{object@inputMatrixPath}.
#'
#' @name parseInteractionMatrixHic
#' @rdname parseInteractionMatrixHic
#'
#' @aliases parseInteractionMatrixHic
#'
#' @param object An \code{HiCDOCExp} object.
#'
#' @return object An \code{HiCDOCExp} object.
#'
#' @export
parseInteractionMatrixHic <- function(object) {
  matrices <- bplapply(object@inputPath, parseHicMatrix, resolution = object@binSize)
  matrices <- map2(matrices, object@replicates, ~ mutate(.x, replicate = .y))
  matrices <- map2(matrices, object@conditions, ~ mutate(.x, condition = .y))
  object@interactions <- bind_rows(matrices)
  #object@interactions <- tibble()
  # for (i in seq_along(object@inputPath)) {
  #     object@interactions <- bind_rows(object@interactions,
  #         parseHicMatrix(object@inputPath[[i]], object@binSize) %>%
  #         mutate(replicate = object@replicates[[i]]) %>%
  #         mutate(condition = object@conditions[[i]]))
  # }
  object@interactions <- object@interactions %>%
    mutate(chromosome = factor(chromosome)) %>%
    mutate(replicate = factor(replicate)) %>%
    mutate(condition = factor(condition)) %>%
    rename(`position 1` = `position.1`) %>%
    rename(`position 2` = `position.2`)
  return(object)
}
