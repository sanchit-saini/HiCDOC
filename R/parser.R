##- Parse input data ---------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Read interaction matrices in tabular 3-column format.
#'
#' This function parses the tsv files provided in
#' \code{object@inputMatrixPath}.
#'
#' @name parseInteractionMatrix3Columns
#' @rdname parseInteractionMatrix3Columns
#'
#' @aliases parseInteractionMatrix3Columns
#'
#' @param object An \code{HiCDOCExp} object.
#'
#' @return object An \code{HiCDOCExp} object.
#'
#' @examples
#' linkToMatrix <-system.file("extdata", "sampleMatrix.tsv", package = "HiCDOC")
#' data   <- makeHiCDOCDataSet(inputPath = linkToMatrix)
#' object <- parseInteractionMatrix3Columns(data)
#' @export
parseInteractionMatrix3Columns <- function(object) {
    object@interactions <- read.table(
        file = object@inputPath,
        sep = '\t',
        header = TRUE,
        comment.char = '#',
        check.names = FALSE
    ) %>% as_tibble()
    if (colnames(object@interactions)[[1]] != "chromosome") {
        stop("First column of the input matrix should be named 'chromosome'.")
    }
    if (colnames(object@interactions)[[2]] != "position 1") {
        stop("Second column of the input matrix should be named 'position 1'.")
    }
    if (colnames(object@interactions)[[3]] != "position 2") {
        stop("Third column of the input matrix should be named 'position 2'.")
    }
    conditions.replicates <-
        colnames(object@interactions)[4:ncol(object@interactions)]
    if (!all(grepl("^replicate (.+?)\\..+$", conditions.replicates))) {
        stop(
            paste(
                "Fourth to last columns of the input matrix",
                "should be named 'replicate x.y'",
                "with x the condition number",
                "and y the replicate number."
            )
        )
    }
    object@conditions <-
        gsub("^replicate (.+?)\\..+$", '\\1', conditions.replicates)
    object@replicates <-
        gsub("^replicate .+?\\.(.+)$", '\\1', conditions.replicates)
    object@interactions %<>%
        tidyr::gather(conditions.replicates,
                     key = condition.replicate,
                     value = value) %>%
        tidyr::separate(condition.replicate,
                         c(NA, "condition", "replicate"),
                         remove = FALSE) %>%
        dplyr::select(-condition.replicate) %>%
        dplyr::rename(position.1 = `position 1`) %>%
        dplyr::rename(position.2 = `position 2`) %>%
        dplyr::mutate(
            chromosome = factor(chromosome),
            condition = factor(condition),
            replicate = factor(replicate)
        )
    return(object)
}


##- h5readCatch --------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Read an HDF5 file, and extract a field.
#' It raises an exception if the field is not present.
#'
#' @name h5readCatch
#' @rdname h5readCatch
#'
#' @aliases h5readCatch
#'
#' @param file The name of a file in HDF5 format.
#' @param name The name a field name.
#'
#' @return The value associated with the field.
h5readCatch <- function(file, name) {
    return(tryCatch(
        h5read(file = file, name = name),
        error = function(e) {
            stop(
                paste0(
                    "The file '",
                    file,
                    "' does not seem to have the .cool format: ",
                    "path '",
                    name,
                    "' is missing."
                ),
                call. = FALSE
            )
        }
    ))
}

##- parseCoolMatrix ----------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Parse a single interaction matrix in .cool format.
#'
#' @name parseCoolMatrix
#' @rdname parseCoolMatrix
#'
#' @aliases parseCoolMatrix
#'
#' @param fileName The name of the matrix file in HDF5 format.
#' @return object An \code{tibble} storing the (sparse) interaction matrix.
parseCoolMatrix <- function(fileName) {
  splitName <- strsplit(fileName, '::', fixed = TRUE)[[1]]
  # Separate file path from URI in case of mcool file
  filePath <- splitName[1]
  fileUri <- ifelse(length(splitName) > 1, splitName[2], '')
  uri <- function(path) { return(paste(fileUri, path, sep='/')) }
  bins <- tibble(
    chromosome = factor(h5readCatch(
      file = filePath,
      name = uri("bins/chrom")
    )),
    start = h5readCatch(
      file = filePath,
      name = uri("bins/start")
    ),
    end = h5readCatch(
      file = filePath,
      name = uri("bins/end")
    )
  )
  step <- bins$end - bins$start
  pc <- sum(step == max(step)) / length(step)
  if (length(step) < 0.9) {
      stop(paste0("Cannot parse cool file ",
                  fileName,
                  ": fixed width only."))
  }
  step <- max(step)

  bins %<>%
      dplyr::select(-(end)) %>%
      dplyr::rename(position = start) %>%
      tibble::rowid_to_column("id") %>%
      dplyr::mutate(id = id - 1)

  data <- tibble(
    id1 = h5readCatch(
      file = filePath,
      name = uri("pixels/bin1_id")
    ),
    id2 = h5readCatch(
      file = filePath,
      name = uri("pixels/bin2_id")
    ),
    value = h5readCatch(
      file = filePath,
      name = uri("pixels/count")
    )
  ) %>%
    dplyr::left_join(bins, by = c("id1" = "id")) %>%
    dplyr::rename(
      chromosome.1 = chromosome,
      position.1 = position
    ) %>%
    dplyr::select(-id1) %>%
    dplyr::left_join(bins, by = c("id2" = "id")) %>%
    dplyr::rename(
      chromosome.2 = chromosome,
      position.2 = position
    ) %>%
    dplyr::select(-id2) %>%
    dplyr::filter(chromosome.1 == chromosome.2) %>%
    dplyr::select(-chromosome.2) %>%
    dplyr::rename(chromosome = chromosome.1) %>%
    dplyr::select(chromosome, position.1, position.2, value)

    return(data)
}

##- mergeMatrices ------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Merge the matrices which have been parsed separately.
#' Store the result into a \code{HiCDOCExp} object.
#'
#' @name mergeMatrices
#' @rdname mergeMatrices
#'
#' @aliases mergeMatrices
#'
#' @param object A \code{HiCDOCExp} object, where the matrices should be
#'          stored.
#' @param matrices The different matrices (one per sample).
#' @return object The \code{HiCDOCExp} object, where the matrices have been
#'           added.
mergeMatrices <- function(object, matrices) {
    for (i in seq_along(matrices)) {
        matrices[[i]]$replicate <- object@replicates[[i]]
        matrices[[i]]$condition <- object@conditions[[i]]
    }
    object@interactions <- dplyr::bind_rows(matrices) %>%
        dplyr::mutate(chromosome = factor(chromosome)) %>%
        dplyr::mutate(replicate = factor(replicate)) %>%
        dplyr::mutate(condition = factor(condition))

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

##- parseHicMatrix -----------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Parse a single interaction matrix in .hic format.
#'
#' @name parseHicMatrix
#' @rdname parseHicMatrix
#'
#' @aliases parseHicMatrix
#'
#' @param fileName   The file name of the .hic file.
#' @param resolution The chosen resolution (should be present in the
#'                     .hic file).
#'
#' @return object An \code{tibble} storing the (sparse) interaction matrix.
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
    matrices <-
        bplapply(object@inputPath, parseHicMatrix, resolution = object@binSize)
    matrices <-
        purrr::map2(matrices, object@replicates, ~ dplyr::mutate(.x, replicate = .y))
    matrices <-
        purrr::map2(matrices, object@conditions, ~ dplyr::mutate(.x, condition = .y))
    object@interactions <- dplyr::bind_rows(matrices) %>% as_tibble()
    #object@interactions <- tibble()
    # for (i in seq_along(object@inputPath)) {
    #         object@interactions <- bind_rows(object@interactions,
    #                 parseHicMatrix(object@inputPath[[i]], object@binSize) %>%
    #                 mutate(replicate = object@replicates[[i]]) %>%
    #                 mutate(condition = object@conditions[[i]]))
    # }
    object@interactions %<>%
        dplyr::mutate(chromosome = factor(chromosome)) %>%
        dplyr::mutate(replicate = factor(replicate)) %>%
        dplyr::mutate(condition = factor(condition))
    return(object)
}
