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
#' @param object An \code{HiCDOCDataSet} object.
#' 
#' @return object An \code{HiCDOCDataSet} object.
#' @examples
#' \dontrun{
#' linkToMatrix <-system.file("extdata", "sampleMatrix.tsv", package = "HiCDOC")
#' data   <- makeHiCDOCDataSet(inputPath = linkToMatrix)
#' object <- parseInteractionMatrix3Columns(data)
#' }
parseInteractionMatrix3Columns <- function(object) {
    object@interactions <- utils::read.table(
        file = object@inputPath,
        sep = '\t',
        header = TRUE,
        comment.char = '#',
        check.names = FALSE
    ) %>% dplyr::as_tibble()
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
            replicate = factor(replicate),
            value = as.numeric(value)
        )
    if(is.integer(object@interactions$value)) 
        object@interactions$value <- as.numeric(object@interactions$value)
    return(object)
}

##- parseCoolMatrix ----------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Parse a single interaction matrix in .(m)cool format.
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
  bins <- dplyr::tibble(
    chromosome = factor(rhdf5::h5read(
      file = filePath,
      name = uri("bins/chrom")
    )),
    start = rhdf5::h5read(
      file = filePath,
      name = uri("bins/start")
    ),
    end = rhdf5::h5read(
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
      dplyr::mutate(id = seq_len(nrow(bins)) - 1) %>%
      dplyr::select(id, chromosome, position)
  rownames(bins) <- NULL

  data <- dplyr::tibble(
    id1 = rhdf5::h5read(
      file = filePath,
      name = uri("pixels/bin1_id")
    ),
    id2 = rhdf5::h5read(
      file = filePath,
      name = uri("pixels/bin2_id")
    ),
    value = rhdf5::h5read(
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
  
    if(is.integer(data$value)) 
      data <- data %>% dplyr::mutate(value=as.numeric(value))

    return(data)
}

##- parseMCoolMatrix ---------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Parse a single interaction matrix in .mcool format.
#'
#' @name parseMCoolMatrix
#' @rdname parseMCoolMatrix
#'
#' @aliases parseMCoolMatrix
#'
#' @param fileName The name of the matrix file in HDF5 format.
#' @param resolution The chosen resolution (in bp)
#' @return object An \code{tibble} storing the (sparse) interaction matrix.
parseMCoolMatrix <- function(fileName, resolution = resolution) {
    return(parseCoolMatrix(paste0(fileName, "::resolutions/", resolution)))
}


##- mergeMatrices ------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Merge the matrices which have been parsed separately.
#' Store the result into a \code{HiCDOCDataSet} object.
#'
#' @name mergeMatrices
#' @rdname mergeMatrices
#'
#' @aliases mergeMatrices
#'
#' @param object A \code{HiCDOCDataSet} object, where the matrices should be
#'          stored.
#' @param matrices The different matrices (one per sample).
#' @return object The \code{HiCDOCDataSet} object, where the matrices have been
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
#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return object An \code{HiCDOCDataSet} object.
parseInteractionMatrixCool <- function(object) {
    matrices <- pbapply::pblapply(object@inputPath, parseCoolMatrix)
    return(mergeMatrices(object, matrices))
}

##- parseInteractionMatrixMCool ----------------------------------------------#
##----------------------------------------------------------------------------#
#' Read interaction matrices in .mcool format.
#'
#' This function parses the .mcool files provided in
#' \code{object@inputMatrixPath} at a given resolution.
#'
#' @name parseInteractionMatrixMCool
#' @rdname parseInteractionMatrixMCool
#'
#' @aliases parseInteractionMatrixMCool
#'
#' @param object     An \code{HiCDOCDataSet} object.
#' @param resolution The chosen resolution, in base pairs.
#' @return object    An \code{HiCDOCDataSet} object.
parseInteractionMatrixMCool <- function(object, resolution) {
    matrices <- pbapply::pblapply(object@inputPath,
                       parseMCoolMatrix,
                       resolution = resolution)
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
#' \code{object@inputPath}.
#'
#' @name parseInteractionMatrixHic
#' @rdname parseInteractionMatrixHic
#'
#' @aliases parseInteractionMatrixHic
#'
#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return object An \code{HiCDOCDataSet} object.
parseInteractionMatrixHic <- function(object) {
    matrices <-
      pbapply::pblapply(object@inputPath, 
                               parseHicMatrix, 
                               resolution = object@binSize)
    matrices <-
        purrr::map2(matrices, object@replicates, 
                    ~ dplyr::mutate(.x, replicate = .y))
    matrices <-
        purrr::map2(matrices, object@conditions, 
                    ~ dplyr::mutate(.x, condition = .y))
    object@interactions <- dplyr::bind_rows(matrices) %>% dplyr::as_tibble()
    object@interactions %<>%
        dplyr::mutate(chromosome = factor(chromosome, 
            levels = 
              gtools::mixedsort(unique(object@interactions$chromosome)))) %>%
        dplyr::mutate(replicate = factor(replicate)) %>%
        dplyr::mutate(condition = factor(condition))
    if(is.integer(object@interactions$value))
        object@interactions$value <- as.numeric(object@interactions$value)
    return(object)
}


##- parseHicPro --------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Parse a single interaction matrix in HiC-Pro format
#'
#' @name parseHicPro
#' @rdname parseHicPro
#'
#' @aliases parseHicPro
#'
#' @param vectFiles 2 length vector, of the links to .matrix and .bed files.
#' @return object An \code{tibble} storing the (sparse) interaction matrix.
parseHicPro <- function(vectFiles) {
    if (length(vectFiles) != 2) {
      stop("vectFiles must be of length 2")
    }
    matrixFile <- vectFiles[1]
    bedFile <- vectFiles[2]
    matrixDf <- utils::read.table(
        matrixFile, 
        header = FALSE, 
        stringsAsFactors = FALSE,
        col.names = c("startIndex", "stopIndex", "value")
    )
    bedDf <- utils::read.table(
        bedFile, 
        header = FALSE, 
        stringsAsFactors = FALSE,
        col.names = c("chromosome", "start", "end", "index")
    )
    
    matrixDf <- dplyr::as_tibble(matrixDf)
    bedDf <- dplyr::as_tibble(bedDf)
    
    # Resolution
    tabDif <- sort(table(abs(bedDf$end-bedDf$start)), decreasing = TRUE)
    resolution <- as.numeric(names(tabDif[1]))
    
    positions <- bedDf %>%
        dplyr::arrange(chromosome, index) %>%
        dplyr::group_by(chromosome) %>%
        dplyr::mutate(bin = index - min(index) + 1) %>%
        dplyr::mutate(end = ifelse(end == lead(start), end - 1, end)) %>%
        dplyr::ungroup() %>%
        dplyr::select(chromosome, start, end, bin)
    
    # Merge data
    matrixDf <- dplyr::left_join(
        matrixDf, 
        bedDf %>% dplyr::select(chr1 = chromosome, 
                                startIndex = index, 
                                position.1 = start),
        by="startIndex"
    )
    matrixDf %<>% dplyr::select(-startIndex)
    matrixDf <- dplyr::left_join(
        matrixDf, 
        bedDf %>% dplyr::select(chr2 = chromosome, 
                                stopIndex = index, 
                                position.2 = start),
        by="stopIndex")
    message("Removing the inter-chromosomes interactions")
    matrixDf %<>% dplyr::filter(chr1 == chr2) %>%
       dplyr::select(chromosome = chr1,
                     position.1, 
                     position.2, 
                     value)
    if(is.integer(matrixDf$value))
        matrixDf$value <- as.numeric(matrixDf$value)
    return(list("matrix" = matrixDf, 
                "resolution" = resolution, 
                "positions" = positions))
}


##- parseInteractionMatrixHicPro ---------------------------------------------#
##----------------------------------------------------------------------------#
#' Read interaction matrices in Hic-Pro format.
#'
#' This function parses the HiC-Pro files provided in
#' \code{object@inputPath}.
#'
#' @name parseInteractionMatrixHicPro
#' @rdname parseInteractionMatrixHicPro
#'
#' @aliases parseInteractionMatrixHicPro
#'
#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return object An \code{HiCDOCDataSet} object.
parseInteractionMatrixHicPro <- function(object) {
    matrices <-pbapply::pblapply(object@inputPath, 
                                      function(x) parseHicPro(x))
    matrices <-
      purrr::map2(matrices, 
                  object@replicates, 
                  ~ dplyr::mutate(.x, replicate = .y, .after = position.2))
    matrices <-
      purrr::map2(matrices, 
                  object@conditions, 
                  ~ dplyr::mutate(.x, condition = .y, .after = position.2))
    object@interactions <- dplyr::bind_rows(matrices) %>% dplyr::as_tibble()
    object@interactions %<>%
      dplyr::mutate(chromosome = factor(chromosome, 
          levels = 
            gtools::mixedsort(unique(object@interactions$chromosome)))) %>%
      dplyr::mutate(replicate = factor(replicate)) %>%
      dplyr::mutate(condition = factor(condition))
    return(object)
}
