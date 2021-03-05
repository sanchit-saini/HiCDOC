## - Parse input data ---------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Read interaction matrices in tabular 3-column format.
#'
#' This function parses the tsv files provided in
#' \code{object@input}.
#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return object An \code{HiCDOCDataSet} object.
#' @keywords internal
#' @noRd
.parseTabular <- function(object) {
    object@interactions <-
        utils::read.table(
            file = object@input,
            sep = "\t",
            header = TRUE,
            comment.char = "#",
            check.names = FALSE
        ) %>%
        dplyr::as_tibble()

    if (colnames(object@interactions)[[1]] != "chromosome") {
        stop(
            "First column of the input matrix should be named 'chromosome'.",
            call. = FALSE
        )
    }
    if (colnames(object@interactions)[[2]] != "position 1") {
        stop(
            "Second column of the input matrix should be named 'position 1'.",
            call. = FALSE
        )
    }
    if (colnames(object@interactions)[[3]] != "position 2") {
        stop(
            "Third column of the input matrix should be named 'position 2'.",
            call. = FALSE
        )
    }

    conditions.replicates <-
        colnames(object@interactions)[4:ncol(object@interactions)]

    if (!all(grepl("^.+?\\..+$", conditions.replicates))) {
        stop(
            "Fourth to last columns of the input matrix ",
            "should be named 'x.y' ",
            "with x the condition number ",
            "and y the replicate number.",
            call. = FALSE
        )
    }

    object@conditions <-
        gsub("^(.+?)\\..+$", "\\1", conditions.replicates)
    object@replicates <-
        gsub("^.+?\\.(.+)$", "\\1", conditions.replicates)

    object@interactions %<>%
        tidyr::gather(
            conditions.replicates,
            key = condition.replicate,
            value = interaction
        ) %>%
        tidyr::separate(
            condition.replicate,
            c("condition", "replicate"),
            sep = "\\.",
            remove = FALSE
        ) %>%
        dplyr::select(-condition.replicate) %>%
        dplyr::rename(position.1 = `position 1`) %>%
        dplyr::rename(position.2 = `position 2`) %>%
        dplyr::mutate(
            chromosome = factor(chromosome),
            condition = factor(condition),
            replicate = factor(replicate),
            interaction = as.numeric(interaction)
        )

    return(object)
}

## - .parseOneCool ----------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Parse a single interaction matrix in .(m)cool format.
#' @param fileName The name of the matrix file in HDF5 format.
#' @return object An \code{tibble} storing the (sparse) interaction matrix.
#' @keywords internal
#' @noRd
.parseOneCool <- function(fileName) {
    splitName <- strsplit(fileName, "::", fixed = TRUE)[[1]]
    # Separate file path from URI in case of mcool file
    filePath <- splitName[1]
    fileURI <- ifelse(length(splitName) > 1, splitName[2], "")
    uri <- function(path) return(paste(fileURI, path, sep = "/"))
    bins <-
        dplyr::tibble(
            chromosome = factor(
                rhdf5::h5read(file = filePath, name = uri("bins/chrom"))
            ),
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
    if (length(step) < 0.9) {
        stop(
            "Cannot parse cool file '",
            fileName,
            "': fixed width only."
        )
    }
    step <- max(step)

    bins %<>%
        dplyr::select(-(end)) %>%
        dplyr::rename(position = start) %>%
        dplyr::mutate(ide = seq_len(nrow(bins)) - 1) %>%
        dplyr::select(ide, chromosome, position)

    rownames(bins) <- NULL

    data <-
        dplyr::tibble(
            id1 = rhdf5::h5read(
                file = filePath,
                name = uri("pixels/bin1_id")
            ),
            id2 = rhdf5::h5read(
                file = filePath,
                name = uri("pixels/bin2_id")
            ),
            interaction = rhdf5::h5read(
                file = filePath,
                name = uri("pixels/count")
            )
        ) %>%
        dplyr::left_join(bins, by = c("id1" = "ide")) %>%
        dplyr::rename(
            chromosome.1 = chromosome,
            position.1 = position
        ) %>%
        dplyr::select(-id1) %>%
        dplyr::left_join(bins, by = c("id2" = "ide")) %>%
        dplyr::rename(
            chromosome.2 = chromosome,
            position.2 = position
        ) %>%
        dplyr::select(-id2) %>%
        dplyr::filter(chromosome.1 == chromosome.2) %>%
        dplyr::select(-chromosome.2) %>%
        dplyr::rename(chromosome = chromosome.1) %>%
        dplyr::select(chromosome, position.1, position.2, interaction)

    return(data)
}

## - .parseOneMCool ---------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Parse a single interaction matrix in .mcool format.
#' @param fileName The name of the matrix file in HDF5 format.
#' @param resolution The chosen resolution (in bp)
#' @return object An \code{tibble} storing the (sparse) interaction matrix.
#' @keywords internal
#' @noRd
.parseOneMCool <- function(fileName, resolution = resolution) {
    return(.parseOneCool(paste0(fileName, "::resolutions/", resolution)))
}


## - .mergeMatrices ------------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Merge the matrices which have been parsed separately, and put the result
#' in the interactions slot of the object.
#' @param object A \code{HiCDOCDataSet} object, where the matrices should be
#'          stored.
#' @param matrices The different matrices (one per sample).
#' @return object The \code{HiCDOCDataSet} object, where the matrices have been
#'           added.
#' @keywords internal
#' @noRd
.mergeMatrices <- function(object, matrices) {
    for (i in seq_along(matrices)) {
        matrices[[i]]$replicate <- object@replicates[[i]]
        matrices[[i]]$condition <- object@conditions[[i]]
    }
    object@interactions <-
        dplyr::bind_rows(matrices) %>%
        dplyr::mutate(chromosome = factor(chromosome)) %>%
        dplyr::mutate(replicate = factor(replicate)) %>%
        dplyr::mutate(condition = factor(condition))
    return(object)
}


## - .parseCool -----------------------------------------------#
## ----------------------------------------------------------------------------#
#' Read interaction matrices in .cool format.
#'
#' This function parses the .cool files provided in
#' \code{object@input}.
#'
#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return object An \code{HiCDOCDataSet} object.
#' @keywords internal
#' @noRd
.parseCool <- function(object) {
    matrices <- pbapply::pblapply(object@input, .parseOneCool)
    return(.mergeMatrices(object, matrices))
}

## - .parseMCool ----------------------------------------------#
## ----------------------------------------------------------------------------#
#' Read interaction matrices in .mcool format.
#'
#' This function parses the .mcool files provided in
#' \code{object@input} at a given resolution.

#' @param object     An \code{HiCDOCDataSet} object.
#' @param resolution The chosen resolution, in base pairs.
#'
#' @return object    An \code{HiCDOCDataSet} object.
#' @keywords internal
#' @noRd
.parseMCool <- function(object, resolution) {
    matrices <-
        pbapply::pblapply(
            object@input,
            .parseOneMCool,
            resolution = resolution
        )
    return(.mergeMatrices(object, matrices))
}

## - .parseOneHiC -----------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Parse a single interaction matrix in .hic format.
#' @param fileName   The file name of the .hic file.
#' @param resolution The chosen resolution (should be present in the
#'                     .hic file).
#'
#' @return object An \code{tibble} storing the (sparse) interaction matrix.
#' @keywords internal
#' @noRd
.parseOneHiC <- function(fileName, resolution = resolution) {
    message(paste0("Parsing .hic file '", fileName, "'."))
    return(parseHiCFile(fileName, resolution))
}

## - .parseHiC ------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Read interaction matrices in .hic format.
#'
#' This function parses the .hic files provided in
#' \code{object@input}.

#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return object An \code{HiCDOCDataSet} object.
#' @keywords internal
#' @noRd
.parseHiC <- function(object) {
    matrices <-
        pbapply::pblapply(
            object@input,
            .parseOneHiC,
            resolution = object@binSize
        )
    matrices <-
        purrr::map2(
            matrices,
            object@replicates,
            ~ dplyr::mutate(.x, replicate = .y)
        )
    matrices <-
        purrr::map2(
            matrices,
            object@conditions,
            ~ dplyr::mutate(.x, condition = .y)
        )
    object@interactions <- dplyr::bind_rows(matrices) %>% dplyr::as_tibble()
    object@interactions %<>%
        dplyr::mutate(
            chromosome = factor(
                chromosome,
                levels = gtools::mixedsort(
                    unique(object@interactions$chromosome)
                )
            )
        ) %>%
        dplyr::mutate(replicate = factor(replicate)) %>%
        dplyr::mutate(condition = factor(condition))
    return(object)
}


## - .parseOneHiCPro --------------------------------------------------------------#
## ----------------------------------------------------------------------------#
#' Parse a single interaction matrix in HiC-Pro format
#'
#' @param vectFiles 2 length vector, of the links to .matrix and .bed files.
#' @return object An \code{tibble} storing the (sparse) interaction matrix.
#' @keywords internal
#' @noRd
.parseOneHiCPro <- function(vectFiles) {
    if (length(vectFiles) != 2) {
        stop("vectFiles must be of length 2")
    }
    matrixFile <- vectFiles[1]
    bedFile <- vectFiles[2]
    dfMatrix <-
        utils::read.table(
            matrixFile,
            header = FALSE,
            stringsAsFactors = FALSE,
            col.names = c("startIndex", "stopIndex", "interaction")
        )
    dfBed <-
        utils::read.table(
            bedFile,
            header = FALSE,
            stringsAsFactors = FALSE,
            col.names = c("chromosome", "start", "end", "index")
        )

    dfMatrix <- dplyr::as_tibble(dfMatrix)
    dfBed <- dplyr::as_tibble(dfBed)

    differences <- sort(table(abs(dfBed$end - dfBed$start)), decreasing = TRUE)
    resolution <- as.numeric(names(differences[1]))

    positions <-
        dfBed %>%
        dplyr::arrange(chromosome, index) %>%
        dplyr::group_by(chromosome) %>%
        dplyr::mutate(bin = index - min(index) + 1) %>%
        dplyr::mutate(end = ifelse(end == dplyr::lead(start), end - 1, end)) %>%
        dplyr::ungroup() %>%
        dplyr::select(chromosome, start, end, bin)

    # Merge data
    dfMatrix <-
        dplyr::left_join(
            dfMatrix,
            dfBed %>% dplyr::select(
                chromosome.1 = chromosome,
                startIndex = index,
                position.1 = start
            ),
            by = "startIndex"
        )
    dfMatrix %<>% dplyr::select(-startIndex)
    dfMatrix <-
        dplyr::left_join(
            dfMatrix,
            dfBed %>% dplyr::select(
                chromosome.2 = chromosome,
                stopIndex = index,
                position.2 = start
            ),
            by = "stopIndex"
        )
    message("Removing inter-chromosome interactions.")
    dfMatrix %<>%
        dplyr::filter(chromosome.1 == chromosome.2) %>%
        dplyr::select(
            chromosome = chromosome.1,
            position.1,
            position.2,
            interaction
        )

    return(list(
        "matrix" = dfMatrix,
        "resolution" = resolution,
        "positions" = positions
    ))
}


## - .parseHiCPro ---------------------------------------------#
## ----------------------------------------------------------------------------#
#' Read interaction matrices in Hic-Pro format.
#'
#' This function parses the HiC-Pro files provided in
#' \code{object@input}.
#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return object An \code{HiCDOCDataSet} object.
#' @keywords internal
#' @noRd
.parseHiCPro <- function(object) {
    matrices <-
        pbapply::pblapply(
            object@input,
            function(x) .parseOneHiCPro(x)
        )
    matrices <-
        purrr::map2(
            matrices,
            object@replicates,
            ~ dplyr::mutate(.x, replicate = .y, .after = position.2)
        )
    matrices <-
        purrr::map2(
            matrices,
            object@conditions,
            ~ dplyr::mutate(.x, condition = .y, .after = position.2)
        )
    object@interactions <- dplyr::bind_rows(matrices) %>% dplyr::as_tibble()
    object@interactions %<>%
        dplyr::mutate(
            chromosome = factor(
                chromosome,
                levels = gtools::mixedsort(
                    unique(object@interactions$chromosome)
                )
            )
        ) %>%
        dplyr::mutate(replicate = factor(replicate)) %>%
        dplyr::mutate(condition = factor(condition))
    return(object)
}
