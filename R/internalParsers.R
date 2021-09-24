#' @description
#' Parses interactions in tabular format and fills the conditions, replicates,
#' and interactions slots of the provided \code{\link{HiCDOCDataSet}}.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param sep
#' The separator of the tabular file.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.parseTabular <- function(object, sep = "\t") {
    
    message("Parsing '", object@input, "'.")
    
    interactions <-
        data.table::fread(
            file = object@input,
            sep = sep,
            header = TRUE,
            # comment.char = "#",
            check.names = FALSE,
            data.table = FALSE
        ) %>%
        dplyr::as_tibble()
    
    if (colnames(interactions)[[1]] != "chromosome") {
        stop(
            "First column of the input file must be named 'chromosome'.",
            call. = FALSE
        )
    }
    if (colnames(interactions)[[2]] != "position 1") {
        stop(
            "Second column of the input file must be named 'position 1'.",
            call. = FALSE
        )
    }
    if (colnames(interactions)[[3]] != "position 2") {
        stop(
            "Third column of the input file must be named 'position 2'.",
            call. = FALSE
        )
    }
    
    columns <- colnames(interactions)[4:ncol(interactions)]
    
    if (!all(grepl("^.+?\\..+$", columns))) {
        stop(
            "Fourth to last column of the input file must be named 'C.R', ",
            "with C the condition number/name and R the replicate number/name.",
            call. = FALSE
        )
    }
    
    object@conditions <- gsub("^(.+?)\\..+$", "\\1", columns)
    object@replicates <- gsub("^.+?\\.(.+)$", "\\1", columns)
    
    object@interactions <-
        tidyr::pivot_longer(
            interactions,
            dplyr::all_of(columns),
            names_to = "condition.replicate",
            values_to = "interaction"
        ) %>%
        dplyr::mutate(CR = strsplit(condition.replicate, ".", fixed = TRUE)) %>%
        dplyr::mutate(condition = purrr::map_chr(CR, 1),
                      replicate = purrr::map_chr(CR, 2)) %>%
        dplyr::select(-CR) %>%
        dplyr::rename(position.1 = `position 1`) %>%
        dplyr::rename(position.2 = `position 2`)
    
    return(object)
}

#' @description
#' Parses a single interactions file in \code{.cool} or \code{.mcool} format.
#'
#' @param path
#' The path to the interactions file.
#' @param binSize
#' The resolution (span of each position in number of bases). Optionally
#' provided to select the appropriate resolution in \code{.mcool} files.
#' Defaults to NULL.
#'
#' @return
#' A tibble of interactions.
#'
#' @keywords internal
#' @noRd
.parseOneCool <- function(path, binSize = NULL) {
    
    message("\nParsing '", path, "'.")
    
    uri <- function(path) {
        if (!is.numeric(binSize)) return(path)
        return(
            paste(
                "resolutions",
                format(binSize, scientific = FALSE),
                path,
                sep = "/"
            )
        )
    }
    
    bins <-
        dplyr::tibble(
            chromosome = factor(
                rhdf5::h5read(file = path, name = uri("bins/chrom"))
            ),
            start = rhdf5::h5read(file = path, name = uri("bins/start")),
            end = rhdf5::h5read(file = path, name = uri("bins/end"))
        )
    
    step <- bins$end - bins$start
    if (length(step) < 0.9) {
        stop("Cannot parse '", path, "': fixed width only.", call. = FALSE)
    }
    step <- max(step)
    
    bins %<>%
        dplyr::select(-end) %>%
        dplyr::rename(position = start) %>%
        dplyr::mutate(ide = seq_len(nrow(bins)) - 1) %>%
        dplyr::select(ide, chromosome, position)
    
    rownames(bins) <- NULL
    
    interactions <-
        dplyr::tibble(
            id1 = rhdf5::h5read(file = path, name = uri("pixels/bin1_id")),
            id2 = rhdf5::h5read(file = path, name = uri("pixels/bin2_id")),
            interaction = rhdf5::h5read(file = path, name = uri("pixels/count"))
        ) %>%
        dplyr::left_join(bins, by = c("id1" = "ide")) %>%
        dplyr::rename(chromosome.1 = chromosome, position.1 = position) %>%
        dplyr::select(-id1) %>%
        dplyr::left_join(bins, by = c("id2" = "ide")) %>%
        dplyr::rename(chromosome.2 = chromosome, position.2 = position) %>%
        dplyr::select(-id2) %>%
        dplyr::filter(chromosome.1 == chromosome.2) %>%
        dplyr::select(-chromosome.2) %>%
        dplyr::rename(chromosome = chromosome.1) %>%
        dplyr::select(chromosome, position.1, position.2, interaction)
    
    return(interactions)
}

#' @description
#' Parses interactions in \code{.cool} or \code{.mcool} format and fills the
#' interactions slot of the provided \code{\link{HiCDOCDataSet}}.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param binSize
#' The resolution (span of each position in number of bases). Optionally
#' provided to select the appropriate resolution in \code{.mcool} files.
#' Defaults to NULL.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.parseCool <- function(object, binSize = NULL) {
    
    interactions <-
        pbapply::pblapply(
            object@input,
            .parseOneCool,
            binSize = binSize
        )
    
    for (i in seq_along(interactions)) {
        interactions[[i]]$replicate <- object@replicates[[i]]
        interactions[[i]]$condition <- object@conditions[[i]]
    }
    
    object@interactions <-
        dplyr::bind_rows(interactions) %>%
        dplyr::as_tibble()
    
    return(object)
}

#' @description
#' Parses a single interactions file in \code{.hic} format. Calls the C++
#' \code{parseHiCFile} parser.
#'
#' @param path
#' The path to the interactions file.
#' @param binSize
#' The resolution (span of each position in number of bases) to select within
#' the \code{.hic} file.
#'
#' @return
#' A dataframe of interactions.
#'
#' @keywords internal
#' @noRd
.parseOneHiC <- function(path, binSize) {
    message("\nParsing '", path, "'.")
    return(parseHiCFile(path, binSize))
}

#' @description
#' Parses interactions in \code{.hic} format and fills the interactions slots of
#' the provided \code{\link{HiCDOCDataSet}}.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param binSize
#' The resolution (span of each position in number of bases) to select within
#' the \code{.hic} files.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.parseHiC <- function(object, binSize) {
    
    interactions <-
        pbapply::pblapply(
            object@input,
            .parseOneHiC,
            binSize = binSize
        ) %>%
        purrr::map2(
            object@replicates,
            ~ dplyr::mutate(.x, replicate = .y)
        ) %>%
        purrr::map2(
            object@conditions,
            ~ dplyr::mutate(.x, condition = .y)
        )
    
    object@interactions <-
        dplyr::bind_rows(interactions) %>%
        dplyr::as_tibble()
    
    return(object)
}

#' @description
#' Parses a single pair of \code{.matrix} and \code{.bed} files.
#'
#' @param matrixPath
#' The path to the interactions matrix file.
#' @param bedPath
#' The path to the bed file.
#'
#' @return
#' A tibble of interactions.
#'
#' @keywords internal
#' @noRd
.parseOneHiCPro <- function(matrixPath, bedPath) {
    
    message("\nParsing '", matrixPath, "' and '", bedPath, "'.")
    
    interactions <-
        data.table::fread(
            matrixPath,
            header = FALSE,
            stringsAsFactors = FALSE,
            col.names = c("startIndex", "stopIndex", "interaction"),
            data.table = TRUE
        )
    
    bed <-
        data.table::fread(
            bedPath,
            header = FALSE,
            stringsAsFactors = FALSE,
            col.names = c("chromosome", "start", "end", "index"),
            data.table = TRUE
        )
    
    interactions <- merge(interactions,
                          bed[,.(chromosome.1 = chromosome, 
                                 startIndex = index, 
                                 position.1 = start)], 
                          all.x=TRUE, 
                          by = "startIndex")
    interactions <- merge(interactions,
                          bed[,.(chromosome.2 = chromosome, 
                                 stopIndex = index, 
                                 position.2 = start)], 
                          all.x=TRUE, 
                          by = "stopIndex")
    interactions <- interactions[chromosome.1 == chromosome.2,
                                 .(chromosome = chromosome.1,
                                   position.1,
                                   position.2,
                                   interaction)]
    
    return(interactions)
}

#' @description
#' Parses interactions in pairs of \code{.matrix} and \code{.bed} files and
#' fills the interactions slots of the provided \code{\link{HiCDOCDataSet}}.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.parseHiCPro <- function(object) {
    
    interactions <-
        pbapply::pblapply(
            object@input,
            function(paths) .parseOneHiCPro(paths[1], paths[2])
        ) %>%
        purrr::map2(
            object@replicates,
            .f = function(x, y) x[, replicate:=y]
        ) %>%
        purrr::map2(
            object@conditions,
            .f = function(x, y) x[, condition:=y]
        )
    
    object@interactions <-
        data.table::rbindlist(interactions) %>%
        dplyr::as_tibble()
    
    return(object)
}
