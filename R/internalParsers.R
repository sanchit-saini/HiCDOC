#' @description
#' Parses interactions in tabular format and fills the conditions, replicates,
#' and interactions slots of the provided \code{\link{HiCDOCDataSet}}.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.parseTabular <- function(object) {

    message(paste0("Parsing '", object@input, "'."))

    interactions <-
        utils::read.table(
            file = object@input,
            sep = "\t",
            header = TRUE,
            comment.char = "#",
            check.names = FALSE
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

    conditions.replicates <- colnames(interactions)[4:ncol(interactions)]

    if (!all(grepl("^.+?\\..+$", conditions.replicates))) {
        stop(
            "Fourth to last column of the input file must be named 'C.R', ",
            "with C the condition number/name and R the replicate number/name.",
            call. = FALSE
        )
    }

    object@conditions <- gsub("^(.+?)\\..+$", "\\1", conditions.replicates)
    object@replicates <- gsub("^.+?\\.(.+)$", "\\1", conditions.replicates)

    object@interactions <-
        tidyr::gather(
            interactions,
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

#' @description
#' Parses a single interactions file in \code{.cool} or \code{.mcool} format.
#'
#' @param path
#' The path to the interactions file.
#' @param resolution
#' The number of bases per bin. Optionally provided to select the appropriate
#' interactions in a \code{.mcool} file. Defaults to NULL.
#'
#' @return
#' A tibble of interactions.
#'
#' @keywords internal
#' @noRd
.parseOneCool <- function(path, resolution = NULL) {

    message(paste0("\nParsing '", path, "'."))

    uri <- function(path) {
        if (!is.numeric(resolution)) return(path)
        return(
            paste(
                "resolutions",
                format(resolution, scientific = FALSE),
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
        dplyr::select(-(end)) %>%
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
#' @param resolution
#' The number of bases per bin. Optionally provided to select the appropriate
#' interactions in \code{.mcool} files. Defaults to NULL.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.parseCool <- function(object, resolution = NULL) {

    interactions <-
        pbapply::pblapply(
            object@input,
            .parseOneCool,
            resolution = resolution
        )

    for (i in seq_along(interactions)) {
        interactions[[i]]$replicate <- object@replicates[[i]]
        interactions[[i]]$condition <- object@conditions[[i]]
    }

    object@interactions <-
        dplyr::bind_rows(interactions) %>%
        dplyr::mutate(
            chromosome = factor(chromosome),
            condition = factor(condition),
            replicate = factor(replicate),
            interaction = as.numeric(interaction)
        )

    return(object)
}

#' @description
#' Parses a single interactions file in \code{.hic} format. Calls the C++
#' \code{parseHiCFile} parser.
#'
#' @param path
#' The path to the interactions file.
#' @param resolution
#' The number of bases per bin. Selects the appropriate interactions in the
#' \code{.hic} file.
#'
#' @return
#' A dataframe of interactions.
#'
#' @keywords internal
#' @noRd
.parseOneHiC <- function(path, resolution) {
    message(paste0("\nParsing '", path, "'."))
    return(parseHiCFile(path, resolution))
}

#' @description
#' Parses interactions in \code{.hic} format and fills the interactions slots of
#' the provided \code{\link{HiCDOCDataSet}}.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param resolution
#' The number of bases per bin. Selects the appropriate interactions in the
#' \code{.hic} files.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.parseHiC <- function(object, resolution) {

    interactions <-
        pbapply::pblapply(
            object@input,
            .parseOneHiC,
            resolution = resolution
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
        dplyr::as_tibble() %>%
        dplyr::mutate(
            chromosome = factor(chromosome),
            condition = factor(condition),
            replicate = factor(replicate),
            interaction = as.numeric(interaction)
        )

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

    message(paste0("\nParsing '", matrixPath, "' and '", bedPath, "'."))

    interactions <-
        utils::read.table(
            matrixPath,
            header = FALSE,
            stringsAsFactors = FALSE,
            col.names = c("startIndex", "stopIndex", "interaction")
        ) %>%
        dplyr::as_tibble()

    bed <-
        utils::read.table(
            bedPath,
            header = FALSE,
            stringsAsFactors = FALSE,
            col.names = c("chromosome", "start", "end", "index")
        ) %>%
        dplyr::as_tibble()

    interactions %<>%
        dplyr::left_join(
            bed %>% dplyr::select(
                chromosome.1 = chromosome,
                startIndex = index,
                position.1 = start
            ),
            by = "startIndex"
        ) %>%
        dplyr::select(-startIndex) %>%
        dplyr::left_join(
            bed %>% dplyr::select(
                chromosome.2 = chromosome,
                stopIndex = index,
                position.2 = start
            ),
            by = "stopIndex"
        ) %>%
        dplyr::filter(chromosome.1 == chromosome.2) %>%
        dplyr::select(
            chromosome = chromosome.1,
            position.1,
            position.2,
            interaction
        )

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
            ~ dplyr::mutate(.x, replicate = .y)
        ) %>%
        purrr::map2(
            object@conditions,
            ~ dplyr::mutate(.x, condition = .y)
        )

    object@interactions <-
        dplyr::bind_rows(interactions) %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(
            chromosome = factor(chromosome),
            condition = factor(condition),
            replicate = factor(replicate),
            interaction = as.numeric(interaction)
        )

    return(object)
}
