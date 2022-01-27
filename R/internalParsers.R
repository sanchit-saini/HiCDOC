#' @description
#' Parses interactions in tabular format and fills the conditions, replicates,
#' and interactions slots of the provided \code{\link{InteractionSet}}.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param input
#' The name of the input file
#' @param conditions
#' The names of the conditions
#' (infered from the table if not specified).
#' @param replicates
#' The names of the replicates
#' (infered from the table if not specified).
#'
#' @return
#' An \code{\link{InteractionSet}}.
#'
#' @keywords internal
#' @noRd
.setFromTabular <- function(tabular, conditions = NULL, replicates = NULL) {

    if (colnames(tabular)[1] != "chromosome") {
        stop(
            "First column of the input file must be named 'chromosome'.",
            call. = FALSE
        )
    }
    if (colnames(tabular)[2] == "position.1") {
        data.table::setnames(tabular, "position.1", "position 1")
    }
    if (colnames(tabular)[2] != "position 1") {
        stop(
            "Second column of the input file must be named 'position 1'.",
            call. = FALSE
        )
    }
    if (colnames(tabular)[3] == "position.2") {
        data.table::setnames(tabular, "position.2", "position 2")
    }
    if (colnames(tabular)[3] != "position 2") {
        stop(
            "Third column of the input file must be named 'position 2'.",
            call. = FALSE
        )
    }
    if (is.null(conditions) != is.null(replicates)) {
        stop(
            "Conditions and replicates should be both NULL, or none.",
            call. = FALSE
        )
    }
    tabular[,chromosome := as.character(chromosome)]
    tabular[,chromosome := factor(chromosome, 
                                       levels=gtools::mixedsort(unique(chromosome)))]
    setorder(tabular, chromosome, `position 1`, `position 2`)
    # Assays part, fill with NA
    assays <- as.matrix(tabular[,4:ncol(tabular), drop=FALSE])

    if (! is.null(conditions)) {
        if ((length(conditions) != ncol(assays)) | (length(replicates) != ncol(assays))) {
            stop(
                "Number of conditions and replicates should match the number of counts in the matrix.",
                call. = FALSE
            )
        }
    } else {
        if (!all(grepl("^.+?\\..+$", colnames(assays)))) {
            stop(
                "Fourth to last column of the input file must be named 'C.R', ",
                "with C the condition number/name and R the replicate number/name.",
                call. = FALSE
            )
        }
    }
    assays[assays == 0] <- NA
    
    # GInteraction part
    tabular <- tabular[,1:3]
    data.table::setnames(tabular, "position 1", "bin.1")
    data.table::setnames(tabular, "position 2", "bin.2")
    
    diag <- (tabular$bin.1 == tabular$bin.2)
    binSize <- modeVector(abs(tabular[!diag,]$bin.1 - 
                                      tabular[!diag,]$bin.2))

    tabular[,bin.1 := bin.1/binSize]
    tabular[,bin.2 := bin.2/binSize]
    
    allRegions <- data.table::melt(tabular[,.(chromosome, bin.1, bin.2)],
                                   id.vars = "chromosome", 
                                   value.name = "indexC")
    allRegions[,variable := NULL]
    allRegions <- unique(allRegions)
    setorder(allRegions, chromosome, indexC)
    
    # Constructing unique index for all chromosomes, taking into account the difference in bins. 
    allRegions[,index := indexC - data.table::shift(indexC, fill = 0), by=.(chromosome)]
    allRegions[, index := cumsum(index)]
    allRegions[, end := (indexC+1) * binSize -1]
    allRegions[, start := (indexC) * binSize]
    data.table::setcolorder(allRegions, 
                            c("chromosome", "start", "end", "index", "indexC"))
    
    tabular <- 
        data.table::merge.data.table(
            tabular, 
            allRegions[,.(chromosome, startIndex = index, bin.1 = indexC)], 
            all.x=TRUE, 
            sort=FALSE,
            by=c("chromosome", "bin.1"))
    tabular <- 
        data.table::merge.data.table(
            tabular, 
            allRegions[,.(chromosome, stopIndex = index, bin.2 = indexC)], 
            all.x=TRUE, 
            sort=FALSE,
            by=c("chromosome", "bin.2"))
    tabular[, bin.1 := NULL]
    tabular[, bin.2 := NULL]
    
    allRegions[,indexC := NULL]
    order1 <- match(tabular$startIndex, allRegions$index)
    order2 <- match(tabular$stopIndex, allRegions$index)
    allRegions <- GenomicRanges::GRanges(allRegions)
    
    gi <- InteractionSet::GInteractions(
        allRegions[order1],
        allRegions[order2],
        regions = allRegions,
        mode="strict")
    
    iset <- InteractionSet::InteractionSet(
        assays = assays,
        interactions = gi)
    if (is.null(conditions)) {
        SummarizedExperiment::colData(iset) <- 
            S4Vectors::DataFrame(
                "condition" = gsub("^(.+?)\\..+$", "\\1", colnames(assays)),
                "replicate" = gsub("^.+?\\.(.+)$", "\\1", colnames(assays))
                )
    } else {
        SummarizedExperiment::colData(iset) <- 
            S4Vectors::DataFrame(
                "condition" = conditions,
                "replicate" = replicates
                )
    }
    
    # Keep only intra-chromosomal interactions
    iset <- iset[InteractionSet::intrachr(iset),]
    
    # Remove zero rows
    zeros <- (rowSums(SummarizedExperiment::assay(iset), na.rm=TRUE) == 0)
    iset <- iset[!zeros,]
    
    return(iset)
}

#' @description
#' Read the file, and fills it using \code{\link{.setFromTabular}}.
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
.parseTabular <- function(input, sep = "\t") {
    
    message("Parsing '", input, "'.")
    
    interactions <-
        data.table::fread(
            file = input,
            sep = sep,
            header = TRUE,
            # comment.char = "#",
            check.names = FALSE,
            data.table = TRUE,
            stringsAsFactors = FALSE
        ) 

    iset   <- .setFromTabular(interactions)
    object <- new("HiCDOCDataSet", iset, input = input)
    
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
.parseOneCool <- function(path, binSize = NA, replicate, condition) {
    
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
        data.table::data.table(
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
    
    bins[,start := as.integer(start)]
    bins[,end := as.integer(end)]
    
    # bins[, end := NULL]bins
    setorder(bins, chromosome, start, end)
    bins[, index := seq_len(nrow(bins))]
    
    interactions <-
        data.table::data.table(
            id1 = rhdf5::h5read(file = path, name = uri("pixels/bin1_id")),
            id2 = rhdf5::h5read(file = path, name = uri("pixels/bin2_id")),
            interaction = rhdf5::h5read(file = path, name = uri("pixels/count"))
        )
    interactions[, id1:=as.integer(id1) + 1]
    interactions[, id2:=as.integer(id2) + 1]
    interactions[, interaction:=as.numeric(interaction)]
    
    order1 <- match(interactions$id1, bins$index)
    order2 <- match(interactions$id2, bins$index)
    allRegions <- GenomicRanges::GRanges(bins)
    
    # GInteractions part
    gi <- InteractionSet::GInteractions(
        allRegions[order1],
        allRegions[order2],
        regions = allRegions,
        mode="strict")
    
    iset <- InteractionSet::InteractionSet(
        assays = as.matrix(interactions$interaction, ncol=1), 
        interactions = gi,  
        colData=S4Vectors::DataFrame("condition" = condition, 
                                     "replicate" = replicate)) 
    
    # Keep only intra-chromosomal interactions
    iset <- iset[InteractionSet::intrachr(iset),]
    
    # Remove zero rows
    zeros <- (rowSums(SummarizedExperiment::assay(iset), na.rm=TRUE) == 0)
    iset <- iset[!zeros,]
    return(iset)
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
.parseCool <- function(object, binSize = NA, replicates, conditions) {
    
    isetCool <-
        pbapply::pbmapply(
            .parseOneCool,
            path      = object@input,
            binSize   = binSize,
            condition = conditions,
	    replicate = replicates
        )
    
    mergedIsetCool <- Reduce(f = .mergeInteractionSet, x = isetCool)
    
    new("HiCDOCDataSet", 
        mergedIsetCool, 
        input = object@input)
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
#' An \code{\link{InteractionSet}}.
#'
#' @keywords internal
#' @noRd
.parseOneHiC <- function(path, binSize, condition, replicate) {
    message("\nParsing '", path, "'.")
    interactions <- parseHiCFile(path, binSize)
    # Automagical stuff to transform Rcpp DataFrame to data.table
    interactions <- data.table::setalloccol(interactions)
    # Create InteractionSet object
    interactions <- .setFromTabular(interactions, condition, replicate)
    return(interactions)
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
.parseHiC <- function(object, binSize, replicates, conditions) {
    
    isetHic <-
        pbapply::pbmapply(
            .parseOneHiC,
            path      = object@input,
            binSize   = binSize,
            condition = conditions,
	        replicate = replicates
        )
    
    mergedIsetHic <- Reduce(f = .mergeInteractionSet, x = isetHic)

    new("HiCDOCDataSet", 
        mergedIsetHic, 
        input = object@input)
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
.parseOneHiCPro <- function(matrixPath, bedPath, replicate, condition) {
    
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
    
    minIndex <- min(bed[,index])
    setorder(bed, chromosome, start, end)
    
    order1 <- match(interactions$startIndex, bed$index)
    order2 <- match(interactions$stopIndex, bed$index)
    allRegions <- GenomicRanges::GRanges(bed)
    
    gi <- InteractionSet::GInteractions(
        allRegions[order1],
        allRegions[order2],
        regions = allRegions,
        mode="strict")
    
    iset <- InteractionSet::InteractionSet(
        assays = as.matrix(interactions$interaction, ncol=1), 
        interactions = gi,  
        colData=S4Vectors::DataFrame("condition" = condition, 
                                     "replicate" = replicate)) 
    
    # Keep only intra-chromosomal interactions
    iset <- iset[InteractionSet::intrachr(iset),]
    
    # Remove zero rows
    zeros <- (rowSums(SummarizedExperiment::assay(iset), na.rm=TRUE) == 0)
    iset <- iset[!zeros,]
    
    return(iset)
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
.parseHiCPro <- function(object, replicates, conditions) {

    matrixPaths <- lapply(object@input, `[[`, 1)
    bedPaths    <- lapply(object@input, `[[`, 2)

    isetHic <-
        # pbapply::pbmapply(
        mapply(
            .parseOneHiCPro,
            matrixPaths,
            bedPaths,
            replicates,
            conditions
        )
    
    mergedIsetHic <- Reduce(f = .mergeInteractionSet, x = isetHic)
    
    object <- new("HiCDOCDataSet", 
                  mergedIsetHic, 
                  input = object@input)
    
    return(object)
}
