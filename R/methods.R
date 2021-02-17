#### - chromosomes ---------------------------------------------------------####
## ----------------------------------------------------------------------------#
#' Accessors for the 'chromosomes' slot of an HiCDOCDataSet object
#'
#' The \code{chromosomes} slot contains the names of the chromosomes.
#'
#' @docType methods
#' @name chromosomes
#' @rdname chromosomes
#' @aliases chromosomes chromosomes,HiCDOCDataSet-method
#' @param object An \code{HiCDOCDataSet} object.
#' @return A character vector
#' @examples
#' exp <- HiCDOCExample()
#' chromosomes(exp)
#' @export
setMethod("chromosomes", "HiCDOCDataSet", function(object) object@chromosomes)

#### - interactions --------------------------------------------------------####
## ----------------------------------------------------------------------------#
#' Accessors for the 'interactions' slot of an HiCDOCDataSet object
#'
#' The \code{interactions} slot contains the (transformed) interaction
#' profiles.
#'
#' @docType methods
#' @name interactions
#' @rdname interactions
#' @aliases interactions interactions,HiCDOCDataSet-method
#' @param object An \code{HiCDOCDataSet} object.
#' @return A tibble
#' @examples
#' exp <- HiCDOCExample()
#' interactions(exp)
#' @export
setMethod(
    f = "interactions",
    signature = "HiCDOCDataSet",
    function(object) {
        if (is.null(object@interactions)) {
            return(NULL)
        }
        interactions <- object@interactions %>%
            dplyr::left_join(object@positions %>%
                dplyr::select(chromosome,
                    bin.1 = bin,
                    start.1 = start
                ),
            by = c("chromosome", "bin.1")
            ) %>%
            dplyr::left_join(object@positions %>%
                dplyr::select(chromosome,
                    bin.2 = bin,
                    start.2 = start
                ),
            by = c("chromosome", "bin.2")
            ) %>%
            dplyr::select(
                chromosome,
                bin.1, bin.2, start.1, start.2,
                condition,
                replicate,
                value
            )
        return(interactions)
    }
)

#### - positions ---------------------------------------------------------####
## ----------------------------------------------------------------------------#
#' Accessors for the 'positions' slot of an HiCDOCDataSet object
#'
#' The \code{positions} slot contains the positions of the bins, by chromosome
#'
#' @docType methods
#' @name positions
#' @rdname positions
#' @aliases positions,HiCDOCDataSet-method
#' @param object An \code{HiCDOCDataSet} object.
#' @return A tibble
#' @examples
#' exp <- HiCDOCExample()
#' positions(exp)
#' @export
setMethod("positions", "HiCDOCDataSet", function(object) object@positions)


#### - differences ---------------------------------------------------------####
## ----------------------------------------------------------------------------#
#' Extracts differentially interacting regions of an HiCDOCDataSet object
#'
#' This function extracts the differentially interacting regions from
#' \code{\link{HiCDOCDataSet}}.
#'
#' @docType methods
#' @name differences
#' @rdname differences
#'
#' @aliases differences differences,HiCDOCDataSet-method
#'
#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return A \code{GenomicRanges} object of the selected differentially
#' interacting regions.
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp)
#' differences(exp)
#' @export
setMethod(
    f = "differences",
    signature = "HiCDOCDataSet",
    function(object) {
        if (is.null(object@differences)) {
            return(NULL)
        }
        if (nrow(object@differences) == 0) {
            message("No 'differences' found.")
            return(GenomicRanges::GRanges())
        } else {
            gr <- object@differences %>%
                dplyr::left_join(object@positions,
                    by = c("chromosome", "bin")
                ) %>%
                dplyr::select(
                    chromosome,
                    start,
                    end,
                    condition.1,
                    condition.2,
                    pvalue, padj, direction
                )
            return(GenomicRanges::makeGRangesFromDataFrame(
                gr,
                keep.extra.columns = TRUE,
                ignore.strand = TRUE
            ))
        }
    }
)


#### - concordances --------------------------------------------------------####
## ----------------------------------------------------------------------------#
#' Extracts concordances
#'
#' This function extracts the concordances from \code{\link{HiCDOCDataSet}}.
#'
#' @docType methods
#' @name concordances
#' @rdname concordances
#'
#' @aliases concordances concordances,HiCDOCDataSet-method
#'
#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return A \code{tibble} object of the concordance
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp)
#' concordances(exp)
#' @export
setMethod(
    f = "concordances",
    signature = "HiCDOCDataSet",
    function(object) {
        if (is.null(object@concordances)) {
            return(NULL)
        }
        concordances <- object@concordances %>%
            dplyr::left_join(object@positions, by = c("chromosome", "bin")) %>%
            dplyr::select(
                chromosome,
                start,
                end,
                condition,
                replicate,
                compartment,
                concordance
            )
        return(concordances)
    }
)


#### - compartments --------------------------------------------------------####
## ----------------------------------------------------------------------------#
#' Extracts compartments
#'
#' This function extracts the compartments from \code{\link{HiCDOCDataSet}}.
#'
#' @docType methods
#' @name compartments
#' @rdname compartments
#'
#' @aliases compartments compartments,HiCDOCDataSet-method
#'
#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return A \code{tibble} object of the compartments
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp)
#' compartments(exp)
#' @export
setMethod(
    f = "compartments",
    signature = "HiCDOCDataSet",
    function(object) {
        if (is.null(object@compartments)) {
            return(NULL)
        }
        grl <- object@compartments %>%
            dplyr::left_join(object@positions,
                by = c("chromosome", "bin")
            ) %>%
            dplyr::select(
                chromosome,
                bin,
                start,
                end,
                condition,
                compartment
            )
        return(grl)
        # grl <- object@compartments %>%
        #     mutate(start = position + 1) %>%
        #     mutate(end = start + object@binSize - 1) %>%
        #     select(-position) %>%
        #     mutate(condition = factor(condition)) %>%
        #     rename(compartment = value) %>%
        #     GenomicRanges::makeGRangesListFromDataFrame(
        #         keep.extra.columns = TRUE,
        #         ignore.strand = TRUE,
        #         split.field = "condition"
        #     )
        # grl2 <- lapply(grl, function(x) { split(x, ~ compartment) })
        # grl3 <- as(lapply(grl2, function(x) {
        #     unlist(as(lapply(names(x), function(y) {
        #         z <- GenomicRanges::reduce(x[[y]])
        #         z$compartment <- factor(y)
        #         return(z)
        #     }), "GRangesList"))
        # }), "GRangesList")
        # return(grl3)
    }
)

#### - centroids ----------------------------------------------------------#####
## ----------------------------------------------------------------------------#
#' Extracts centroids
#'
#' This function extracts the centroids from \code{\link{HiCDOCDataSet}}.
#'
#' @docType methods
#' @name centroids
#' @rdname centroids
#'
#' @aliases centroids centroids,HiCDOCDataSet-method
#'
#' @param object An \code{HiCDOCDataSet} object.
#'
#' @return A \code{tibble} object of the centroids
#'
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp)
#' centroids(exp)
#' @export
setMethod("centroids", "HiCDOCDataSet", function(object) object@centroids)


#### - show ----------------------------------------------------------------####
## ----------------------------------------------------------------------------#
#' Show the components of a HiCDOC objects
#'
#' This function informatively display object contents.
#'
#' @docType methods
#' @name show
#' @rdname show
#' @aliases show show,HiCDOCDataSet-method
#' @param object An \code{HiCDOCDataSet} object.
#' @return The \code{show} method informatively display object contents.
#' @export
setMethod(
    f = "show",
    signature = "HiCDOCDataSet",
    function(object) {
        nbCond <- length(unique(object@conditions))
        nbRep <- length(unique(object@replicates))
        cat("Object of class HiCDOCDataSet.\n", "HiCDOC Experiment with:\n")
        cat(
            length(object@chromosomes),
            "chromosomes:",
            object@chromosomes, "\n"
        )
        cat(
            object@totalReplicates, "replications in",
            length(unique(object@conditions)), "conditions\n"
        )
    }
)


#### - parameters ---------------------------------------------------------####
## ----------------------------------------------------------------------------#
#' Accessors for the 'parameters' slot of an HiCDOCDataSet object
#'
#' The \code{parameters} slot holds the parameter values
#' used in an experiment as a named \code{list}. Default values
#' exist for parameters, but these can also be supplied as input
#' values in the \code{useParameters} argument of the \code{\link{HiCDOC}}
#' function or using the assignment function \code{\link{parameters<-}}.
#'
#' Parameters in a HiCDOC experiment.
#'
#' \subsection{Global parameters}{
#'    \describe{
#'        \item{\code{minLengthChr}}{The minimum chromosome
#'                size (in number of bins), to be kept by the function
#'                \code{filterSmallChromosomes()}. Default to 100.}
#'        \item{\code{weakPosThreshold}}{To be kept by the function
#'                \code{filterWeakPositions()}, the bins (positions)
#'                of a chromosome must have a mean value greater than
#'                \code{weakPosThreshold}, on all replicates and all
#'                conditions. The mean is computed on the row of the
#'                reconstructed full interaction matrix for 1 chromosome,
#'                1 condition and 1 replicate. Default to 0.}
#'        \item{\code{sparseThreshold}}{To be kept by the function
#'                \code{filterSparseChromosomes()}, the sparsity (percentage
#'                of empty cells) of the interactions matrix must be lower
#'                than the threshold, on all replicates and all conditions.
#'                Default to 0.95.}
#'       \item{\code{sampleSize}}{The number of bins used when sampling
#'                all the bins in \code{normalizeDistanceEffect()}.
#'                Default to 20000.}
#'       \item{\code{kMeansIterations}}{The maximum number of 2-means
#'                iterations.}
#'       \item{\code{kMeansDelta}}{The stop criterion of convergence of
#'                the 2-means method.}
#'       \item{\code{kMeansRestarts}}{The maximum number of restarts
#'                for the 2-means.}
#'    }
#' }
#' @docType methods
#' @name parameters
#' @rdname parameters
#' @aliases parameters parameters,HiCDOCDataSet-method
#' parameters<- parameters<-,HiCDOCDataSet-method
#' @param object An \code{HiCDOCDataSet} object.
#' @param value A named \code{list} containing valid parameters. See details.
#' @return The named list of the parameters used in the analysis.
#' @seealso
#' \code{useParameters} argument in \code{\link{HiCDOC}} function.
#' @examples
#' exp <- HiCDOCExample()
#' parameters(exp)
#' parameters(exp) <- list("weakPosThreshold" = 100)
#' @export
setMethod(
    f = "parameters",
    signature = "HiCDOCDataSet",
    definition = function(object) {
        testSlotsHiCDOC(object, "parameters")

        parameters <- object@parameters
        kmparam <- which(grepl("kMeans", names(parameters)))

        cat("------------ Global parameters ------------\n")
        matparam <- as.matrix(parameters[-kmparam])
        colnames(matparam)[1] <- "Value"
        print(matparam)

        cat("\n----- Constrained K-means parameters ----- \n")
        matparam <- as.matrix(parameters[kmparam])
        cat("\n")
        colnames(matparam)[1] <- "Value"
        print(matparam)
        return(invisible(parameters))
    }
)


#' @name parameters
#' @rdname parameters
#' @exportMethod "parameters<-"
setReplaceMethod(
    "parameters",
    signature(object = "HiCDOCDataSet", value = "ANY"),
    function(object, value) {
        ## - checking input value ---------------------------------#
        ## --------------------------------------------------------#
        defaultParNames <- names(HiCDOCDefaultParameters)
        currentPar <- object@parameters
        if (!is(value, "list")) {
            stop(
                "'value' must be a named list. See ",
                "help(parameters) for details.",
                call. = FALSE
            )
        }

        valueNames <- names(value)

        if (any(duplicated(valueNames))) {
            stop(
                "duplicate name parameters in 'value'. See ",
                "help(parameters) for details.",
                call. = FALSE
            )
        }

        if (!all(valueNames %in% defaultParNames)) {
            stop(
                "'value' must be a named list of valid ",
                "parameters. See help(parameters) for details.",
                call. = FALSE
            )
        }

        ## - replacing individual parameters
        currentPar[valueNames] <- value

        ## - end check -------------------------------------------#

        object@parameters <- currentPar
        return(object)
    }
)
