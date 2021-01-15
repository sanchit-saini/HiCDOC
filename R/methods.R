####- chromosomes ---------------------------------------------------------####
##----------------------------------------------------------------------------#
#' Accessors for the 'chromosomes' slot of an HiCDOCDataSet object
#'
#' The \code{chromosomes} slot contains the names of the chromosomes, eventually
#' filtred after \code{filterSmallChromosomes()}
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
#'
#' @export
setMethod(
    f = "chromosomes",
    signature = "HiCDOCDataSet",
    definition = function(object) {
        testSlotsHiCDOC(object, "chromosomes")
        return(object@chromosomes)
    }
)

####- interactions --------------------------------------------------------####
##----------------------------------------------------------------------------#
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
#'
#' @export
setMethod(
    f = "interactions",
    signature = "HiCDOCDataSet",
    definition = function(object) {
        testSlotsHiCDOC(object, "interactions")
        interactions <- object@interactions %>%
                dplyr::left_join(object@positions %>%
                                     select(chromosome, 
                                            bin.1 = bin, 
                                            start.1 = start), 
                                 by=c("chromosome", "bin.1")) %>%
                dplyr::left_join(object@positions %>%
                                     select(chromosome, 
                                            bin.2 = bin, 
                                            start.2 = start), 
                                 by=c("chromosome", "bin.2")) %>%
                dplyr::select(chromosome, 
                              start.1, start.2,
                              condition, 
                              replicate,
                              value)
        return(interactions)
    }
)


####- differences ---------------------------------------------------------####
##----------------------------------------------------------------------------#
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
#'
#' @export
setMethod(
    f = "differences",
    signature = "HiCDOCDataSet",
    definition = function(object) {
        testSlotsHiCDOC(object, "differences")
        if (nrow(object@differences) == 0) {
            message("No 'differences' found.")
            return(GenomicRanges::GRanges())
        } else {
            gr <- object@differences %>%
                dplyr::left_join(object@positions, by=c("chromosome", "bin")) %>%
                dplyr::select(chromosome, 
                              start, 
                              end, 
                              condition.1,
                              condition.2,
                              pvalue, padj, direction) 
            return (GenomicRanges::makeGRangesFromDataFrame(
                gr,
                keep.extra.columns = TRUE,
                ignore.strand = TRUE
            ))
        }
    }
)


####- concordances --------------------------------------------------------####
##----------------------------------------------------------------------------#
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
#'
#' @export
setMethod(
    f = "concordances",
    signature = "HiCDOCDataSet",
    definition = function(object) {
        testSlotsHiCDOC(object, "concordances")
        concordances <- object@concordances %>% 
            dplyr::left_join(object@positions, by=c("chromosome", "bin")) %>%
            dplyr::select(chromosome, 
                   start, 
                   end, 
                   condition, 
                   replicate, 
                   compartment, 
                   concordance)
        return(concordances)
    }
)


####- compartments --------------------------------------------------------####
##----------------------------------------------------------------------------#
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
#'
#' @export
setMethod(
    f = "compartments",
    signature = "HiCDOCDataSet",
    definition = function(object) {
        if (is.null(object@compartments)) {
            message(
                "No 'compartments' slot found in the ",
                "HiCDOCDataSet object. Run HiCDOC first."
            )
        } else {
            grl <- object@compartments %>%
                    dplyr::left_join(object@positions, by=c("chromosome", "bin")) %>%
                    dplyr::select(chromosome, 
                                  bin, 
                                  start, 
                                  end, 
                                  condition, 
                                  compartment)
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
    }
)

####- centroids ----------------------------------------------------------#####
##----------------------------------------------------------------------------#
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
#'
#' @export
setMethod(
    f = "centroids",
    signature = "HiCDOCDataSet",
    definition = function(object) {
        testSlotsHiCDOC(object, "centroids")
        return(object@centroids)
    }
)



####- parameters ---------------------------------------------------------####
##----------------------------------------------------------------------------#
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
#'        \item{\code{minDepth}}{The cutoff to filter the base-level
#'                coverage. Bases where at least one sample has (normalized)
#'                coverage greater than \code{minDepth} be been retained.
#'                Default to \code{10}.}
#'       \item{\code{sampleSize}}{The number of bins used when sampling
#'                all the bins.}
#'       \item{\code{loessSpan}}{The optimal span value used for the
#'                diagonal normalization.}
#'       \item{\code{kMeansIterations}}{The maximum number of 2-means
#'                iterations.}
#'       \item{\code{kMeansDelta}}{The stop criterion of convergence of
#'                the 2-means method.}
#'       \item{\code{kMeansRestarts}}{The maximum number of restarts
#'                for the 2-means.}
#'       \item{\code{minLengthChr}}{The minimum chromosome
#'                size (in number of bins), to be kept.}
#'       \item{\code{weakPosThreshold}}{To be kept, the bins (positions)
#'                of a chromosome must have a mean value greater than
#'                \code{weakPosThreshold}, on all replicates and all
#'                conditions. The mean is computed on the row of the
#'                reconstructed full interaction matrix for 1 chromosome,
#'                1 condition and 1 replicate.}
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
#' exp <- HiCDOC(exp)
#' print(parameters(exp))
#'
#' @export
setMethod(
    f = "parameters",
    signature = "HiCDOCDataSet",
    definition = function(object) {
        testSlotsHiCDOC(object, "parameters")
        # object@parameters
        # class(object@parameters) <- "HiCDOCParameters"
        return(object@parameters)
    }
)


#' @name parameters
#' @rdname parameters
#' @exportMethod "parameters<-"
setReplaceMethod(
    "parameters",
    signature(object = "HiCDOCDataSet", value = "ANY"),
    function(object, value) {

        ##- checking input value ---------------------------------#
        ##--------------------------------------------------------#
        defaultParNames <- names(HiCDOCDefaultParameters)

        if (!is.null(object@parameters)) {
            HiCDOCDefaultParameters <- object@parameters
        }

        if (!is(value, "list")) {
            print(value)
            print(typeof(value))
            print(class(value))
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

        ##- individual parameters
        HiCDOCDefaultParameters[valueNames] <- value
        checkParameters(HiCDOCDefaultParameters)

        ##- end check -------------------------------------------#

        object@parameters <- HiCDOCDefaultParameters
        object
    }
)

##- show ---------------------------------------------------------------------#
##----------------------------------------------------------------------------#
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
    definition = function(object) {
        nbCond <- length(unique(object@conditions))
        nbRep <- length(unique(object@replicates))
        cat("Object of class HiCDOCDataSet.\n", "HiCDOC Experiment with:\n")
        cat(length(object@chromosomes), "chromosomes:", object@chromosomes, "\n")
        cat(object@totalReplicates, "replications in",
            length(unique(object@conditions)), "conditions\n")
    }
)


##- print method for parameters ----------------------------------------------#
##----------------------------------------------------------------------------#
#' Dispatch print method for the parameters used by an \code{HiCDOCDataSet}
#' object.
#'
#' @docType methods
#' @name parameters
#' @rdname parameters
#' @aliases parameters parameters,HiCDOCDataSet-method
#' @param x The first element of the parameters used by an \code{HiCDOCDataSet}
#' object
#' @param ... The other elements of the parameters
#' @examples
#' exp <- HiCDOCExample()
#' exp <- HiCDOC(exp)
#' print(parameters(exp))
#'
#' @export
printHiCDOCParameters <- function(object) {
    parameters  <- unlist(object@parameters)
    kmparam <- which(grepl("kMeans", names(parameters)))
    cat("\n Global parameters: \n", "------------------ \n")
    df <- data.frame(value = parameters[-kmparam])
    print(df)

    cat("\n Constrained K-means parameters: \n", "---------------------- \n")
    df <- data.frame(value = parameters[kmparam])
    print(df)
}
