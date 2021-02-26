#### - chromosomes ------------------------------------------------------####
## -------------------------------------------------------------------------#
#' @describeIn HiCDOCDataSet access and print the chromosomes
#' @seealso
#' \code{\link[HiCDOC]{parameters}} function to access the parameters of a 
#' HiCDOCDataSet object.
#' @param object a HiCDOCDataSet object.
#' @return For \code{chromosomes}: a character vector, with the name 
#' of chromosomes
#' @examples 
#' object <- HiCDOCExample()
#' ## Print the object
#' object
#' 
#' ## Accessors
#' chromosomes(object)
#' conditions(object)
#' replicates(object)
#' interactions(object)
#' positions(object)
#' parameters(object)
#' 
#' ## Accessors after a run of detectCompartments
#' object <- HiCDOC(object)
#' differences(object)
#' concordances(object)
#' compartments(object)
#' centroids(object)
#' @export
setMethod("chromosomes", "HiCDOCDataSet", function(object) object@chromosomes)

#### - conditions -------------------------------------------------------####
## --------------------------------------------------------------------------#
#' @describeIn HiCDOCDataSet access and print the conditions
#' @export
#' @return For \code{conditions}: a character vector, with the name of 
#' conditions, with repetition over the replicates.
setMethod("conditions", "HiCDOCDataSet", function(object) object@conditions)

#### - replicates -------------------------------------------------------####
## -------------------------------------------------------------------------#
#' @describeIn HiCDOCDataSet access and print the replicates
#' @export
#' @return For \code{replicates}: a character vector, with the name of 
#' conditions, with repetition over the conditions
setMethod("replicates", "HiCDOCDataSet", function(object) object@replicates)

#### - interactions ------------------------------------------------------####
## --------------------------------------------------------------------------#
#' @describeIn HiCDOCDataSet access and print the interactions matrix
#' @export
#' @return For \code{interactions}: a tibble of the interactions on all 
#' the choromosomes, conditions and replicates
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
## --------------------------------------------------------------------------#
#' @describeIn HiCDOCDataSet access and print the positions
#' @export
#' @return For \code{positions}: a tibble of the positions corresponding
#' to the bins, for each chromosome
setMethod("positions", "HiCDOCDataSet", function(object) object@positions)


#### - differences -------------------------------------------------------####
## --------------------------------------------------------------------------#
#' @describeIn HiCDOCDataSet access and print the differences
#' @export
#' @return For \code{differences}: a tibble of the differences found by
#' \code{detectCompartments}
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


#### - concordances ------------------------------------------------------####
## --------------------------------------------------------------------------#
#' @describeIn HiCDOCDataSet access and print the concordances
#' @export
#' @return For \code{concordances}: a tibble of the concordances found by
#' \code{detectCompartments}
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


#### - compartments ------------------------------------------------------####
## --------------------------------------------------------------------------#
#' @describeIn HiCDOCDataSet access and print the compartments
#' @export
#' @return For \code{compartments}: a tibble of the compartments found by
#' \code{detectCompartments}
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

#### - centroids --------------------------------------------------------#####
## --------------------------------------------------------------------------#
#' @describeIn HiCDOCDataSet access and print the centroids
#' @export
#' @return For \code{centroids}: a tibble of the centroids found by
#' \code{detectCompartments}
setMethod("centroids", "HiCDOCDataSet", function(object) object@centroids)


#### - show --------------------------------------------------------------####
## --------------------------------------------------------------------------#
#' @describeIn HiCDOCDataSet access and print the object
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


#### - parameters --------------------------------------------------------####
## --------------------------------------------------------------------------#
#' @details 
#' The \code{parameters} slot holds the parameter values
#' used in an experiment as a named \code{list}. 
#' 
#' Default values exist for parameters, but these can also be supplied 
#' as input values using the assignment function \code{parameters<-}
#' or by changing the default parameters in the functions using the 
#' corresponding parameter, see below. 
#'
#' \subsection{Parameters in a HiCDOC experiment:}{
#'    \describe{
#'        \item{\code{minLengthChr}}{The minimum chromosome
#'                size (in number of bins), to be kept by the function
#'                \code{\link{filterSmallChromosomes}}. Default to 100.}
#'        \item{\code{weakPosThreshold}}{To be kept by the function
#'                \code{\link{filterWeakPositions}}, the bins (positions)
#'                of a chromosome must have a mean value greater than
#'                \code{weakPosThreshold}, on all replicates and all
#'                conditions. The mean is computed on the row of the
#'                reconstructed full interaction matrix for 1 chromosome,
#'                1 condition and 1 replicate. Default to 0.}
#'        \item{\code{sparseThreshold}}{To be kept by the function
#'                \code{\link{filterSparseChromosomes}}, the sparsity 
#'                (percentage of empty cells) of the interactions matrix 
#'                must be lower than the threshold, on all replicates and 
#'                all conditions. Default to 0.95.}
#'       \item{\code{sampleSize}}{The number of bins used when sampling
#'                all the bins in \code{\link{normalizeDistanceEffect}}.
#'                Default to 20000.}
#'       \item{\code{kMeansIterations}}{The maximum number of 2-means
#'                iterations, used in \code{\link{detectCompartments}} 
#'                function.}
#'       \item{\code{kMeansDelta}}{The stop criterion of convergence of
#'                the 2-means method, used in \code{\link{detectCompartments}} 
#'                function.}
#'       \item{\code{kMeansRestarts}}{The maximum number of restarts
#'                for the 2-means, used in \code{\link{detectCompartments}} 
#'                function.}
#'    }
#' }
#' @docType methods
#' @name parameters
#' @title Print or change the parameters of an HiCDOCDataSet object
#' @description parameters Access and print the parameters.
#' @aliases parameters HiCDOCDefaultParameters parameters<-
#' @export
#' @usage 
#' ## S4 replacement method for signature 'HiCDOCDataSet'
#' parameters(object)
#' @param object a HiCDOCDataSet object
#' @examples
#' object <- HiCDOCExample()
#' 
#' # Default parameters - when object is constructed
#' parameters(object) 
#' 
#' # Changing a parameter
#' parameters(object) <- list("minLengthChr" = 50)
#' parameters(object)
NULL

#' @describeIn HiCDOCDataSet access and print the parameters
#' @return For \code{paramters}: a named list of the parameters use in the 
#' the filters functions  and \code{detectCompartments}
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

# @describeIn parameters Change the values of parameters.
#' @rdname parameters
#' @param value a named list with the new parameters values. The names should
#' be in the \code{HiCDOCDefaultParameters} names.
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
