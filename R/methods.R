#' @title
#' Methods to manipulate a \code{\link{HiCDOCDataSet}}.
#'
#' @docType methods
#'
#' @name
#' HiCDOCDataSet-methods
#'
#' @aliases
#' chromosomes conditions replicates resolution interactions compartments
#' concordances differences show
#'
#' @description
#' Retrieve information and results from a \code{\link{HiCDOCDataSet}}.
#'
#' @usage
#' chromosomes(object)
#' positions(object)
#' conditions(object)
#' replicates(object)
#' resolution(object)
#' interactions(object)
#' compartments(object)
#' concordances(object)
#' differences(object)
#' show(object)
#' @param object
#' a HiCDOCDataSet object
#' @examples
#' object <- HiCDOCDataSetExample()
#' chromosomes(object)
#' conditions(object)
#' object <- HiCDOC(object)
#' differences(object)
#' @return
#' A character vector (for \code{chromosomes}, \code{conditions},
#' \code{replicates}), an integer(for \code{resolution}), a tibble
#' (for \code{interactions}), or a GRanges object (for \code{compartments},
#' \code{concordances}, \code{differences}).
NULL

#' @describeIn HiCDOCDataSet-methods
#' Retrieves the vector of chromosome names.
#' @usage
#' NULL
#' @export
setMethod("chromosomes", "HiCDOCDataSet", function(object) object@chromosomes)

#' @describeIn HiCDOCDataSet-methods
#' Retrieves the genomic positions corresponding to bins for each chromosome.
#' @usage
#' NULL
#' @export
setMethod("positions", "HiCDOCDataSet", function(object) object@positions)

#' @describeIn HiCDOCDataSet-methods
#' Retrieves the vector of condition names.
#' @usage
#' NULL
#' @export
setMethod("conditions", "HiCDOCDataSet", function(object) object@conditions)

#' @describeIn HiCDOCDataSet-methods
#' Retrieves the vector of replicate names.
#' @usage
#' NULL
#' @export
setMethod("replicates", "HiCDOCDataSet", function(object) object@replicates)

#' @describeIn HiCDOCDataSet-methods
#' Retrieves the resolution (span of each position in number of bases).
#' @usage
#' NULL
#' @export
setMethod("resolution", "HiCDOCDataSet", function(object) object@binSize)

#' @describeIn HiCDOCDataSet-methods
#' Retrieves a tibble of the interactions.
#' @usage
#' NULL
#' @export
setMethod("interactions", "HiCDOCDataSet", function(object) {

    if (is.null(object@interactions)) return(NULL)

    interactions <-
        object@interactions %>%
        dplyr::left_join(
            object@positions %>% dplyr::select(
                chromosome,
                bin.1 = bin,
                position.1 = start,
            ),
            by = c("chromosome", "bin.1")
        ) %>%
        dplyr::left_join(
            object@positions %>% dplyr::select(
                chromosome,
                bin.2 = bin,
                position.2 = start,
            ),
            by = c("chromosome", "bin.2")
        ) %>%
        dplyr::select(
            chromosome,
            position.1,
            position.2,
            condition,
            replicate,
            interaction
        )

    return(interactions)
})

#' @describeIn HiCDOCDataSet-methods
#' Retrieves a \code{GenomicRange} of the compartment of every position
#' in every condition.
#' @usage
#' NULL
#' @export
setMethod("compartments", "HiCDOCDataSet", function(object) {

    if (is.null(object@compartments)) return(NULL)

    compartments <-
        object@compartments %>%
        dplyr::left_join(object@positions, by = c("chromosome", "bin")) %>%
        dplyr::select(
            chromosome,
            start,
            end,
            condition,
            compartment
        ) %>%
        dplyr::arrange(chromosome, condition, start, end) %>%
        dplyr::mutate(
            consecutive = (
                start - dplyr::lag(end) == 1 &
                dplyr::lag(chromosome) == chromosome &
                dplyr::lag(condition) == condition &
                dplyr::lag(compartment) == compartment
            ),
            consecutive = tidyr::replace_na(consecutive, TRUE),
            switching = dplyr::if_else(consecutive, 0, 1),
            group = cumsum(switching)
        ) %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(start = min(start), end = max(end)) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select(-consecutive, -switching, -group) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    return(compartments)
})

#' @describeIn HiCDOCDataSet-methods
#' Retrieves a \code{GenomicRange} of the significant compartment differences
#' between conditions, and their p-values.
#' @usage
#' NULL
#' @export
setMethod("differences", "HiCDOCDataSet", function(object, threshold = NULL) {

    if (is.null(object@differences)) return(NULL)

    if (
        !is.null(threshold) &&
        (!is.numeric(threshold) || length(threshold) > 1)
    ) {
        stop(
            "'threshold' must be a number.",
            call. = FALSE
        )
    }

    differences <-
        object@differences %>%
        dplyr::left_join(
            object@positions,
            by = c("chromosome", "bin")
        ) %>%
        dplyr::mutate(
            significance = dplyr::case_when(
                pvalue.adjusted <= 0.0001 ~ "****",
                pvalue.adjusted <= 0.001 ~ "***",
                pvalue.adjusted <= 0.01 ~ "**",
                pvalue.adjusted <= 0.05 ~ "*",
                TRUE ~ ""
            )
        ) %>%
        dplyr::select(
            chromosome,
            start,
            end,
            condition.1,
            condition.2,
            pvalue,
            pvalue.adjusted,
            direction,
            significance
        )

    if (!is.null(threshold)) {
        differences %<>% dplyr::filter(pvalue.adjusted <= threshold)
    }

    if (nrow(differences) == 0) {
        message(
            "No differences found",
            if (!is.null(threshold))
            paste0(" with adjusted p-value <= ", threshold),
            "."
        )
        return(GenomicRanges::GRanges())
    }

    genomicRange <-
        GenomicRanges::makeGRangesFromDataFrame(
            differences,
            keep.extra.columns = TRUE
        )

    return(genomicRange)
})

#' @describeIn HiCDOCDataSet-methods
#' Retrieves a \code{GenomicRange} of the concordance (confidence in assigned
#' compartment) of every position in every replicate.
#' @usage
#' NULL
#' @export
setMethod("concordances", "HiCDOCDataSet", function(object) {

    if (is.null(object@concordances)) return(NULL)

    concordances <-
        object@concordances %>%
        dplyr::left_join(object@positions, by = c("chromosome", "bin")) %>%
        dplyr::select(
            chromosome,
            start,
            end,
            condition,
            replicate,
            compartment,
            concordance
        ) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    return(concordances)
})

#' @describeIn HiCDOCDataSet-methods
#' Describes the object and its methods.
#' @usage
#' NULL
#' @export
setMethod("show", "HiCDOCDataSet", function(object) {
    cat(
        "Object of class 'HiCDOCDataSet'\n\n",
        "- Chromosomes:\n  ",
        if (is.null(object@chromosomes) || length(object@chromosomes) == 0)
        "None"
        else
        paste(object@chromosomes, collapse = ", "),
        "\n\n",
        "- Replicates:\n",
        if (is.null(object@replicates) || length(object@replicates) == 0)
        "  None\n"
        else
        paste0(
            "  condition ",
            object@conditions,
            ", replicate ",
            object@replicates,
            "\n"
        ),
        "\n",
        "- Resolution:\n  ",
        if (is.null(object@binSize))
        "None"
        else
        object@binSize,
        "\n\n",
        "- Methods:\n",
        "  chromosomes(object)\n",
        "  conditions(object)\n",
        "  replicates(object)\n",
        "  resolution(object)\n",
        "  interactions(object)\n",
        "  compartments(object)\n",
        "  differences(object)\n",
        "  concordances(object)\n",
        "  parameters(object)\n",
        "  parameters(object) <- list()\n\n",
        sep = ""
    )
})

#' @title
#' Access the parameters of a \code{\link{HiCDOCDataSet}}.
#'
#' @docType methods
#'
#' @name
#' HiCDOCDataSet-parameters
#'
#' @aliases
#' parameters parameters<- defaultHiCDOCParameters
#'
#' @description
#' Retrieves or sets parameters used for filtering, normalization, and
#' prediciton of compartments.
#'
#' @details
#' A \code{\link{HiCDOCDataSet}}'s parameters are automatically set to default
#' values retrieved from \code{\link{defaultHiCDOCParameters}}. They are
#' accessed by filtering, normalization, and compartment detection functions.
#' If those functions are called with custom arguments, the object's
#' parameters are updated to record the actual parameters used. If the
#' object's parameters are customized before calling the functions, the
#' custom parameters will be used.
#'
#' \subsection{All parameters are listed here:}{
#'     \describe{
#'         \item{\code{smallChromosomeThreshold}}{
#'             The minimum length (number of positions) for a chromosome to be
#'             kept when filtering with \code{\link{filterSmallChromosomes}}.
#'             Defaults to
#'             \code{defaultHiCDOCParameters$smallChromosomeThreshold} = 100.
#'         }
#'         \item{\code{sparseReplicateThreshold}}{
#'             The minimum percentage of non-zero interactions for a chromosome
#'             replicate to be kept when filtering with
#'             \code{\link{filterSparseReplicates}}. If a chromosome replicate's
#'             percentage of non-zero interactions is lower than this value, it
#'             is removed. Defaults to
#'             \code{defaultHiCDOCParameters$smallChromosomeThreshold} = 0.05.
#'         }
#'         \item{\code{weakPositionThreshold}}{
#'             The minimum average interaction for a position to be kept when
#'             filtering with \code{\link{filterWeakPositions}}. If a position's
#'             average interaction with the entire chromosome is lower than this
#'             value in any of the replicates, it is removed from all replicates
#'             and conditions. Defaults to
#'             \code{defaultHiCDOCParameters$smallChromosomeThreshold} = 1.
#'         }
#'         \item{\code{loessSampleSize}}{
#'             The number of positions used as a sample to estimate the effect
#'             of distance on proportion of interactions when normalizing with
#'             \code{\link{normalizeDistanceEffect}} Defaults to
#'             \code{defaultHiCDOCParameters$loessSampleSize} = 20000.
#'         }
#'         \item{\code{kMeansDelta}}{
#'             The convergence stop criterion for the clustering when detecting
#'             compartments with \code{\link{detectCompartments}}. When the
#'             centroids' distances between two iterations is lower than this
#'             value, the clustering stops. Defaults to
#'             \code{defaultHiCDOCParameters$kMeansDelta} = 0.0001.
#'         }
#'         \item{\code{kMeansIterations}}{
#'             The maximum number of iterations during clustering when detecting
#'             compartments with \code{\link{detectCompartments}}. Defaults to
#'             \code{defaultHiCDOCParameters$kMeansIterations} = 50.
#'         }
#'         \item{\code{kMeansRestarts}}{
#'             The amount of times the clustering is restarted when detecting
#'             compartments with \code{\link{detectCompartments}}. For each
#'             restart, the clustering iterates until convergence or reaching
#'             the maximum number of iterations. The clustering that minimizes
#'             inner-cluster variance is selected. Defaults to
#'             \code{defaultHiCDOCParameters$kMeansRestarts} = 20.
#'         }
#'     }
#' }
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @examples
#' object <- HiCDOCDataSetExample()
#'
#' # Retrieve parameters
#' parameters(object)
#'
#' # Set parameters
#' parameters(object) <- list("smallChromosomeThreshold" = 50)
#' parameters(object) <- list(
#'     "weakPositionThreshold" = 10,
#'     "kMeansRestarts" = 30
#' )
#'
#' @usage
#' parameters(object)
#' parameters(object) <- list()
NULL

#' @describeIn HiCDOCDataSet-parameters
#' Retrieves the parameters used for filtering, normalization, and prediction of
#' compartments. See
#' \code{\link{filterSmallChromosomes}},
#' \code{\link{filterSparseReplicates}},
#' \code{\link{filterWeakPositions}},
#' \code{\link{normalizeDistanceEffect}}, and
#' \code{\link{detectCompartments}},
#' for details on how these parameters are used.
#' @usage
#' NULL
#' @export
setMethod("parameters", "HiCDOCDataSet", function(object) object@parameters)

#' @describeIn HiCDOCDataSet-parameters
#' Sets the parameters used for filtering, normalization, and prediciton of
#' compartments. See
#' \code{\link{filterSmallChromosomes}},
#' \code{\link{filterSparseReplicates}},
#' \code{\link{filterWeakPositions}},
#' \code{\link{normalizeDistanceEffect}}, and
#' \code{\link{detectCompartments}},
#' for details on how these parameters are used.
#' @usage
#' NULL
#' @export
setReplaceMethod("parameters", "HiCDOCDataSet", function(object, value) {

    defaultParameterNames <- names(defaultHiCDOCParameters)

    if (!is(value, "list")) {
        stop(
            "'parameters' must be a named list.\n",
            "No parameters were updated. ",
            "See 'help(parameters)' for details.",
            call. = FALSE
        )
    }

    parameterNames <- names(value)

    duplicatedParameterNames <-
        unique(parameterNames[duplicated(parameterNames)])

    if (length(duplicatedParameterNames) > 0) {
        stop(
            "Duplicate parameter",
            if (length(duplicatedParameterNames) != 1) "s",
            " provided: ",
            paste(duplicatedParameterNames, collapse = ", "),
            "\nNo parameters were updated. ",
            "See 'help(parameters)' for details.",
            call. = FALSE
        )
    }

    invalidParameterNames <-
        parameterNames[!(parameterNames %in% defaultParameterNames)]

    if (length(invalidParameterNames) > 0) {
        stop(
            "Invalid parameter",
            if (length(invalidParameterNames) != 1) "s",
            " provided: ",
            paste(invalidParameterNames, collapse = ", "),
            "\nNo parameters were updated. ",
            "See 'help(parameters)' for details.",
            call. = FALSE
        )
    }

    object@parameters[parameterNames] <- value

    return(object)
})
