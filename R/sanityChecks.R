#' @title
#' Check centroid PCA
#'
#' @description
#' Check whether centroids are correctly placed on a PCA.
#'
#' @param chromosomeName
#' A chromosome name or index in \code{chromosomes(object)}.
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A \code{data.table}, with 3 elements: the name of the chromosome,
#' whether the centroids of each compartment have the same sign 
#' only works when the number of conditions is 2), and whether the
#' variance explained by the first axis is greater than a threshold
#' (here, 75%).
#'
#' @examples
#' data(exampleHiCDOCDataSetProcessed)
#' .checkPca("X", exampleHiCDOCDataSetProcessed)
#'
#' @export
.checkPca <- function(chromosomeName, object) {
    compartments <- as.character(unique(object@centroids$compartment))
    pcaData <- .computePca(object, chromosomeName)
    pca     <- pcaData$PCA
    propvar <- pcaData$propvar
    f <- function (c) {
        return(length(unique(sign(pca[compartment == c, ]$PC1))) == 1)
    }
    centroid <- all(vapply(compartments, f, FUN.VALUE = TRUE))
    pc1 <- (propvar[[1]] >= object@parameters$PC1CheckThreshold)
    return(data.table(chromosome     = chromosomeName,
                      centroid.check = centroid,
                      PC1.check      = pc1))
}


#' @title
#' Check compartment assignment
#'
#' @description
#' Check "A" compartments are different from "B" compartments.
#'
#' @param chromosome
#' A chromosome name or index in \code{chromosomes(object)}.
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A \code{data.table}, with 2 elements: the name of the chromosome,
#' whether the check passed (a Wilcoxon test, with p-value less than 5%).
#'
#' @examples
#' data(exampleHiCDOCDataSetProcessed)
#' .checkCompartmentAssignment("X", exampleHiCDOCDataSetProcessed)
#'
#' @export
.checkCompartmentAssignment <- function(chromosomeName, object) {
    compartments <- object@compartments[chromosome == chromosomeName, ]
    selfInteractionRatios <- object@selfInteractionRatios[chromosome == chromosomeName, ]
    compartments <- data.table::merge.data.table(compartments,
                                                 selfInteractionRatios,
                                                 by = c("index", "condition"))
    t            <- wilcox.test(ratio ~ compartment, data = compartments)
    return(data.table(chromosome       = chromosomeName,
                      assignment.check = (t$p.value <= 0.05)))
}


#' @title
#' Check whether compartments satisfy sanity checks.
#'
#' @description
#' Check centroid PCA and compartment assignments.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A \code{data.table}, with 3 elements: the name of the chromosome,
#' whether the centroids of each compartment have the same sign 
#' only works when the number of conditions is 2), and whether the
#' variance explained by the first axis is greater than a threshold
#' (here, 75%), and whether the check passed (a Wilcoxon test, with p-value
#' less than 5%).
#'
#' @examples
#' data(exampleHiCDOCDataSetProcessed)
#' .checkResults(exampleHiCDOCDataSetProcessed)
#'
#' @export
.checkResults <- function(object) {
    pcaChecks         <- lapply(chromosomes(object), .checkPca, object)
    pcaChecks         <- data.table::rbindlist(pcaChecks)
    compartmentChecks <- lapply(chromosomes(object),
                                .checkCompartmentAssignment, object)
    compartmentChecks <- data.table::rbindlist(compartmentChecks)
    checks            <- pcaChecks
    checks[, assignment.check := compartmentChecks$assignment.check]
    object@checks <- checks
    return(object)
}
