###############################################################################
#
# HiCDOC-package organization of R files
#
# AllClasses .... class definitions and object constructors
# AllGenerics ... the generics defined in srnadiff-package
# methods ....... the S4 methods (accessors for slots and replace methods)
# helper ........ computeNormFactors, computeLogFoldChange, cvgNormalization,
#                 reconcileRegions, checkParameters, IRlist2GR
# RcppExports ... the R wrappers for the C++ functions (auto)
#
#
# General outline of the internal function calls.
# Note: not all of these functions are exported.
#
# --------------------------------
# HiCDOCDataSet
# |- functions
# --------------------------------
#
# --------------------------------
# HiCDOC
# |- HiCDOCCore
#    |- page
#       |- function
# --------------------------------
#
###############################################################################
#' @docType package
#' @name HiCDOC
#' @title Find compartments with significantly differential binding using
#' replicated Hi-C data
#' @aliases HiCDOC-package
#' @import methods
#' @import ggplot2
#' @importFrom GenomicRanges start end distance
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom stats D loess median quantile
#' @useDynLib HiCDOC
"_PACKAGE"
