###############################################################################
#
# HiCDOC-package organization of R files
#
# AllClasses .... class definitions and object constructors
# AllGenerics ... the generics defined in srnadiff-package
# methods ....... the S4 methods (accessors for slots and replace methods)
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
#' @title Find compartments with significantly differential binding using
#' replicated Hi-C data
#' @import methods
#' @import ggplot2
#' @importFrom GenomicRanges start end distance
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom stats D loess median quantile
"_PACKAGE"
