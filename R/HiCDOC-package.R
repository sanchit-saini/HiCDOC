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

#' Title...
#'
#' \code{HiCODC} is a package...
#' TODO
#'
#' @examples
#' #TODO
#'
#' @docType package
#' @name HiCDOC
#' @aliases HiCDOC-package
#'
#' @author Cyril Kurylo and Matthias Zytnicki and Sylvain Foissac and
#' Élise Maigné
#'
#' @import methods
# @import dplyr
# @import tidyr
# @import HiCcompare
#' @import ggplot2
# @importFrom Rcpp evalCpp sourceCpp
#' @importFrom GenomicRanges start end distance
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib HiCDOC
#'
#' @keywords package
NULL
