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
# HiCDOCExp
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
#'
#' @docType package
#' @name HiCDOC
#' @aliases HiCDOC-package
#'
#' @author Cyril Kurylo and Matthias Zytnicki and Sylvain Foissac
#'
#' @import methods
#' @import rhdf5
#' @import devtools
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import readr
#' @import purrr
#' @import reshape2
#' @import HiCcompare
#' @importFrom multiHiCcompare make_hicexp cyclic_loess hic_table
#' @import ggplot2
#' @import ggExtra
#' @import hexbin
#' @import gtools
#' @import GenomicRanges
#' @import BiocParallel
#' @importFrom BiocManager version
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom GridExtra arrangeGrob
#' @importFrom magrittr %<>%
#' @useDynLib HiCDOC
#'
#' @keywords package
NULL
