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
#' @import rhdf5
#' @import devtools
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import readr
#' @import reshape2
#' @import HiCcompare
#' @importFrom multiHiCcompare make_hicexp cyclic_loess hic_table
#' @import ggplot2
#' @import hexbin
#' @import gtools
#' @import BiocParallel
#' @importFrom BiocManager version
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom GenomicRanges makeGRangesFromDataFrame makeGRangesListFromDataFrame start end distance
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggpubr as_ggplot text_grob
#' @importFrom magrittr %<>%
#' @importFrom purrr map map2 flatten_int flatten_chr map_dfr pmap_dfr
#' @importFrom plyr compact
#' @importFrom stats D ecdf loess loess.control median optimize p.adjust prcomp predict quantile setNames update
#' @importFrom utils modifyList read.table
#' @importFrom gtable gtable_filter
#' @importFrom ggExtra ggMarginal
#' @importFrom parallel mcmapply detectCores
#' @useDynLib HiCDOC
#'
#' @keywords package
NULL
