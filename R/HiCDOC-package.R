#' @title
#' A/B compartment detection and differential analysis
#' @docType package
#' @import methods
#' @import zlibbioc
#' @import ggplot2
#' @importFrom GenomicRanges start end distance GRanges
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom stats D loess median quantile
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom DescTools GCD
#' @import data.table
#' @import InteractionSet
#' @useDynLib HiCDOC
#' @aliases HiCDOC-package
"_PACKAGE"
