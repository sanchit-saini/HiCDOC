#' @title
#' A/B compartment detection and differential analysis
#' @docType package
#' @import methods
#' @import zlibbioc
#' @import ggplot2
#' @importFrom GenomicRanges start end distance GRanges match union seqnames
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom stats loess loess.control ecdf predict median update optimize prcomp p.adjust
#' @importFrom S4Vectors DataFrame mcols split runLength Rle
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom SummarizedExperiment assay colData
#' @importFrom gtools mixedsort
#' @importFrom pbapply pbmapply
#' @importFrom BiocParallel bpparam bpmapply bplapply
#' @importFrom BiocGenerics cbind
#' @importFrom rhdf5 h5read
#' @importFrom ggpubr as_ggplot
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggExtra ggMarginal
#' @import data.table
#' @import InteractionSet
#' @useDynLib HiCDOC
#' @aliases HiCDOC-package
"_PACKAGE"
