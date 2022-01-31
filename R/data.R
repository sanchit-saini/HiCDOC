#' @title
#' Example HiCDOCDataSet.
#'
#' @description
#' A S4 HiCDOCDataSet object with 4 chromosomes, 3 conditions and 3 replicates.
#'
#' @format
#' S4 HiCDOCDataSet object with the following characteristics:
#' \describe{
#'   \item{chromosomes}{4 chromosomes: W, X, Y, Z}
#'   \item{conditions}{3 conditions: 1, 2, 3}
#'   \item{replicates}{3 replicates: R1, R2, R3}
#'   \item{binSize}{A resolution of 137 bases}
#' }
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}.
"exampleHiCDOCDataSet"


#' @title
#' Example HiCDOCDataSet, filtered, normalized and with compartements detected.
#'
#' @description
#' A S4 HiCDOCDataSet object with 3 chromosomes, 3 conditions and 3 replicates.
#' Can be retrieved by running :
#' \code{set.seed(123); 
#' exampleHiCDOCDataSetProcessed <- HiCDOC(exampleHiCDOCDataSet)}
#'
#' @format
#' S4 HiCDOCDataSet object with the following characteristics:
#' \describe{
#'   \item{chromosomes}{4 chromosomes: X, Y, Z}
#'   \item{conditions}{3 conditions: 1, 2, 3}
#'   \item{replicates}{3 replicates: R1, R2, R3}
#'   \item{binSize}{A resolution of 137 bases}
#' }
#'
#' @return
#' A \code{\link{HiCDOCDataSet}}, already filtered and normalized.
"exampleHiCDOCDataSetProcessed"
