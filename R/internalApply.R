# code inspired from multiHiCcompare package (fonction smartApply)
.internalLapply <- function(parallel, ...) {
    if (parallel) {
        if (!is.null(BiocParallel::bpparam()$RNGseed)) {
            warning(
                "The use of RNGseed may not be ensured ",
                "See ?detectCompartments for more details",
                call. = FALSE
            )
        }
        BiocParallel::bplapply(..., BPPARAM = BiocParallel::bpparam())
    } else {
        pbapply::pblapply(...)
    }
}
