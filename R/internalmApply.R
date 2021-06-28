# code inspired from multiHiCcompare package (fonction smartApply)
.internalmApply <- function(parallel, n, ...) {
    if (parallel) {
        if(!is.null(BiocParallel::bpparam()$RNGseed)) {
             warning("The use of RNGseed may not be ensured ",
             "See ?detectCompartments for more details",
            call. = FALSE)
        }
        BiocParallel::bpmapply(...,
                                SIMPLIFY = FALSE,
                               BPPARAM = BiocParallel::bpparam())
    } else {
        pbapply::pbmapply(..., SIMPLIFY = FALSE)
    }
}
