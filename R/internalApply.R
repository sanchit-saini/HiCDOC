# code inspired from multiHiCcompare package (fonction smartApply)
.internalApply <- function(parallel, ..., type="mapply") {
    if (parallel) {
        if(!is.null(BiocParallel::bpparam()$RNGseed)) {
             warning("The use of RNGseed may not be ensured ",
             "See ?detectCompartments for more details",
            call. = FALSE)
        }
        if(type == "mapply") {
            BiocParallel::bpmapply(...,
                                   SIMPLIFY = FALSE,
                                   BPPARAM = BiocParallel::bpparam())
        } else {
            BiocParallel::bplapply(...,
                                   BPPARAM = BiocParallel::bpparam())
        }
    } else {
        if(type == "mapply") {
            pbapply::pbmapply(..., SIMPLIFY = FALSE)
        } else {
            pbapply::pblapply(...)
        }
    }
}
