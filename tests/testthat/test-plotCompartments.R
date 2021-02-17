test_that("plotCompartments behaves as expected", {
    object <- HiCDOCExample()
    expect_error(pp <- plotCompartments(object), 
                 "Please run 'detectCompartments' first.")
    set.seed(3215)
    object <- detectCompartments(object)
    expect_error(plotCompartments(object), 'chromosomeId" est manquant')
    expect_error(plotCompartments(object, 3), 'Unknown')
    
    pp <- plotCompartments(object, 1)
    expect_is(pp, "ggplot")
    expect_identical(pp$labels, list("x" = "position",
                                     "fill" = "compartment",
                                     "y" = "count",
                                     "weight" = "weight"))
    expect_is(pp$layers[[1]]$geom, "GeomBar")
    # No error when printed
    expect_error(print(pp), NA)
})
