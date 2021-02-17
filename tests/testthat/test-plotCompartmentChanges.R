test_that("plotCompartmentChanges behaves as expected", {
    object <- HiCDOCExample()
    expect_error(pp <- plotCompartmentChanges(object), 
                 "Please run 'detectCompartments' first.")
    set.seed(3215)
    object <- detectCompartments(object)
    expect_error(plotCompartmentChanges(object), 'chromosomeId" est manquant')
    expect_error(plotCompartmentChanges(object, 3), 'Unknown')
    
    pp <- plotCompartmentChanges(object, 1)
    expect_is(pp, "ggplot")
    expect_is(pp$labels, "list")
    expect_equal(length(pp$labels), 0)
    expect_is(pp$layers[[1]]$geom, "GeomDrawGrob")
    # No error when printed
    expect_error(print(pp), NA)
})
