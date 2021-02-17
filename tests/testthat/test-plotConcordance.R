test_that("plotConcordance behaves as expected", {
    object <- HiCDOCExample()
    expect_error(pp <- plotConcordance(object), 
                 "Please run 'detectCompartments' first.")
    set.seed(3215)
    object <- detectCompartments(object)
    expect_error(plotConcordance(object), 'chromosomeId" est manquant')
    expect_error(plotConcordance(object, 3), 'Unknown')
    
    pp <- plotConcordance(object, 1)
    expect_is(pp, "ggplot")
    expect_identical(pp$labels, 
                     list("caption" = "No change is significant (pAdj < 5 %)",
                          "x" = "position",
                          "y" = "concordance",
                          "colour" = "replicate",
                          "yintercept" = "yintercept"))
    expect_is(pp$layers[[1]]$geom, "GeomLine")
    # No error when printed
    expect_error(print(pp), NA)
})
