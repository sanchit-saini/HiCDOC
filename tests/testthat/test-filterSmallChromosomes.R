object <- HiCDOCExample()

test_that("filterSmallChromosomes behave as expected if no filter", {
    # No filter on the example dataset
    expect_message(object <- filterSmallChromosomes(object), 
                   "Keeping only the chromosomes with 100 bins or more")
    expect_equal(length(object@chromosomes), 2)
    expect_equal(nrow(object@interactions), 86736)
})

test_that("filterSmallChromosomes behave as expected with filter", {
    # Filter on 1 chromosome
    expect_message(object <- filterSmallChromosomes(object, 115), 
                   "Keeping only the chromosomes with 115 bins or more")
    expect_equal(length(object@chromosomes), 1)
    expect_identical(object@chromosomes,  "17")
    expect_equal(nrow(object@interactions), 48768)
    expect_identical(object@parameters$minLengthChr, 115)
})
