test_that("detectCompartments behaves as expected", {
    obj <- HiCDOCExample()
    # Detect Compartments
    set.seed(3215) # Test with 123 : no significant differences
    expect_message(obj <- detectCompartments(obj), "Done.")
    # Keep obj format
    expect_is(obj@interactions$chromosome, "factor")
    expect_is(obj@interactions$bin.1, "integer")
    expect_is(obj@interactions$bin.2, "integer")
    expect_is(obj@interactions$condition, "factor")
    expect_is(obj@interactions$replicate, "factor")
    expect_is(obj@interactions$value, "numeric")
    # Create new objects in correct format
    expect_is(obj@distances, "tbl_df")
    expect_is(obj@diagonalRatios, "tbl_df")
    expect_is(obj@compartments, "tbl_df")
    expect_is(obj@concordances, "tbl_df")
    expect_is(obj@differences, "tbl_df")
    expect_is(obj@centroids, "tbl_df")
    
    # Differences
    expect_equal(nrow(obj@differences), 12)
    expect_equal(nrow(obj@differences %>% dplyr::filter(padj<=0.05)), 6)
    expect_is(obj@differences$chromosome, "factor")
    expect_is(obj@differences$condition.1, "factor")
    expect_is(obj@differences$condition.2, "factor")
    
    # Centroids
    expect_equal(nrow(obj@centroids), 8)
    expect_equal(lapply(obj@centroids$centroid, length), 
                 list(127, 127, 127, 127, 112, 112, 112, 112))
    expect_equal(lapply(obj@centroids$centroid, mean), 
                 list(316.5268, 282.0699, 182.8555, 232.813, 
                      361.9592, 448.6829, 252.8538, 358.9238), 
                 tolerance = 1e-04)
    expect_is(obj@centroids$chromosome, "factor")
    expect_is(obj@centroids$condition, "factor")
    expect_is(obj@centroids$compartment, "factor")
    expect_is(obj@centroids$centroid, "list")
    
    # Compartments
    expect_equal(nrow(obj@compartments), 478)
    expect_equal(nrow(obj@compartments %>% 
                          dplyr::filter(compartment=="B")), 162)
    expect_is(obj@compartments$chromosome, "factor")
    expect_is(obj@compartments$condition, "factor")
    expect_is(obj@compartments$compartment, "factor")
    expect_is(obj@compartments$bin, "integer")
    
    # Concordance
    expect_equal(nrow(obj@concordances), 1434)
    expect_equal(nrow(obj@concordances %>% 
                          dplyr::filter(compartment==1)), 798)
    expect_equal(mean(obj@concordances$concordance), 
                 -0.004183894, tolerance = 1e-05)
    expect_is(obj@concordances$chromosome, "factor")
    expect_is(obj@concordances$bin, "integer")
    expect_is(obj@concordances$condition, "factor")
    expect_is(obj@concordances$replicate, "factor")
    expect_is(obj@concordances$compartment, "numeric")
    expect_is(obj@concordances$concordance, "numeric")
    
    # Distances
    expect_equal(nrow(obj@distances), 2868)
    expect_equal(mean(obj@distances$distance), 7274.914, tolerance = 1e-04)
    expect_is(obj@distances$chromosome, "factor")
    expect_is(obj@distances$bin, "integer")
    expect_is(obj@distances$condition, "factor")
    expect_is(obj@distances$replicate, "factor")
    expect_is(obj@distances$compartment, "factor")
    expect_is(obj@distances$distance, "numeric")

    # DiagonalRatios
    expect_equal(nrow(obj@diagonalRatios), 1434)
    expect_equal(mean(obj@diagonalRatios$value), 4323.244, tolerance = 1e-04)
    expect_is(obj@diagonalRatios$chromosome, "factor")
    expect_is(obj@diagonalRatios$bin, "integer")
    expect_is(obj@diagonalRatios$condition, "factor")
    expect_is(obj@diagonalRatios$replicate, "factor")
    expect_is(obj@diagonalRatios$value, "numeric")
})

test_that("detectCompartments behaves as expected in parallel", {
    obj <- HiCDOCExample()
    # Detect Compartments
    multiParam <- BiocParallel::MulticoreParam(workers = 3, progressbar = TRUE)
    BiocParallel::register(multiParam, default = TRUE)
    expect_message(obj <- detectCompartments(obj, parallel = TRUE), 
                   "Done.")
    # Keep obj format
    expect_is(obj@interactions$chromosome, "factor")
    expect_is(obj@interactions$bin.1, "integer")
    expect_is(obj@interactions$bin.2, "integer")
    expect_is(obj@interactions$condition, "factor")
    expect_is(obj@interactions$replicate, "factor")
    expect_is(obj@interactions$value, "numeric")
    # Create new objects in correct format
    expect_is(obj@distances, "tbl_df")
    expect_is(obj@diagonalRatios, "tbl_df")
    expect_is(obj@compartments, "tbl_df")
    expect_is(obj@concordances, "tbl_df")
    expect_is(obj@differences, "tbl_df")
    expect_is(obj@centroids, "tbl_df")
    
    # Differences
    expect_is(obj@differences$chromosome, "factor")
    expect_is(obj@differences$condition.1, "factor")
    expect_is(obj@differences$condition.2, "factor")
    
    # Centroids
    expect_is(obj@centroids$chromosome, "factor")
    expect_is(obj@centroids$condition, "factor")
    expect_is(obj@centroids$compartment, "factor")
    expect_is(obj@centroids$centroid, "list")
    
    # Compartments
    expect_is(obj@compartments$chromosome, "factor")
    expect_is(obj@compartments$condition, "factor")
    expect_is(obj@compartments$compartment, "factor")
    expect_is(obj@compartments$bin, "integer")
    
    # Concordance
    expect_is(obj@concordances$chromosome, "factor")
    expect_is(obj@concordances$bin, "integer")
    expect_is(obj@concordances$condition, "factor")
    expect_is(obj@concordances$replicate, "factor")
    expect_is(obj@concordances$compartment, "numeric")
    expect_is(obj@concordances$concordance, "numeric")
    
    # Distances
    expect_is(obj@distances$chromosome, "factor")
    expect_is(obj@distances$bin, "integer")
    expect_is(obj@distances$condition, "factor")
    expect_is(obj@distances$replicate, "factor")
    expect_is(obj@distances$compartment, "factor")
    expect_is(obj@distances$distance, "numeric")
    
    # DiagonalRatios
    expect_is(obj@diagonalRatios$chromosome, "factor")
    expect_is(obj@diagonalRatios$bin, "integer")
    expect_is(obj@diagonalRatios$condition, "factor")
    expect_is(obj@diagonalRatios$replicate, "factor")
    expect_is(obj@diagonalRatios$value, "numeric")
})
