data(HiCDOCDataSetExample)

test_that("detectCompartments behaves as expected", {
    # Detect Compartments
    set.seed(3215) # Test with 123 : no significant differences
    expect_message(
        object <- detectCompartments(HiCDOCDataSetExample, parallel = FALSE),
        "Detecting significant differences."
    )
    # Keep object format
    expect_is(object@interactions$chromosome, "factor")
    expect_is(object@interactions$bin.1, "integer")
    expect_is(object@interactions$bin.2, "integer")
    expect_is(object@interactions$condition, "factor")
    expect_is(object@interactions$replicate, "factor")
    expect_is(object@interactions$interaction, "numeric")
    # Create new objects in correct format
    expect_is(object@distances, "tbl_df")
    expect_is(object@selfInteractionRatios, "tbl_df")
    expect_is(object@compartments, "tbl_df")
    expect_is(object@concordances, "tbl_df")
    expect_is(object@differences, "tbl_df")
    expect_is(object@centroids, "tbl_df")

    # Differences
    expect_equal(nrow(object@differences), 52)
    expect_equal(
        nrow(object@differences %>% dplyr::filter(pvalue.adjusted <= 0.05)),
        0
    )
    expect_is(object@differences$chromosome, "factor")
    expect_is(object@differences$condition.1, "factor")
    expect_is(object@differences$condition.2, "factor")

    # Centroids
    expect_equal(nrow(object@centroids), 20)
    expect_equal(
        sapply(object@centroids$centroid, length),
        c(rep(80, 6), rep(120,6), rep(160, 6), rep(200,2))
    )
    expect_equal(
        sapply(object@centroids$centroid, mean),
        c(
            637.7590, 562.3542, 357.6060, 431.8505, 582.9829, 702.6890,
            631.2621, 572.7756, 571.6183, 514.5372, 914.0674, 1036.9587,
            688.8065, 688.7135, 446.8057, 457.7582, 931.4028, 857.0596,
            676.4497, 803.4541
        ),
        tolerance = 1e-04
    )
    expect_is(object@centroids$chromosome, "factor")
    expect_is(object@centroids$condition, "factor")
    expect_is(object@centroids$compartment, "factor")
    expect_is(object@centroids$centroid, "list")

    # Compartments
    expect_equal(nrow(object@compartments), 1280)
    expect_equal(
        nrow(object@compartments %>% dplyr::filter(compartment == "B")),
        745
    )
    expect_is(object@compartments$chromosome, "factor")
    expect_is(object@compartments$condition, "factor")
    expect_is(object@compartments$compartment, "factor")
    expect_is(object@compartments$bin, "integer")

    # Concordance
    expect_is(object@concordances$chromosome, "factor")
    expect_is(object@concordances$bin, "integer")
    expect_is(object@concordances$condition, "factor")
    expect_is(object@concordances$replicate, "factor")
    expect_is(object@concordances$compartment, "factor")
    expect_is(object@concordances$concordance, "numeric")
    expect_equal(nrow(object@concordances), 2720)
    expect_equal(
      nrow(object@concordances %>% dplyr::filter(compartment == "A")),
      1162
    )
    expect_equal(
      100 * mean(object@concordances$concordance),
      0.02854235,
      tolerance = 1e-05
    )

    # Distances
    expect_is(object@distances$chromosome, "factor")
    expect_is(object@distances$bin, "integer")
    expect_is(object@distances$condition, "factor")
    expect_is(object@distances$replicate, "factor")
    expect_is(object@distances$compartment, "factor")
    expect_is(object@distances$distance, "numeric")
    expect_equal(nrow(object@distances), 5440)
    expect_equal(mean(object@distances$distance), 4542.623, tolerance = 1e-04)

    # SelfInteractionRatios
    expect_equal(nrow(object@selfInteractionRatios), 2714)
    expect_equal(
      mean(object@selfInteractionRatios$ratio),
      361.8327,
      tolerance = 1e-04
    )
    expect_is(object@selfInteractionRatios$chromosome, "factor")
    expect_is(object@selfInteractionRatios$bin, "integer")
    expect_is(object@selfInteractionRatios$condition, "factor")
    expect_is(object@selfInteractionRatios$replicate, "factor")
    expect_is(object@selfInteractionRatios$ratio, "numeric")
})

test_that("detectCompartments behaves as expected in parallel", {
    # Detect Compartments
    multiParam <- BiocParallel::MulticoreParam(
        workers = 3,
        progressbar = TRUE,
        RNGseed = 12345
    )
    BiocParallel::register(multiParam, default = TRUE)
    expect_message(
        object <- detectCompartments(HiCDOCDataSetExample, parallel = TRUE),
        "Detecting significant differences."
    )
    # Keep object format
    expect_is(object@interactions$chromosome, "factor")
    expect_is(object@interactions$bin.1, "integer")
    expect_is(object@interactions$bin.2, "integer")
    expect_is(object@interactions$condition, "factor")
    expect_is(object@interactions$replicate, "factor")
    expect_is(object@interactions$interaction, "numeric")
    # Create new objects in correct format
    expect_is(object@distances, "tbl_df")
    expect_is(object@selfInteractionRatios, "tbl_df")
    expect_is(object@compartments, "tbl_df")
    expect_is(object@concordances, "tbl_df")
    expect_is(object@differences, "tbl_df")
    expect_is(object@centroids, "tbl_df")

    # Differences
    expect_equal(nrow(object@differences), 52)
    expect_equal(
      nrow(object@differences %>% dplyr::filter(pvalue.adjusted <= 0.05)),
      0
    )
    expect_is(object@differences$chromosome, "factor")
    expect_is(object@differences$condition.1, "factor")
    expect_is(object@differences$condition.2, "factor")

    # Centroids
    expect_equal(nrow(object@centroids), 20)
    expect_equal(
      sapply(object@centroids$centroid, length),
      c(rep(80, 6), rep(120,6), rep(160, 6), rep(200,2))
    )
    expect_equal(
      sapply(object@centroids$centroid, mean),
      c(
          637.7590, 562.3542, 431.8505, 357.6060, 702.6890, 582.9829,
          572.7756, 631.2621, 571.6183, 514.5372, 1036.9587, 914.0674,
          688.8065, 688.7135, 457.7582, 446.8057, 931.4028, 857.0596,
          803.4541, 676.4497
      ),
      tolerance = 1e-04
    )
    expect_is(object@centroids$chromosome, "factor")
    expect_is(object@centroids$condition, "factor")
    expect_is(object@centroids$compartment, "factor")
    expect_is(object@centroids$centroid, "list")

    # Compartments
    expect_equal(nrow(object@compartments), 1280)
    expect_equal(
      nrow(object@compartments %>% dplyr::filter(compartment == "B")),
      745
    )
    expect_is(object@compartments$chromosome, "factor")
    expect_is(object@compartments$condition, "factor")
    expect_is(object@compartments$compartment, "factor")
    expect_is(object@compartments$bin, "integer")

    # Concordance
    expect_is(object@concordances$chromosome, "factor")
    expect_is(object@concordances$bin, "integer")
    expect_is(object@concordances$condition, "factor")
    expect_is(object@concordances$replicate, "factor")
    expect_is(object@concordances$compartment, "factor")
    expect_is(object@concordances$concordance, "numeric")
    expect_equal(nrow(object@concordances), 2720)
    expect_equal(
      nrow(object@concordances %>% dplyr::filter(compartment == "A")),
      1162
    )
    expect_equal(
      100 * mean(object@concordances$concordance),
      0.02854235,
      tolerance = 1e-05
    )

    # Distances
    expect_is(object@distances$chromosome, "factor")
    expect_is(object@distances$bin, "integer")
    expect_is(object@distances$condition, "factor")
    expect_is(object@distances$replicate, "factor")
    expect_is(object@distances$compartment, "factor")
    expect_is(object@distances$distance, "numeric")
    expect_equal(nrow(object@distances), 5440)
    expect_equal(mean(object@distances$distance), 4542.623, tolerance = 1e-04)

    # SelfInteractionRatios
    expect_equal(nrow(object@selfInteractionRatios), 2714)
    expect_equal(
      mean(object@selfInteractionRatios$ratio),
      361.8327,
      tolerance = 1e-04
    )
    expect_is(object@selfInteractionRatios$chromosome, "factor")
    expect_is(object@selfInteractionRatios$bin, "integer")
    expect_is(object@selfInteractionRatios$condition, "factor")
    expect_is(object@selfInteractionRatios$replicate, "factor")
    expect_is(object@selfInteractionRatios$ratio, "numeric")
})
