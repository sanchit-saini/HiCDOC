# Get the positions where we observe a A->B or B->A compartment change
# The object should be a HiCDOCDataSet, with compartments called.
.getCompartmentBoundaries <- function (object, distance) {
    # Convert compartments to data.tables
    comp <- data.table::as.data.table(as(compartments(object), "data.frame"))
    # Compute ref sizes
    size <- comp[, .(size = max(end)), by = seqnames]
    # Remove unchecked information
    comp <- comp[centroid.check & PC1.check & assignment.check, ]
    # For each bin, get the compartment for each condition
    comp <- data.table::dcast(comp, seqnames + start + end ~ condition, value.var = "compartment")
    # Discard lines with NAs
    comp <- comp[stats::complete.cases(comp), ]
    # Compute number of predicted compartments per bins
    # (could probably implement something better)
    nComp <- apply(comp[, 4:ncol(comp)], 1, FUN = function(x) { length(unique(x)) })
    # Keep bins with only one predicted compartment
    comp <- comp[nComp == 1]
    # Keep ref, start, end, and compartment
    comp <- comp[, 1:4]
    colnames(comp)[4] <- "compartment"
    # Compare bin compartment with previous bin compartment
    comp <- comp[, prevComp := data.table::shift(compartment, 1), by = seqnames]
    comp <- comp[compartment != prevComp, ]
    # Add the size of the refs
    comp <- data.table::merge.data.table(comp, size, by = "seqnames")
    # Keep bins that are not too close to the ref ends
    comp <- comp[(start >= distance) & (end <= size - distance), ]
    return(comp)
}

# Replace a line/column with the previous one
# The line should be shifted in order to have the diagonal at the right place
.modifyMatrix <- function (cm, index) {
    prevIndex <- index - 1
    # Extract the line that should be used as template
    # I do not know why, but specifying namespace "S4Vectors" for as.matrix() is compulsory.
    newLine <- S4Vectors::as.matrix(cm[prevIndex, ])
    newLine <- c(newLine[[1]], newLine)
    newLine <- head(newLine, -1)
    # Paste template
    # I do not know why, but S4Vectors::as.matrix()<- does not work (not implemented).
    as.matrix(cm[index, ]) <- newLine
    as.matrix(cm[, index]) <- newLine
    #if (! isSymmetric(S4Vectors::as.matrix(cm))) stop ("Output matrix is not symmetric.")
    return(cm)
}

# Write the regions as BED format
.writeRegionsToFile <- function (regions, ref, outputFileName) {
    regions <- data.table::as.data.table(as(regions, "data.frame"))
    regions <- regions[seqnames == ref, ]
    # BED format is 0-based for start, and 1-based for end
    regions[, start := start - 1]
    regions <- regions[, .(seqnames, start, end, id)]
    data.table::fwrite(regions, outputFileName, sep = "\t", col.names = FALSE)
}

# Write one sample in ContactMatrix to HiC-Pro format
.writeMatrixToFile <- function (cm, outputFileName) {
    cm <- data.table::as.data.table(Matrix::mat2triplet(S4Vectors::as.matrix(cm)))
    cm <- cm[!is.na(x), ]
    # Make sure that the matrix is triangular
    cm <- cm[i >= j, ]
    data.table::fwrite(cm, outputFileName, sep = "\t", col.names = FALSE)
}

# Write the modified bin
.writeBinToFile <- function (ref, start, end, outputFileName) {
    f <- file(outputFileName)
    writeLines(paste(ref, start, end, sep = "\t"), f)
    close(f)
}

# Write all the (possibly modified) samples to HiC-Pro format
.writeToFile <- function (object, cms, ref, start, end, outputDirName) {
    if (! dir.exists(outputDirName)) {
        dir.create(outputDirName, recursive = TRUE)
    }
    .writeRegionsToFile(regions(object), ref, file.path(outputDirName, "coords.bed"))
    .writeBinToFile(ref, start, end, file.path(outputDirName, "bin.txt"))
    conditions <- sampleConditions(object)
    replicates <- sampleReplicates(object)
    for (i in seq_along(conditions)) {
        outputFileName <- file.path(outputDirName, paste0("cond_", conditions[[i]], "_rep_", replicates[[i]], ".pro"))
        .writeMatrixToFile(cms[[i]], outputFileName)
    }
}

modifyCompartment <- function (object, distance, outputDirName) {
    # Keep a copy
    objectRaw <- object
    # Call the compartments
    object <- HiCDOC(object)
    # Get all the compartment boundaries
    compartmentBoundaries <- .getCompartmentBoundaries(object, distance)
    if (nrow(compartmentBoundaries) == 0) {
        stop("There is no compartment boundary.")
    }
    # Select a random compartment boundary
    randomIndex <- trunc(runif(1, min = 0, max = nrow(compartmentBoundaries)))
    compartmentBoundary <- compartmentBoundaries[randomIndex, ]
    ref   <- as.character(compartmentBoundary$seqnames)
    start <- compartmentBoundary$start
    end   <- compartmentBoundary$end
    # Restrict the dataset to this chromosome
    selectedRefIndices <- (seqnames(anchors(objectRaw)$first) ==
                           seqnames(anchors(objectRaw)$second)) &
                          (seqnames(anchors(objectRaw)$first) == ref)
    # This only selects the right anchors, not the regions
    objectRaw <- objectRaw[selectedRefIndices]
    # This restricts the regions
    objectRaw <- InteractionSet::reduceRegions(objectRaw)
    # Transform all the samples to ContactMatrix
    allSampleIds <- seq_along(sampleConditions(objectRaw))
    allIndices <- seq_along(regions(objectRaw))
    cms <- lapply(allSampleIds,
                  function (sampleId) {
                      InteractionSet::inflate(objectRaw,
                                              rows = allIndices,
                                              columns = allIndices,
                                              sample = sampleId)
                  } )
    # Select the samples of the first condition
    modifiedCondition <- sampleConditions(objectRaw)[[1]]
    modifiedSampleIds <- which(sampleConditions(objectRaw) == modifiedCondition)
    # Find the index that corresponds to the selected bin
    selectedIndex <- which(GenomicRanges::start(regions(objectRaw)) == start)
    # Modify the boundary of the samples of the first condition
    modifiedCMs <- lapply(cms[modifiedSampleIds], .modifyMatrix, index = selectedIndex)
    cms[modifiedSampleIds] <- modifiedCMs
    # Write the results to files in HiC-Pro format
    .writeToFile(objectRaw, cms, ref, start, end, outputDirName)
}
