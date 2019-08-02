#' @export
normalizeKnightRuiz <- function(object) {
  input <- object@interactionMatrix %>%
    mutate(bin1 = `position 1` / object@binSize + 1) %>%
    mutate(bin2 = `position 2` / object@binSize + 1)

  outputTidy <- tibble()
  # outputTidy <- tibble(chromosome = double(),
  #                      replicate = character(),
  #                      values = double(),
  #                      bin1 = integer(),
  #                      bin2 = integer())
  for (chr in object@chromosomes) {
    inputChromosome <- filter(input, chromosome == chr) %>%
      spread(replicate, value)
    inputReplicate  <- tibble(bin1 = inputChromosome$bin1,
                              bin2 = inputChromosome$bin2,
                              data = 0)
    n <- max(inputChromosome$bin1, inputChromosome$bin2)
    message(paste0("Chromosome ", chr, ", of dim. ", n))
    for (replicate in object@replicates) {
      message(paste0("  Replicate ", replicate))
      inputReplicate$data <- inputChromosome[[replicate]]
      mat <- matrix(0, nrow = n, ncol = n)
      tmp <- as.matrix(inputReplicate)
      mat[ tmp[, 1:2] ] <- tmp[, 3]
      mat[ rev(tmp[, 1:2]) ] <- tmp[, 3]
      nullRows <- which(colSums(mat) == 0)
      mat[nullRows, nullRows] <- 1
      matKR <- KRnorm(mat)
      matKR[nullRows, ] <- 0
      matKR[, nullRows] <- 0
      vecKR <- as.vector(t(matKR))
      vecKR[is.na(vecKR)] <- 0
      tmpOutput <- tibble(chromosome = chr,
                          bin1 = rep(seq(n), each = n),
                          bin2 = rep(seq(n), times = n),
                          replicate = replicate,
                          value = vecKR)
      outputTidy %<>% bind_rows(tmpOutput)
    }
  }
  outputTidy %<>%
    filter(bin1 <= bin2) %>%
    filter(value != 0.0) %>%
    mutate(`position 1` = (`bin1` - 1) * object@binSize,
                         `position 2` = (`bin2` - 1) * object@binSize) %>%
    select(-c(`bin1`, `bin2`))

  object@interactionMatrix <- outputTidy

  return(object)
}
