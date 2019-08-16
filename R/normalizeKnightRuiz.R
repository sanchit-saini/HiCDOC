KR <- function(A, tol = 1e-6, delta = 0.1, Delta = 3) {
  n <- nrow(A)
  e <- matrix(1, nrow = n, ncol = 1)
  x0 <- e
  # inner stopping criterior
  g <- 0.9
  etamax <- 0.1
  eta <- etamax
  stop_tol <- tol * .5
  x <- x0
  rt <- tol^2
  v <- x * (A %*% x)
  rk <- 1 - v
  rho_km1 <- drop(t(rk) %*% rk)
  rout <- rho_km1
  rold <- rout
  MVP <- 0 # We'll count matrix vector products.
  i <- 0 # Outer iteration count.
  while (rout > rt) { # Outer iteration
    i <- i + 1
    k <- 0
    y <- e
    innertol <- max(c(eta^2 * rout, rt))
    while (rho_km1 > innertol) { #Inner iteration by CG
      k <- k + 1
      if (k == 1) {
        Z <- rk / v
        p <- Z
        rho_km1 <- drop(t(rk) %*% Z)
      }
      else {
        beta <- rho_km1 / rho_km2
        p <- Z + beta * p
      }

      # update search direction efficiently
      w <- x * (A %*% (x * p)) + v * p
      alpha <- rho_km1 / drop(t(p) %*% w)
      ap <- alpha * p
      # test distance to boundary of cone
      ynew <- y + ap
      if (min(ynew) <= delta) {
        if (delta == 0) break()
        ind <- which(ap < 0)
        gamma <- min((delta - y[ind]) / ap[ind])
        y <- y + gamma * ap
        break()
      }
      if (max(ynew) >= Delta) {
        ind <- which(ynew > Delta)
        gamma <- min((Delta - y[ind]) / ap[ind])
        y <- y + gamma * ap
        break()
      }
      y <- ynew
      rk <- rk - alpha * w
      rho_km2 <- rho_km1
      Z <- rk / v
      rho_km1 <- drop(t(rk) %*% Z)
    }
    x <- x * y
    v <- x * (A %*% x)
    rk <- 1 - v
    rho_km1 <- drop(t(rk) %*% rk)
    rout <- rho_km1
    MVP <- MVP + k + 1
    # Update inner iteration stopping criterion.
    rat <- rout / rold
    rold <- rout
    res_norm <- sqrt(rout)
    eta_o <- eta
    eta <- g * rat
    if (g * eta_o^2 > 0.1) {
      eta <- max(c(eta, g * eta_o^2))
    }
    eta <- max(c(min(c(eta, etamax)), stop_tol / res_norm))
  }

  result <- t(t(x[,1] * A) * x[,1])
  return(result)
}

#' @export
normalizeKnightRuiz <- function(object) {
  input <- object@interactionMatrix %>%
    mutate(bin1 = `position 1` / object@binSize + 1) %>%
    mutate(bin2 = `position 2` / object@binSize + 1)

  outputTidy <- tibble()
  for (chr in object@chromosomes) {
    inputChromosome <- filter(input, chromosome == chr) %>%
      spread(replicate, value)
    inputReplicate  <- tibble(bin1 = inputChromosome$bin1,
                              bin2 = inputChromosome$bin2,
                              data = 0)
    n <- max(inputChromosome$bin1, inputChromosome$bin2)
    message(paste0("Chromosome ", chr, ", of dim. ", n))
    for (conditionId in c(1, 2)) {
      for (replicateId in seq.int(object@nReplicatesPerCond[conditionId])) {
        replicate <- paste0("replicate ", conditionId, ".", replicateId)
        message(paste0("  Replicate ", replicate))
        inputReplicate$data <- inputChromosome[[replicate]]
        mat <- matrix(0, nrow = n, ncol = n)
        tmp <- as.matrix(inputReplicate)
        mat[ tmp[, 1:2] ] <- tmp[, 3]
        mat <- mat + t(mat) - diag(diag(mat))
        if (!isSymmetric(mat)) {
          stop("Matrix is not symmetric.")
        }
        nullRows <- which((colSums(mat) == 0) | (rowSums(mat) == 0))
        if (length(nullRows) > 0) {
          message(paste0("    ",
                         length(nullRows),
                         " rows/columns are empty."))
        }
        diag(mat)[nullRows] <- 1
        matKR <- KR(mat)
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
