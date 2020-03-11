KR <- function(A, tol = 1e-6, delta = 0.1, Delta = 3) {

  # Remove any cols/rows of 0s
  zeros <- unique(which(colSums(A) == 0), which(rowSums(A) == 0))
  if (length(zeros) > 0) {
    A = A[-zeros, -zeros]
    message('Cols/Rows removed: ')
    message(paste(" ", zeros, sep = " "))
  }

  n <- nrow(A)
  e <- matrix(1, nrow = n, ncol = 1)
  x0 <- e
  # Inner stopping criterior
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
    while (rho_km1 > innertol) { # Inner iteration by CG
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

      # Update search direction efficiently
      w <- x * (A %*% (x * p)) + v * p
      alpha <- rho_km1 / drop(t(p) %*% w)
      ap <- alpha * p
      # Test distance to boundary of cone
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

  # Refill cols/rows of 0s
  if (length(zeros) > 0) {
    refilled <- matrix(0, nrow=n+length(zeros), ncol=n+length(zeros))
    refilled[-zeros, -zeros] <- result
    return (refilled)
  }

  return(result)
}

#' @export
normalizeKnightRuiz <- function(object) {

  input <- object@interactions %>% mutate(
    bin1 = `position 1` / object@binSize + 1,
    bin2 = `position 2` / object@binSize + 1
  )

  object@interactions <- tibble()

  for (chromosome in object@chromosomes) {

    message("Chromosome: ", chromosome)
    chromosomeInteractions <- input[input$chromosome == chromosome,]
    totalBins <- max(chromosomeInteractions$bin1, chromosomeInteractions$bin2)
    if (totalBins == -Inf) next

    for (conditionId in unique(object@conditions)) {

      replicates <- object@replicates[which(object@conditions == conditionId)]

      for (replicateId in replicates) {
        message("Replicate: ", conditionId, ".", replicateId)
        replicateInteractions <- chromosomeInteractions %>%
          filter(condition == conditionId & replicate == replicateId) %>%
          select(bin1, bin2, value)
        normalizedInteractions <- KR(
          sparseInteractionsToMatrix(replicateInteractions, totalBins)
        )
        object@interactions %<>% bind_rows(
          tibble(
            chromosome = chromosome,
            `position 1` = (rep(seq(totalBins), each = totalBins) - 1) * object@binSize,
            `position 2` = (rep(seq(totalBins), times = totalBins) - 1) * object@binSize,
            condition = conditionId,
            replicate = replicateId,
            value = as.vector(t(normalizedInteractions))
          ) %>% filter(`position 1` <= `position 2` & value != 0)
        )
      }
    }
  }

  object@interactions %<>% mutate(
    chromosome = factor(chromosome),
    condition = factor(condition),
    replicate = factor(replicate)
  )

  return(object)
}
