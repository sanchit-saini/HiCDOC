# From multiHiCcompare package
# perform cyclic loess on a table
.cloess <- function(tab,
                    iterations,
                    verbose,
                    span,
                    degree = 1,
                    loess.criterion = "gcv") {
    # make matrix of IFs
    IF_mat <- as.matrix(tab[,-c("chr", "region1", "region2", "D"),
                            with = FALSE])
    # make index matrix
    idx_mat <- IF_mat
    idx_mat[idx_mat != 0] <- 1
    # log the matrix
    IF_mat <- log2(IF_mat + 1)
    n <- ncol(IF_mat)
    # begin cyclic loess
    for (i in seq(from = 1, to = iterations)) {
        for (j in seq(from = 1, to = (n - 1))) {
            for (k in seq(from = (j + 1), to = n)) {
                M <- IF_mat[, k] - IF_mat[, j]
                D <- tab$D
                if (is.na(span)) {
                    l <- .loess.as(
                        x = D,
                        y = M,
                        degree = degree,
                        criterion = loess.criterion,
                        control = stats::loess.control(
                            surface = "interpolate",
                            statistics = "approximate",
                            trace.hat = "approximate"
                        )
                    )
                } else {
                    l <- .loess.as(
                        x = D,
                        y = M,
                        degree = degree,
                        user.span = span,
                        criterion = loess.criterion,
                        control = stats::loess.control(
                            surface = "interpolate",
                            statistics = "approximate",
                            trace.hat = "approximate"
                        )
                    )
                }
                # calculate gcv and AIC
                traceL <- l$trace.hat
                sigma2 <- sum(l$residuals ^ 2) / (l$n - 1)
                aicc <-
                    log(sigma2) + 1 + 2 * (2 * (traceL + 1)) / (l$n - traceL - 2)
                gcv <- l$n * sigma2 / (l$n - traceL) ^ 2
                # print the span picked by gcv
                if (verbose) {
                    message("Span for loess: ", l$pars$span)
                    message("GCV for loess: ", gcv)
                    message("AIC for loess: ", aicc)
                }
                # adjust IFs
                IF_mat[, j] <- IF_mat[, j] + l$fitted / 2
                IF_mat[, k] <- IF_mat[, k] - l$fitted / 2
            }
        }
    }
    # anti-log IFs
    IF_mat <- (2 ^ IF_mat) - 1
    # reset zeros
    IF_mat <- IF_mat * idx_mat
    # set negative values to 0
    IF_mat[IF_mat < 0] <- 0
    # fix any potential Infs or NaN's
    IF_mat[is.nan(IF_mat)] <- 0
    IF_mat[is.infinite(IF_mat)] <- 0
    # recombine table
    tab <- cbind(tab[, 1:4, with = FALSE], IF_mat)
    return(tab)
}


# From multiHiCcompare package
# loess with Automatic Smoothing Parameter Selection adjusted possible
# range of smoothing originally from fANCOVA package
.loess.as <-
    function(x,
             y,
             degree = 1,
             criterion = c("aicc", "gcv"),
             family = c("gaussian",
                        "symmetric"),
             user.span = NULL,
             plot = FALSE,
             ...) {
        criterion <- match.arg(criterion)
        family <- match.arg(family)
        x <- as.matrix(x)
        
        data.bind <- data.frame(x = x, y = y)
        if (ncol(x) == 1) {
            names(data.bind) <- c("x", "y")
        } else {
            names(data.bind) <- c("x1", "x2", "y")
        }
        
        opt.span <- function(model,
                             criterion = c("aicc", "gcv"),
                             span.range = c(0.01, 0.9)) {
            as.crit <- function(x) {
                span <- x$pars$span
                traceL <- x$trace.hat
                sigma2 <- sum(x$residuals ^ 2) / (x$n - 1)
                aicc <-
                    log(sigma2) + 1 + 2 * (2 * (traceL + 1)) / (x$n - traceL -
                                                                    2)
                gcv <- x$n * sigma2 / (x$n - traceL) ^ 2
                result <- list(span = span,
                               aicc = aicc,
                               gcv = gcv)
                return(result)
            }
            criterion <- match.arg(criterion)
            fn <- function(span) {
                mod <- stats::update(model, span = span)
                as.crit(mod)[[criterion]]
            }
            result <- optimize(fn, span.range)
            return(list(
                span = result$minimum,
                criterion = result$objective
            ))
        }
        
        if (ncol(x) == 1) {
            if (is.null(user.span)) {
                fit0 <-
                    stats::loess(
                        y ~ x,
                        degree = degree,
                        family = family,
                        data = data.bind,
                        ...
                    )
                span1 <- opt.span(fit0, criterion = criterion)$span
            } else {
                span1 <- user.span
            }
            fit <-
                stats::loess(
                    y ~ x,
                    degree = degree,
                    span = span1,
                    family = family,
                    data = data.bind,
                    ...
                )
        } else {
            if (is.null(user.span)) {
                fit0 <- stats::loess(y ~ x1 + x2,
                                     degree = degree,
                                     family = family,
                                     data.bind,
                                     ...)
                span1 <- opt.span(fit0, criterion = criterion)$span
            } else {
                span1 <- user.span
            }
            fit <-
                stats::loess(
                    y ~ x1 + x2,
                    degree = degree,
                    span = span1,
                    family = family,
                    data = data.bind,
                    ...
                )
        }
        if (plot) {
            if (ncol(x) == 1) {
                m <- 100
                x.new <- seq(min(x), max(x), length.out = m)
                fit.new <- stats::predict(fit, data.frame(x = x.new))
                plot(x,
                     y,
                     col = "lightgrey",
                     xlab = "x",
                     ylab = "m(x)",
                     ...)
                lines(x.new, fit.new, lwd = 1.5, ...)
            } else {
                m <- 50
                x1 <- seq(min(data.bind$x1), max(data.bind$x1), len = m)
                x2 <- seq(min(data.bind$x2), max(data.bind$x2), len = m)
                x.new <- expand.grid(x1 = x1, x2 = x2)
                fit.new <- matrix(predict(fit, x.new), m, m)
                persp(
                    x1,
                    x2,
                    fit.new,
                    theta = 40,
                    phi = 30,
                    ticktype = "detailed",
                    xlab = "x1",
                    ylab = "x2",
                    zlab = "y",
                    col = "lightblue",
                    expand = 0.6
                )
            }
        }
        return(fit)
    }
