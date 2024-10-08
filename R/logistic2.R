# @rdname logistic6_new
logistic2_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  # 2-parameter log-logistic curve is tricky because according to our
  # parameterization we have two options:
  #
  # x^eta / (x^eta + phi^eta)
  #
  # or
  #
  # 1 - x^eta / (x^eta + phi^eta)
  #
  # this duality creates some problem while initializing the model because we
  # do not know yet if the curve is increasing or decreasing.
  #
  # We will try our best to infer it from the data.

  stats <- suff_stats(x, y, w)
  m <- nrow(stats)

  # depending on the quality of data, the guessing can be hard or easy
  #
  # we use `lowess` as a non-linear non-parametric smoother to get an idea of
  # the trend.
  xx <- sqrt(stats[, 2]) * stats[, 1]
  yy <- sqrt(stats[, 2]) * stats[, 3]
  fit <- lowess(xx, yy, f = 0.9)

  a <- 0
  d <- 1
  if (fit$y[1] > fit$y[m]) {
    a <- 1
    d <- -1
  }

  start <- if (!is.null(start)) {
    if (length(start) != 2) {
      stop("'start' must be of length 2", call. = FALSE)
    }

    if (start[1] <= 0) {
      stop("parameter 'eta' cannot be negative nor zero", call. = FALSE)
    }

    c(a, d, log(start[1]), start[2])
  } else {
    c(a, d, NA_real_, NA_real_)
  }

  object <- structure(
    list(
      x = x,
      y = y,
      w = w,
      n = length(y),
      stats = stats,
      constrained = FALSE,
      start = start,
      max_iter = max_iter,
      m = m
    ),
    class = "logistic2"
  )

  if (!is.null(lower_bound) || !is.null(upper_bound)) {
    object$constrained <- TRUE

    if (is.null(lower_bound)) {
      lower_bound <- rep(-Inf, 2)
    } else {
      if (length(lower_bound) != 2) {
        stop("'lower_bound' must be of length 2", call. = FALSE)
      }

      lower_bound[1] <- if (lower_bound[1] > 0) {
        log(lower_bound[1])
      } else {
        -Inf
      }
    }

    if (is.null(upper_bound)) {
      upper_bound <- rep(Inf, 2)
    } else {
      if (length(upper_bound) != 2) {
        stop("'upper_bound' must be of length 2", call. = FALSE)
      }

      if (upper_bound[1] <= 0) {
        stop("'upper_bound[1]' cannot be negative nor zero.", call. = FALSE)
      }

      upper_bound[1] <- log(upper_bound[1])
    }

    object$lower_bound <- lower_bound
    object$upper_bound <- upper_bound
  }

  object
}

#' 2-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
#' are free to vary (therefore the name) while vector `c(alpha, delta)` is
#' constrained to be either `c(0, 1)` (monotonically increasing curve) or
#' `c(1, -1)` (monotonically decreasing curve).
#'
#' This function allows values other than `{0, 1, -1}` for `alpha` and `delta`
#' but will coerce them to their proper constraints.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(alpha, delta, eta, phi)`. `alpha` can only be equal to 0 or 1 while
#'   `delta` can only be equal to 1 or -1.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
#'
#' @export
logistic2_fn <- function(x, theta) {
  alpha <- 0
  delta <- 1
  if (theta[2] < 0) {
    alpha <- 1
    delta <- -1
  }

  eta <- theta[3]
  phi <- theta[4]

  alpha + delta / (1 + exp(-eta * (x - phi)))
}

#' @export
fn.logistic2 <- function(object, x, theta) {
  logistic2_fn(x, c(object$start[1:2], theta))
}

#' @export
fn.logistic2_fit <- function(object, x, theta) {
  # within a fit, parameter theta is known exactly
  alpha <- theta[1]
  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  alpha + delta / (1 + exp(-eta * (x - phi)))
}

#' 2-parameter logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
#' are free to vary (therefore the name) while vector `c(alpha, delta)` is
#' constrained to be either `c(0, 1)` (monotonically increasing curve) or
#' `c(1, -1)` (monotonically decreasing curve).
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the two parameters in the form
#'   `c(eta, phi)`.
#' @param delta value of delta parameter (either 1 or -1).
#'
#' @return Gradient or Hessian evaluated at the specified point.
#'
#' @export
logistic2_gradient <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  b <- exp(-eta * (x - phi))

  f <- 1 + b

  q <- (x - phi) * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  G <- matrix(1, nrow = k, ncol = 2)

  G[, 1] <- t
  G[, 2] <- u

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(G)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the gradient at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    G[is_nan] <- 0
  }

  sign(delta) * G
}

#' @export
gradient.logistic2_fit <- function(object, x) {
  theta <- object$coefficients
  logistic2_gradient(x, theta[3:4], theta[2])
}

#' @export
#'
#' @rdname logistic2_gradient
logistic2_hessian <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  b <- exp(-eta * (x - phi))

  f <- 1 + b

  q <- (x - phi) * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- q * t * (2 / f - 1 / b)
  H[, 2, 1] <- (1 / eta + f * (2 - f / b) * t) * u

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- (2 / f - 1 / b) * r * u

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(H)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the Hessian at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    H[is_nan] <- 0
  }

  sign(delta) * H
}

#' @export
#'
#' @rdname logistic2_gradient
logistic2_gradient_hessian <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  b <- exp(-eta * (x - phi))

  f <- 1 + b

  q <- (x - phi) * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  G <- matrix(1, nrow = k, ncol = 2)

  G[, 1] <- t
  G[, 2] <- u

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- q * t * (2 / f - 1 / b)
  H[, 2, 1] <- (1 / eta + f * (2 - f / b) * t) * u

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- (2 / f - 1 / b) * r * u

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(G)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the gradient at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    G[is_nan] <- 0
  }

  is_nan <- is.nan(H)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the Hessian at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    H[is_nan] <- 0
  }

  list(G = sign(delta) * G, H = sign(delta) * H)
}

#' 2-parameter logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
#' are free to vary (therefore the name) while vector `c(alpha, delta)` is
#' constrained to be either `c(0, 1)` (monotonically increasing curve) or
#' `c(1, -1)` (monotonically decreasing curve).
#'
#' This set of functions use a different parameterization from
#' \code{link[drda]{logistic2_gradient}}. To avoid the non-negative
#' constraints of parameters, the gradient and Hessian computed here are for
#' the function with `eta2 = log(eta)`.
#'
#' Note that argument `theta` is on the original scale and not on the log scale.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the two parameters in the form
#'   `c(eta, phi)`.
#' @param delta value of delta parameter (either 1 or -1).
#'
#' @return Gradient or Hessian of the alternative parameterization evaluated at
#'   the specified point.
#'
#' @export
logistic2_gradient_2 <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + b

  q <- y * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  G <- matrix(1, nrow = k, ncol = 2)

  G[, 1] <- eta * t
  G[, 2] <- u

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(G)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the gradient at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    G[is_nan] <- 0
  }

  sign(delta) * G
}

#' @export
#'
#' @rdname logistic2_gradient_2
logistic2_hessian_2 <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + b

  q <- y * b
  r <- -eta * b

  u <- (r / f) / f

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- -y * (1 + eta * (2 / f - 1 / b) * q) * u
  H[, 2, 1] <- (1 + eta * (2 / f - 1 / b) * q) * u

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- (2 / f - 1 / b) * r * u

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(H)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the Hessian at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    H[is_nan] <- 0
  }

  sign(delta) * H
}

#' @export
#'
#' @rdname logistic2_gradient_2
logistic2_gradient_hessian_2 <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + b

  q <- y * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  G <- matrix(1, nrow = k, ncol = 2)

  G[, 1] <- eta * t
  G[, 2] <- u

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- -y * (1 + eta * (2 / f - 1 / b) * q) * u
  H[, 2, 1] <- (1 + eta * (2 / f - 1 / b) * q) * u

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- (2 / f - 1 / b) * r * u

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(G)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the gradient at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    G[is_nan] <- 0
  }

  is_nan <- is.nan(H)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the Hessian at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    H[is_nan] <- 0
  }

  list(G = sign(delta) * G, H = sign(delta) * H)
}

# 2-parameter logistic function
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# 2-parameter logistic function.
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name) while vector `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization with `log(eta)`.
#
# @param object object of class `logistic2`.
# @param theta numeric vector with the parameters in the form `c(eta, phi)`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
#
#' @export
gradient_hessian.logistic2 <- function(object, theta) {
  logistic2_gradient_hessian_2(object$stats[, 1], theta, object$start[2])
}

# Residual sum of squares
#
# Evaluate the residual sum of squares (RSS) against the mean of a
# 2-parameter logistic model.
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name), while `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization with `log(eta)`.
#
# @param object object of class `logistic2`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(p)` to evaluate the RSS associated to a particular
#   parameter choice `p`.
#
#' @export
rss.logistic2 <- function(object) {
  function(theta) {
    theta[1] <- exp(theta[1])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.logistic2
#
#' @export
rss_fixed.logistic2 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 2)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[1] <- exp(theta[1])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# Residual sum of squares
#
# Evaluate the gradient and Hessian of the residual sum of squares (RSS)
# against the mean of a 2-parameter logistic model.
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name), while `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization with `log(eta)`.
#
# @param object object of class `logistic2`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#   the RSS associated to a particular parameter choice `theta`.
#
#' @export
rss_gradient_hessian.logistic2 <- function(object) {
  function(theta) {
    theta[1] <- exp(theta[1])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(object$m, 2, 2))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

# @rdname rss_gradient_hessian.logistic2
#
#' @export
rss_gradient_hessian_fixed.logistic2 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 2)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[1] <- exp(theta[1])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(object$m, 2, 2))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)

    list(
      G = apply(gradient[, idx, drop = FALSE], 2, sum),
      H = apply(hessian[, idx, idx, drop = FALSE], 2:3, sum)
    )
  }
}

# Maximum likelihood estimators
#
# Given a set of parameters, compute the maximum likelihood estimates of the
# lower and upper horizontal asymptotes.
#
# @param object object of class `logistic2`.
# @param theta vector of parameters.
#
# @return Numeric vector `theta`.
mle_asy.logistic2 <- function(object, theta) {
  names(theta) <- NULL
  theta
}

# Initialize vector of parameters
#
# Given the sufficient statistics, try to guess a good approximation to the
# Maximum Likelihood estimator of the four parameters of the logistic function.
#
# @param object object of class `logistic2`.
#
# @return Numeric vector of length 2 with a (hopefully) good starting point.
#
#' @importFrom stats lm
#'
#' @noRd
init.logistic2 <- function(object) {
  stats <- object$stats
  rss_fn <- rss(object)

  theta <- if (any(is.na(object$start))) {
    # data might not be compatible with a 2-parameter log-logistic function
    idx <- (stats[, 3] >= 0) & (stats[, 3] <= 1)

    if (mean(idx) <= 0.3) {
      warning(
        paste(
          "Too many values outside the range (0, 1).",
          "Is the data properly scaled?"
        ),
        call. = FALSE
      )
    }

    xx <- stats[idx, 1]
    yy <- stats[idx, 3]

    # we cannot guarantee that all values are between zero and one, so we assume
    # a 4-parameter logistic model as a starting point
    #
    # y = a + (b - a) / (1 + exp(-e * (x - p)))
    # w = (y - a) / (b - a) = 1 / (1 + exp(-e * (x - p)))
    #
    # by construction w is defined in (0, 1).
    #
    # z = log(w / (1 - w)) = - e * p + e * x
    #
    # fit a linear model `z ~ u0 + u1 x` and set `eta = u1` and `phi = -u0 / u1`
    #
    # we add a very small number to avoid the logarithm of zero.
    zv <- (yy + 1.0e-8) / (1 + 2.0e-8)
    zv <- log(zv) - log1p(-zv)
    tmp <- lm(zv ~ xx)

    # the curve can either increase of decrease depending on the `alpha` and
    # `delta` parameter. However, we want `eta` to be positive. If `eta` is
    # negative we simply change its sign and switch curve direction.
    log_eta <- log(abs(tmp$coefficients[2]))
    phi <- -tmp$coefficients[1] / tmp$coefficients[2]

    c(log_eta, phi)
  } else {
    object$start
  }

  best_rss <- rss_fn(theta)

  # this is a space-filling design using a max entropy grid
  v <- 250
  param_set <- matrix(
    c(
      # log_eta
      -9.82, -9.78, -9.74, -9.71, -9.64, -9.55, -9.55, -9.51, -9.48, -9.46,
      -9.25, -8.95, -8.83, -8.75, -8.73, -8.7, -8.6, -8.58, -8.57, -8.56, -8.55,
      -8.34, -8.28, -8.19, -8.18, -8.08, -8.06, -8.01, -7.95, -7.8, -7.74,
      -7.73, -7.67, -7.56, -7.5, -7.46, -7.4, -7.25, -7.18, -7.17, -7.17, -7.05,
      -7.05, -7.04, -6.95, -6.92, -6.92, -6.92, -6.89, -6.66, -6.63, -6.56,
      -6.54, -6.38, -6.32, -6.26, -6.26, -6.19, -6.14, -6.14, -6.08, -5.99,
      -5.9, -5.89, -5.88, -5.88, -5.77, -5.75, -5.72, -5.65, -5.54, -5.51,
      -5.36, -5.35, -5.32, -5.28, -5.1, -5.08, -5.03, -5, -4.98, -4.94, -4.93,
      -4.92, -4.9, -4.9, -4.82, -4.73, -4.67, -4.48, -4.46, -4.43, -4.37, -4.32,
      -4.29, -4.27, -4.23, -4.17, -4.17, -4.06, -4.02, -3.94, -3.94, -3.92,
      -3.92, -3.87, -3.77, -3.7, -3.65, -3.59, -3.57, -3.57, -3.55, -3.3, -3.3,
      -3.28, -3.27, -3.18, -3.15, -3.09, -3.04, -2.97, -2.93, -2.91, -2.88,
      -2.88, -2.85, -2.82, -2.77, -2.72, -2.64, -2.64, -2.61, -2.61, -2.54,
      -2.3, -2.25, -2.21, -2.18, -2.17, -2.15, -2.09, -2, -1.97, -1.93, -1.88,
      -1.84, -1.83, -1.82, -1.79, -1.6, -1.6, -1.57, -1.56, -1.55, -1.52, -1.36,
      -1.35, -1.23, -1.23, -1.14, -1.13, -1.11, -1.02, -1, -0.99, -0.95, -0.89,
      -0.78, -0.73, -0.62, -0.61, -0.58, -0.5, -0.46, -0.41, -0.39, -0.23, -0.2,
      -0.17, -0.15, -0.14, -0.1, -0.05, 0.03, 0.04, 0.1, 0.1, 0.16, 0.19, 0.31,
      0.36, 0.36, 0.37, 0.49, 0.49, 0.51, 0.61, 0.63, 0.75, 0.81, 0.81, 0.84,
      0.87, 0.93, 1.12, 1.24, 1.27, 1.32, 1.49, 1.52, 1.56, 1.58, 1.68, 1.72,
      1.75, 1.92, 2, 2.09, 2.1, 2.1, 2.18, 2.19, 2.23, 2.23, 2.46, 2.53, 2.55,
      2.67, 2.71, 2.72, 2.79, 2.86, 2.88, 2.94, 3.03, 3.04, 3.15, 3.21, 3.52,
      3.54, 3.63, 3.66, 3.72, 3.74, 3.74, 3.75, 3.88, 3.9, 3.91,
      # phi
      -10.72, 7.43, -4.6, -6.85, -0.37, -2.6, -8.5, 3.38, 9.46, -13.59, -16.39,
      -3.49, -14.78, 0.74, -10.22, -6.25, -1.7, -12.81, 4.79, -8.57, 3.19, 1.81,
      -7.21, -4.09, -11.67, 0.15, -16.69, -2.81, -9.15, 6.01, -1.19, -5.17,
      -14.25, -7.96, 2.49, -4.26, -6.55, -11.14, 1.26, 6.58, -12.95, -2.64,
      4.84, -9.6, 8.9, 3.17, -7.53, -18.06, 0.26, -4.35, -10.42, -6.47, -13.03,
      -1.33, -14.2, 2.12, 6.76, -8.75, 0.46, 4.01, -2.73, -11.54, -15.22, -4.41,
      -17.91, -1.17, 9, -10.2, -6.52, 1.4, -7.86, -12.49, -2.92, -13.63, 5.29,
      -0.02, -6.47, -15.27, -10.03, -8.35, 6.57, -1.93, 8.68, -4.72, 3.96, 2.14,
      -17.65, -13.44, -11.54, -0.08, 5.82, -9.7, 8.17, -15.14, -3.17, -17.6,
      1.89, -7.01, -3.98, -1.19, 3.47, 9.37, -5.03, -11.52, -14.03, 0.64, 6.51,
      -8.92, 4.92, 7.72, -2.57, -18.31, -0.95, -11.32, -16.32, 1.84, -9.45,
      -7.06, 2.95, 4.18, -4.35, 7.3, 9.69, -14.05, 8.25, -12.54, -8.7, -2.83,
      -5.73, -0.28, 1.32, 6.36, -1.26, -11.39, -4.62, -9.87, -14.09, -18.06,
      7.83, 4.47, -15.74, -8.37, -12.23, -2.55, 2.96, -5.32, 0.52, -4.18, -6.71,
      5.89, -13.33, -1.48, -9.77, -18.29, 7.12, 9.22, -16.75, 2.13, -5.8, -7.78,
      -14.42, 3.55, -2.58, -4.49, -11.28, -0.51, -9.87, 5.09, 7.42, -18.31,
      -1.37, -12.31, 3.21, 0.9, -7.48, -2.86, -16.28, -11.05, 6.88, -5.19,
      -14.1, -3.68, -9.66, 1.38, -15.52, 9.36, -1.68, 3.95, -6.08, -17.54, 5.98,
      2.22, -0.54, -13.56, -7.35, -11.3, 0.78, -9.39, -3.76, -1.75, 9.39, -5.18,
      8.01, 4.42, -15.16, -8.3, -6.64, -10.14, 2.86, -3.93, -12.67, 0.16,
      -11.42, -14.22, -2.44, 6.6, 3.69, -6.1, 9.84, -17.15, -0.97, 1.77, -10.32,
      -3.67, 8.7, 5.01, -4.73, -7.33, -8.86, -2.51, -13.49, 4.2, -1.12, -6.06,
      -12.29, -15.43, -10.48, -16.71, 0.97, 2.51, 7.63, -3.8, -5.32, 5.64,
      -7.34, -1.17, -10.01, -14.72, -13.37, -8.87
    ),
    ncol = v, byrow = TRUE
  )

  theta_tmp <- matrix(nrow = 2, ncol = v)
  rss_tmp <- rep(10000, v)

  for (i in seq_len(v)) {
    theta_tmp[, i] <- param_set[, i]
    rss_tmp[i] <- rss_fn(param_set[, i])
  }

  # update the total iteration count
  # 1: initial crude estimation
  # 2: flat line approximation
  # v: total amount of grid points tested
  niter <- 2 + v

  # many points might have the same approximate RSS therefore we must choose the
  # most likely ones
  ord <- order(
    round(rss_tmp, digits = 4),
    apply(theta_tmp, 2, function(x) sqrt(crossprod(x)))
  )

  # select the best solution and other not so good solutions as starting points
  theta_1 <- theta_tmp[, ord[1]]
  theta_2 <- theta_tmp[, ord[5]]
  theta_3 <- theta_tmp[, ord[8]]

  if (object$constrained) {
    # fix the candidates to be within the constraints
    theta <- pmax(
      pmin(theta, object$upper_bound, na.rm = TRUE),
      object$lower_bound, na.rm = TRUE
    )
    theta_1 <- pmax(
      pmin(theta_1, object$upper_bound, na.rm = TRUE),
      object$lower_bound, na.rm = TRUE
    )
    theta_2 <- pmax(
      pmin(theta_2, object$upper_bound, na.rm = TRUE),
      object$lower_bound, na.rm = TRUE
    )
    theta_3 <- pmax(
      pmin(theta_3, object$upper_bound, na.rm = TRUE),
      object$lower_bound, na.rm = TRUE
    )
  }

  start <- cbind(theta, theta_1, theta_2, theta_3)

  tmp <- fit_nlminb(object, start, object$max_iter - niter)

  if (!is.infinite(tmp$rss) && (tmp$rss < best_rss)) {
    theta <- tmp$theta
    best_rss <- tmp$rss
  }

  niter <- niter + tmp$niter

  names(theta) <- NULL
  names(niter) <- NULL

  list(theta = theta, niter = niter)
}

# 2-parameter logistic fit
#
# Fit a 2-parameter logistic function to observed data with a Maximum
# Likelihood approach.
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name), while `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# @param object object of class `logistic2`.
#
# @return A list with the following components:
#   \describe{
#     \item{converged}{boolean. `TRUE` if the optimization algorithm converged,
#       `FALSE` otherwise.}
#     \item{iterations}{total number of iterations performed by the
#       optimization algorithm}
#     \item{constrained}{boolean. `TRUE` if optimization was constrained,
#       `FALSE` otherwise.}
#     \item{estimated}{boolean vector indicating which parameters were
#       estimated from the data.}
#     \item{coefficients}{maximum likelihood estimates of the model
#       parameters.}
#     \item{rss}{minimum value found of the residual sum of squares.}
#     \item{df.residual}{residual degrees of freedom.}
#     \item{fitted.values}{fitted mean values.}
#     \item{residuals}{residuals, that is response minus fitted values.}
#     \item{weights}{vector of weights used for the fit.}
#   }
#
#' @export
fit.logistic2 <- function(object) {
  solution <- find_optimum(object)

  # bring the parameters back to their natural scale
  theta <- c(object$start[1:2], exp(solution$optimum[1]), solution$optimum[2])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(c(FALSE, TRUE), c(2, 2)),
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = length(object$y) - 2,
    fitted.values = logistic2_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("logistic2_fit", "logistic")

  result
}

# @rdname fit.logistic2
#
#' @export
fit_constrained.logistic2 <- function(object) {
  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 2, 3)

  for (i in seq_len(2)) {
    lb_is_inf <- is.infinite(object$lower_bound[i])
    ub_is_inf <- is.infinite(object$upper_bound[i])

    if (object$lower_bound[i] == object$upper_bound[i]) {
      constraint[i, 2] <- TRUE
    } else if (!lb_is_inf || !ub_is_inf) {
      constraint[i, 3] <- TRUE
    } else {
      constraint[i, 1] <- TRUE
    }
  }

  known_param <- ifelse(constraint[, 2], object$lower_bound, NA_real_)

  solution <- find_optimum_constrained(object, constraint, known_param)

  # bring the parameters back to their natural scale
  theta <- object$lower_bound
  theta[!constraint[, 2]] <- solution$optimum
  theta <- c(object$start[1:2], exp(theta[1]), theta[2])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = c(FALSE, FALSE, estimated),
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = length(object$y) - sum(estimated),
    fitted.values = logistic2_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("logistic2_fit", "logistic")

  result
}

# 2-parameter logistic fit
#
# Evaluate the Fisher information matrix at the maximum likelihood estimate.
#
# @details
# Let `mu(x; theta)` be the 2-parameter logistic function. We assume that our
# observations `y` are independent and such that
# `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
# distribution `N(0, 1)`.
#
# The 2-by-2 (symmetric) Fisher information matrix is the expected value of
# the negative Hessian matrix of the log-likelihood function.
#
# @param object object of class `logistic2`.
# @param theta numeric vector with the model parameters.
# @param sigma estimate of the standard deviation.
#
# @return Fisher information matrix evaluated at `theta`.
#
#' @export
fisher_info.logistic2 <- function(object, theta, sigma) {
  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]
  z <- fn(object, x, theta[3:4]) - y

  gh <- logistic2_gradient_hessian(x, theta[3:4], theta[2])

  # in case of theta being the maximum likelihood estimator, this gradient G
  # should be zero. We compute it anyway because we likely have rounding errors
  # in our estimate.
  G <- matrix(0, nrow = object$m, ncol = 2)
  G[, 1] <- w * z * gh$G[, 1]
  G[, 2] <- w * z * gh$G[, 2]

  G <- apply(G, 2, sum)

  H <- array(0, dim = c(object$m, 2, 2))

  H[, , 1] <- w * (z * gh$H[, , 1] + gh$G[, 1] * gh$G)
  H[, , 2] <- w * (z * gh$H[, , 2] + gh$G[, 2] * gh$G)

  H <- apply(H, 2:3, sum)

  mu <- fn(object, object$x, theta[3:4])
  z <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - sum(object$w > 0)

  fim <- rbind(cbind(H, -2 * G / sigma), c(-2 * G / sigma, z)) / sigma^2

  lab <- c(names(theta)[-seq_len(2)], "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

# 2-parameter logistic fit
#
# Find the dose that produced the observed response.
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name), while `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# This function evaluates the inverse function of `f(x; theta)`, that is
# if `y = fn(x; theta)` then `x = inverse_fn(y; theta)`.
#
#' @export
inverse_fn.logistic2_fit <- function(object, y) {
  inverse_fn.logistic4_fit(object, y)
}

# 2-parameter logistic fit
#
# Evaluate at a particular point the gradient of the inverse logistic function.
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name), while `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# This function evaluates the gradient of the inverse function.
#
#' @export
inverse_fn_gradient.logistic2_fit <- function(object, y) {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]

  z <- delta / (y - alpha)
  u <- 1 / (z - 1)

  G <- matrix(1, nrow = length(y), ncol = 2)
  G[, 1] <- -log(u) / eta^2

  G
}
