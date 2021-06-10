## ----echo = FALSE, include = FALSE----------------------------------
options(prompt = "R> ", continue = "+ ", width = 70, useFancyQuotes = FALSE)

knitr::render_sweave()
knitr::opts_chunk$set(
  fig.pos = "ht!", comment = "", dev.args = list(pointsize = 14), cache = TRUE
)

colpal <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db",
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
  "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d"
)

set.seed(4336273)

## ----logistic5, echo = FALSE, fig.height = 10, fig.width = 10, fig.cap = "Generalized (5-parameter) logistic function with various choices of parameters."----
fnl5 <- function(x, y) {
  a <- y[1]
  b <- y[2]
  e <- y[3]
  p <- y[4]
  n <- y[5]

  a + (b - a) * (1 + n * exp(-e * (x - p)))^(-1 / n)
}

old_par <- par(mfrow = c(2, 2))

curve(
  fnl5(x, c(0.1, 0.9, -1, 0, 1)), from = -10, to = 10, col = colpal[1],
  ylab = expression(paste(mu, "(x;", psi, ")")), xlab = "x", ylim = c(0, 1),
  main = expression(paste("Varying ", phi, " parameter", sep = ""))
)
mtext(expression(paste(alpha, "= 0.1,", beta, "= 0.9,", eta, "= -1,", nu, "= 1")))
curve(fnl5(x, c(0.1, 0.9, -1, 5, 1)), add = TRUE, col = colpal[2])
curve(fnl5(x, c(0.1, 0.9, -1, 2, 1)), add = TRUE, col = colpal[3])
curve(fnl5(x, c(0.1, 0.9, -1, -2, 1)), add = TRUE, col = colpal[4])
curve(fnl5(x, c(0.1, 0.9, -1, -5, 1)), add = TRUE, col = colpal[5])
abline(h = c(0.1, 0.9), lty = 2)
legend(
  "topright", title = expression(phi), legend = c(5, 2, 0, -2, -5),
  col = colpal[c(2:3, 1, 4:5)], lty = 1, bg = "white"
)

curve(
  fnl5(x, c(0.05, 0.95, -1, 0, 1)), from = -10, to = 10, col = colpal[1],
  ylab = expression(paste(mu, "(x;", psi, ")")), xlab = "x", ylim = c(0, 1),
  main = expression(paste("Varying ", nu, " parameter", sep = ""))
)
mtext(expression(paste(alpha, "= 0.05,", beta, "= 0.95,", eta, "= -1,", phi, "= 0")))
curve(fnl5(x, c(0.05, 0.95, -1, 0, 0.1)), add = TRUE, col = colpal[2])
curve(fnl5(x, c(0.05, 0.95, -1, 0, 0.5)), add = TRUE, col = colpal[3])
curve(fnl5(x, c(0.05, 0.95, -1, 0, 2.5)), add = TRUE, col = colpal[4])
curve(fnl5(x, c(0.05, 0.95, -1, 0, 5)), add = TRUE, col = colpal[5])
abline(h = c(0.05, 0.95), lty = 2)
legend(
  "topright", title = expression(nu), legend = c(0.1, 0.5, 1, 2.5, 5),
  col = colpal[c(2:3, 1, 4:5)], lty = 1, bg = "white"
)

curve(
  fnl5(x, c(0, 1, -1, 0, 1)), from = -10, to = 10, col = colpal[1],
  ylab = expression(paste(mu, "(x;", psi, ")")), xlab = "x", ylim = c(0, 1),
  main = expression(paste("Varying ", eta, " parameter", sep = ""))
)
mtext(expression(paste(alpha, "= 0,", beta, "= 1,", phi, "= 0,", nu, "= 1")))
curve(fnl5(x, c(0, 1, -0.25, 0, 1)), add = TRUE, col = colpal[2])
curve(fnl5(x, c(0, 1, -0.5, 0, 1)), add = TRUE, col = colpal[3])
curve(fnl5(x, c(0, 1, -2, 0, 1)), add = TRUE, col = colpal[4])
curve(fnl5(x, c(0, 1, -5, 0, 1)), add = TRUE, col = colpal[5])
abline(h = c(0, 1), lty = 2)
legend(
  "topright", title = expression(eta), legend = c(-0.25, -0.5, -1, -2, -5),
  col = colpal[c(2:3, 1, 4:5)], lty = 1, bg = "white"
)

curve(
  fnl5(x, c(0, 1, 1, 0, 1)), from = -10, to = 10, col = colpal[1],
  ylab = expression(paste(mu, "(x;", psi, ")")), xlab = "x", ylim = c(0, 1),
  main = expression(paste("Varying ", eta, " parameter", sep = ""))
)
mtext(expression(paste(alpha, "= 0,", beta, "= 1,", phi, "= 0,", nu, "= 1")))
curve(fnl5(x, c(0, 1, 0.25, 0, 1)), add = TRUE, col = colpal[2])
curve(fnl5(x, c(0, 1, 0.5, 0, 1)), add = TRUE, col = colpal[3])
curve(fnl5(x, c(0, 1, 2, 0, 1)), add = TRUE, col = colpal[4])
curve(fnl5(x, c(0, 1, 5, 0, 1)), add = TRUE, col = colpal[5])
abline(h = c(0, 1), lty = 2)
legend(
  "bottomright", title = expression(eta), legend = c(0.25, 0.5, 1, 2, 5),
  col = colpal[c(2:3, 1, 4:5)], lty = 1, bg = "white"
)

par(old_par)

## ----objfnproblem, echo = FALSE, fig.height = 6, fig.width = 10, fig.cap = "Problematic real data (cell line: BT-20, compound: BI-2536, dataset: CTRPv2) \\citep{rees_2016_ncb_correlating, seashore-ludlow_2015_cd_harnessing, basu_2013_cell_interactive}). A) 4-parameter logistic function as fitted by the BFGS algorithm. Starting point \\(\\boldsymbol{\\psi} = (\\alpha, \\beta, \\eta, \\phi)^{T} = (0, 1, -1, 0)^{T}\\). B) Contour plot of the residual sum of squares \\(g(\\boldsymbol{\\psi})\\) with respect to parameters \\(\\eta\\) and \\(\\phi\\). Fixed parameters \\(\\alpha = 0\\) and \\(\\beta = 1\\)."----
fig_x <- c(
  -6.90775527898214, -6.21460809842219, -5.49676830527187, -4.81589121730374,
  -4.13516655674236, -3.44201937618241, -2.7333680090865, -2.04022082852655,
  -1.34707364796661, -0.653926467406664, 0, 0.741937344729377,
  1.43508452528932, 2.11625551480255, 2.83321334405622, 3.49650756146648
)

fig_y <- c(
  0.9953, 1.074, 0.6401, 0.5836,
  0.5796, 0.6442, 0.5219, 0.625,
  0.5991, 0.652, 0.6246, 0.6743,
  0.577, 0.6559, 0.5197, 0.1061
)

fig_theta_drda <- c(
  0.566384615286767, 1.03465000008177, -30.6814903350415, -5.55144408247884
)

fig_theta_optim <- c(
  -3.20708018302608, 4.37542010565651, -0.0217950331820467, -0.598560483094873
)

fig_fn <- function(x, y) {
  y[1] + (y[2] - y[1]) / (1 + exp(-y[3] * (x - y[4])))
}

fig_rss <- function(x) {
  mu <- fig_fn(fig_x, c(0, 1, x[1], x[2]))
  sum((fig_y - mu)^2) / 2
}

N <- 400
eta_set <- seq(-2, 0, length.out = N)
phi_set <- seq(-20, 20, length.out = N)
rss_val <- matrix(
  apply(expand.grid(eta_set, phi_set), 1, fig_rss), nrow = N, ncol = N
)

old_par <- par(mfrow = c(1, 2))

plot(
  fig_x, fig_y, type = "p", xlab = "log(dose)", ylab = "Percent viability",
  ylim = c(0, 1.2)
)
curve(
  fig_fn(x, fig_theta_drda),
  add = TRUE, lty = 2, lwd = 2, col = "#EE6677FF"
)
curve(
  fig_fn(x, fig_theta_optim),
  add = TRUE, lty = 2, lwd = 2, col = "#4477AAFF"
)
legend(
  "bottomleft", legend = c("True estimate", "BFGS"), lty = 2, lwd = 2,
  bg = "white", col = c("#EE6677FF", "#4477AAFF"), bty = "n"
)
title("A)", adj = 0)

contour(
  x = eta_set, y = phi_set, z = rss_val,
  levels = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 3.0, 3.4),
  xlab = expression(eta), ylab = expression(phi),
)
title("B)", adj = 0)

par(old_par)

## -------------------------------------------------------------------
library(drda)

## -------------------------------------------------------------------
dose <- rep(c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), each = 3)
relative_viability <- c(
  0.877362, 0.812841, 0.883113, 0.873494, 0.845769, 0.999422, 0.888961,
  0.735539, 0.842040, 0.518041, 0.519261, 0.501252, 0.253209, 0.083937,
  0.000719, 0.049249, 0.070804, 0.091425, 0.041096, 0.000012, 0.092564
)

## -------------------------------------------------------------------
fit <- drda(relative_viability ~ dose, is_log = FALSE)

## -------------------------------------------------------------------
log_dose <- log(dose)
test_data <- data.frame(d = dose, x = log_dose, y = relative_viability)
# the following calls are equivalent
fit <- drda(relative_viability ~ log_dose)
fit <- drda(y ~ d, data = test_data, is_log = FALSE)
fit <- drda(y ~ x, data = test_data)

## -------------------------------------------------------------------
summary(fit)

## -------------------------------------------------------------------
coef(fit)
fit$coefficients

## -------------------------------------------------------------------
sigma(fit)
fit$sigma

## -------------------------------------------------------------------
deviance(fit)

## -------------------------------------------------------------------
residuals(fit)

## -------------------------------------------------------------------
logLik(fit)

## -------------------------------------------------------------------
predict(fit)

## -------------------------------------------------------------------
predict(fit, x = log(c(0.002, 0.2, 2)))

## -------------------------------------------------------------------
fit_logi2 <- drda(y ~ x, data = test_data, mean_function = "logistic2")
anova(fit_logi2)

## -------------------------------------------------------------------
fit_logi4 <- drda(y ~ x, data = test_data, mean_function = "logistic4")
fit_gompe <- drda(y ~ x, data = test_data, mean_function = "gompertz")
anova(fit_logi2, fit_logi4, fit_gompe)

## -------------------------------------------------------------------
weights <- c(
  0.990868, 1.095238, 0.974544, 0.973318, 1.107001, 1.012844, 1.052806,
  1.019427, 1.032544, 0.919827, 0.971385, 0.959019, 1.037789, 1.006835,
  0.969383, 0.935633, 1.016597, 1.011085, 0.982307, 1.066032, 0.959870
)
fit_weights <- drda(y ~ x, data = test_data, weights = weights)
summary(fit_weights)

## -------------------------------------------------------------------
weights(fit_weights)

## -------------------------------------------------------------------
residuals(fit_weights, type = "weighted")

## -------------------------------------------------------------------
lb <- c(0, 1, -5, -Inf)
ub <- c(0, 1,  5,  Inf)
fit_cnstr <- drda(
  y ~ x, data = test_data, lower_bound = lb, upper_bound = ub
)
summary(fit_cnstr)

## -------------------------------------------------------------------
fit_cnstr <- drda(
  y ~ x, data = test_data, lower_bound = lb, upper_bound = ub,
  start = c(0, 1, -0.6, -2), max_iter = 10000
)
summary(fit_cnstr)

## ----plot_logi5, fig.pos = "H", fig.height = 6, fig.width = 10, fig.cap = ""----
fit_logi5 <- drda(y ~ x, data = test_data, mean_function = "logistic5")
plot(fit_logi5)

## ----plot_multi, fig.pos = "H", fig.height = 6, fig.width = 10, fig.cap = ""----
plot(
  fit_logi2, fit_logi4, fit_gompe,
  base = "10", level = 0.9,
  xlim = c(-10, 5), ylim = c(-0.1, 1.1),
  xlab = "Dose", ylab = "Relative viability",
  legend = c("2-param logistic", "4-param logistic", "Gompertz")
)

## -------------------------------------------------------------------
naac(fit_logi4)

## -------------------------------------------------------------------
naac(fit_logi4, xlim = c(-2, 2), ylim = c(0.1, 0.9))

