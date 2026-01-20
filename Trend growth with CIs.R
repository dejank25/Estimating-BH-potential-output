# ============================================================
# Confidence intervals for trend growth g (annualized, %)
# Robust method: PARAMETRIC BOOTSTRAP over variance parameters
# Works across KFAS versions (no KFS tol/method/nsim quirks)
#
# Requires in workspace:
#   fit      : output from fitSSM(...)
#   updatefn : your update function (par -> model with updated Q/H)
#   dates    : yearqtr vector
# ============================================================

stopifnot(exists("fit"), exists("updatefn"), exists("dates"))

# Point estimate
kfs0 <- KFS(fit$model, smoothing = "state")
g_hat_ann <- 4 * as.numeric(kfs0$alphahat[,2])
Tn <- length(g_hat_ann)
stopifnot(length(dates) == Tn)

# --- Bootstrap settings ---
set.seed(123)
B <- 300   # 300 is usually enough; raise to 500 if you want smoother bands

par_hat <- fit$optim.out$par

# Cov matrix of estimates (BFGS approx). If this is NULL, we fall back to diagonal.
Vhat <- fit$optim.out$hessian
if (!is.null(Vhat)) {
  # invert Hessian to get covariance (approx)
  Vcov <- tryCatch(solve(Vhat), error = function(e) NULL)
} else {
  Vcov <- NULL
}
if (is.null(Vcov) || any(!is.finite(Vcov))) {
  message("Warning: covariance not available; using diagonal fallback.")
  Vcov <- diag( (0.1 * abs(par_hat) + 0.05)^2 )
}

# Cholesky for sampling
L <- chol((Vcov + t(Vcov)) / 2 + diag(1e-10, length(par_hat)))

# Storage
g_draws <- matrix(NA_real_, nrow = Tn, ncol = B)

for (b in 1:B) {
  # draw parameters in unconstrained (log-variance) space
  z <- rnorm(length(par_hat))
  par_b <- as.numeric(par_hat + L %*% z)
  
  # build model with these parameters (NO re-optimization)
  mod_b <- updatefn(par_b, fit$model)
  
  # smooth
  kfs_b <- KFS(mod_b, smoothing = "state")
  
  # store annualized g
  g_draws[, b] <- 4 * as.numeric(kfs_b$alphahat[,2])
}

# Quantile CI
g_lo <- apply(g_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
g_hi <- apply(g_draws, 1, quantile, probs = 0.975, na.rm = TRUE)

# Plot (robust numeric x-axis)
# --- Robust x and y limits ---
x_num <- as.numeric(dates)
xlim_num <- range(x_num) + c(-0.5, 0.5)

# y-limits must include BOTH CI bands
ylim_num <- range(c(g_hat_ann, g_lo, g_hi), na.rm = TRUE)
# add a little padding
pad <- 0.05 * diff(ylim_num)
ylim_num <- ylim_num + c(-pad, pad)

plot(x_num, g_hat_ann,
     type = "l",
     lwd = 2,
     xaxt = "n",
     xlim = xlim_num,
     ylim = ylim_num,
     main = "Trend growth g (annualized, %) with 95% bootstrap CI",
     xlab = "",
     ylab = "%")

lines(x_num, g_lo, lty = 2)
lines(x_num, g_hi, lty = 2)

yrs <- seq(floor(min(x_num)), ceiling(max(x_num)))
axis(1, at = yrs, labels = yrs)
