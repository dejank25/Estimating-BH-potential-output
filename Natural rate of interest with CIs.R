# ============================================================
# NATURAL RATE (BiH) with VISIBLE / MEANINGFUL 95% CI
# Fix: set a realistic measurement error for the spread equation (h_s),
# so rho uncertainty is not ~0 and CI is visible.
#
# You only need to:
#   1) Use the UPDATED "par0 + updatefn" block below inside your main CB-LW script
#      (it fixes h_s instead of estimating it).
#   2) Then run the CI plotting block at the end.
#
# Recommended: spread_me_sd_bps = 15 (i.e., 15 bps measurement noise)
# Try 10–25 bps as robustness.
# ============================================================

# ---------------------------
# A) IN YOUR MAIN MODEL: replace ONLY par0 + updatefn + fitSSM call
#    (keep everything else: Zmat/Tmat/Rmat/model0 exactly as you have it)
# ---------------------------

# --- FIXED spread measurement sd (in basis points) ---
spread_me_sd_bps <- 25                 # <-- change 10/15/20/25 for robustness
spread_me_sd_pp  <- spread_me_sd_bps/100   # bps -> percentage points
h_s_fixed_var     <- spread_me_sd_pp^2

# Initial variances (log-scale) for estimation:
# We will STILL estimate all Q's and (h_y, h_pi),
# but we will FIX h_s at h_s_fixed_var (do not estimate it).
par0 <- log(c(
  q_y    = 1e-2,
  q_g    = 5e-3,
  q_xbar = 1e-4,
  q_xdev = 1e-1,
  q_rho  = 5e-3,
  q_pi   = 5e-3,
  h_y    = 1e-2,
  h_pi   = 1e-1
  # h_s removed from estimation on purpose
))

# Build initial Q/H (must match dimensions)
Q0 <- diag(exp(par0[1:6]), 6, 6); dimnames(Q0) <- NULL
H0 <- diag(c(exp(par0[7]), exp(par0[8]), h_s_fixed_var), 3, 3); dimnames(H0) <- NULL
storage.mode(Q0) <- "double"; storage.mode(H0) <- "double"

# IMPORTANT: rebuild model0 using H0/Q0 if your model0 currently used old H0/Q0
# If your model0 is already built, just overwrite:
model0$Q[,,1] <- Q0
model0$H[,,1] <- H0

# Update function: estimates Q + (h_y, h_pi), keeps h_s fixed
updatefn <- function(par, model) {
  qv <- exp(par[1:6])
  hy <- exp(par[7])
  hp <- exp(par[8])
  
  model$Q[,,1] <- diag(qv, 6, 6)
  model$H[,,1] <- diag(c(hy, hp, h_s_fixed_var), 3, 3)  # <-- FIXED h_s
  model
}

fit <- fitSSM(
  inits    = par0,
  model    = model0,
  updatefn = updatefn,
  method   = "BFGS",
  control  = list(maxit = 800, reltol = 1e-9)
)

kfs <- KFS(fit$model, smoothing = "state")

cat("\nCB-LW estimated with FIXED spread measurement noise:\n")
cat("spread_me_sd_bps =", spread_me_sd_bps, "bps\n")
cat("Converged:", fit$optim.out$convergence == 0, "\n")


# ---------------------------
# B) NATURAL RATE + 95% CI PLOT (run after the model above)
#    Uses:
#      - D$r_star (EA neutral proxy)
#      - kfs$alphahat and kfs$V for rho uncertainty
# ---------------------------

suppressPackageStartupMessages({ library(zoo) })

stopifnot(exists("kfs"))
stopifnot(exists("D"))
stopifnot("r_star" %in% names(D))

AH <- kfs$alphahat
n  <- nrow(AH)

# Aligned dates
dts <- if (exists("dates") && length(dates) >= n) dates[seq_len(n)] else D$date[seq_len(n)]
dts <- as.yearqtr(dts)

# Components
rho_hat   <- as.numeric(AH[, 5])            # rho (pp) = state 5 in your main model
rstar_hat <- as.numeric(D$r_star)[seq_len(n)]
rbih_hat  <- rstar_hat + rho_hat

# ---- CI for BiH natural rate ----
# Use rho uncertainty from main model (visible now because h_s fixed realistically)
rho_var <- as.numeric(kfs$V[5, 5, 1:n])
rho_se  <- sqrt(pmax(rho_var, 0))

# Option 1 (clean): CI driven by rho only (EA anchor treated as observed proxy)
z <- 1.96
rbih_lo <- rbih_hat - z * rho_se
rbih_hi <- rbih_hat + z * rho_se

# Optional: print CI width sanity
cat("\nCI half-width (pp): min/median/max = ",
    min(z*rho_se, na.rm=TRUE),
    median(z*rho_se, na.rm=TRUE),
    max(z*rho_se, na.rm=TRUE), "\n")

# ---- Plot: point lines + shaded CI band ----
yl <- range(c(rbih_lo, rbih_hi, rstar_hat), na.rm=TRUE)

plot(dts, rbih_hat, type="n", ylim=yl,
     main=paste0("Natural real rate r* for BiH (CB): r*_EA + rho (95% CI)\n",
                 "Spread ME sd fixed at ", spread_me_sd_bps, " bps"),
     xlab="", ylab="% points")

# CI shading (grey polygon)
polygon(c(dts, rev(dts)),
        c(rbih_lo, rev(rbih_hi)),
        col=grDevices::adjustcolor("grey70", alpha.f=0.6),
        border=NA)

# Lines
lines(dts, rbih_hat, lwd=2)          # BiH point
lines(dts, rstar_hat, lty=2, lwd=2)  # EA point (anchor)
abline(h=0, lty=3, col="grey60")

legend("bottomright",
       legend=c("r*_BiH (point)", "r*_BiH 95% CI"),
       lty=c(1, 2, NA),
       lwd=c(2, 2, NA),
       pch=c(NA, NA, 15),
       pt.cex=1.2,
       col=c("black","black","grey70"),
       bty="n")
