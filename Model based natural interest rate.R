# ============================================================
# MODEL-IMPLIED natural real rate (LW-style):
#   r*_t = gamma_r * g_t + z_t
#
# Uses states from your CB-LW KFAS model:
#   state 2 = g_t (trend growth, in 100*log units per quarter)
#   state 5 = rho_t (risk premium, pp)  [in your current code: rho_pp]
#
# IMPORTANT:
# - In your main model you CENTERED the spread: spr_pp_c = spr_pp - mean(spr_pp)
#   => rho_t is a *deviation* around zero.
#   For z_t in r* = gamma*g + z, you usually want the *level* premium:
#        z_t = rho_dev_t + mean(spread in pp)
#
# Output:
#   - rstar_model (point)
#   - optional 95% CI from smoothed state variance
# ============================================================

suppressPackageStartupMessages({ library(zoo) })

stopifnot(exists("kfs"))
stopifnot(exists("D"))
stopifnot("spread_raw" %in% names(D))   # from your Croatia CDS merge (bps typically)

AH <- kfs$alphahat
n  <- nrow(AH)

# --- aligned dates ---
dts <- if (exists("dates") && length(dates) >= n) dates[seq_len(n)] else D$date[seq_len(n)]
dts <- as.yearqtr(dts)

# --- USER CHOICE: gamma in r* = gamma*g + z ---
# NOTE: this is NOT your inflation-import gamma. Use a different name.
gamma_r <- 1.0   # <-- set this to your chosen/calibrated value (e.g., 1, 2, 0.5)

# --- 1) Extract g_t and annualize to percent ---
# In your model: g_trend = state 2 (quarterly change in y* where y=100*log)
# Annualized percent: 4 * g_trend
g_q        <- as.numeric(AH[, 2])
g_ann_pct  <- 4 * g_q

# --- 2) Extract rho_t (deviation, pp) ---
rho_dev_pp <- as.numeric(AH[, 5])

# --- 3) Build z_t (risk premium level, pp) ---
# spread_raw is likely in bps -> convert to pp
spr_pp <- as.numeric(D$spread_raw)[seq_len(n)] / 100

# mean spread in pp over the estimation sample (same centering you used)
spr_mean_pp <- mean(spr_pp, na.rm = TRUE)

# level premium:
z_pp <- rho_dev_pp + spr_mean_pp

# --- 4) Model-implied natural rate ---
rstar_model <- gamma_r * g_ann_pct + z_pp   # % points

# ============================================================
# 95% CI (from smoothed state covariance V_t)
# ============================================================
add_ci <- TRUE

if (add_ci) {
  # KFAS: kfs$V[i,j,t] is smoothed state covariance
  g_var   <- as.numeric(kfs$V[2, 2, 1:n])          # var(g_q)
  rho_var <- as.numeric(kfs$V[5, 5, 1:n])          # var(rho_dev)
  
  # cov(g_q, rho_dev) if you want exact variance:
  g_rho_cov <- as.numeric(kfs$V[2, 5, 1:n])
  
  # Var(4*g_q) = 16*Var(g_q)
  var_g_ann <- 16 * pmax(g_var, 0)
  
  # z_pp = rho_dev + constant => Var(z_pp) = Var(rho_dev)
  var_z <- pmax(rho_var, 0)
  
  # Var(gamma*g_ann + z) = gamma^2 Var(g_ann) + Var(z) + 2*gamma*Cov(g_ann,z)
  # Cov(g_ann, z) = Cov(4*g_q, rho_dev) = 4*Cov(g_q, rho_dev)
  var_rstar <- (gamma_r^2) * var_g_ann + var_z + 2 * gamma_r * (4 * g_rho_cov)
  
  se_rstar <- sqrt(pmax(var_rstar, 0))
  zcrit <- 1.96
  lo <- rstar_model - zcrit * se_rstar
  hi <- rstar_model + zcrit * se_rstar
}

# --- Plot ---
yl <- if (add_ci) range(c(lo, hi, rstar_model), na.rm=TRUE) else range(rstar_model, na.rm=TRUE)

plot(dts, rstar_model, type="l", lwd=2, ylim=yl,
     main = bquote("Model-implied natural real rate"~r[t]^"*"
                   ~": "~r[t]^"*"==gamma[r]*g[t]+z[t]),
     xlab="", ylab="% points")

if (add_ci) {
  ord <- order(dts)
  polygon(
    x = c(as.numeric(dts[ord]), rev(as.numeric(dts[ord]))),
    y = c(lo[ord], rev(hi[ord])),
    col = "grey85", border = NA
  )
  lines(dts, rstar_model, lwd=2)
  lines(dts, lo, lty=3)
  lines(dts, hi, lty=3)
}

abline(h=0, lty=3, col="grey60")

legend("topleft",
       legend = if (add_ci)
         c(expression(r[t]^"* (point)"), "95% CI")
       else
         c(expression(r[t]^"* (point)")),
       lty = if (add_ci) c(1, NA) else c(1),
       lwd = if (add_ci) c(2, 8)  else c(2),
       pch = if (add_ci) c(NA, 15) else c(NA),
       col = if (add_ci) c("black","grey85") else c("black"),
       bty="n")

# --- quick sanity prints ---
cat("\n--- SANITY CHECKS ---\n")
cat("gamma_r =", gamma_r, "\n")
cat("g_ann_pct:  min/mean/max =", min(g_ann_pct,na.rm=TRUE),
    mean(g_ann_pct,na.rm=TRUE), max(g_ann_pct,na.rm=TRUE), "\n")
cat("rho_dev_pp: min/mean/max =", min(rho_dev_pp,na.rm=TRUE),
    mean(rho_dev_pp,na.rm=TRUE), max(rho_dev_pp,na.rm=TRUE), "\n")
cat("spr_mean_pp (added back) =", spr_mean_pp, "\n")
cat("rstar_model: min/mean/max =", min(rstar_model,na.rm=TRUE),
    mean(rstar_model,na.rm=TRUE), max(rstar_model,na.rm=TRUE), "\n")
