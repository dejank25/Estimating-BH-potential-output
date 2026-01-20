# ============================================================
# MODEL-BASED natural real rate for BiH (from your KFAS CB-LW):
#   r*_BiH,model(t) = r*_EA(t) + rho_t
#
# Uses:
#   - D$r_star   : EA neutral real rate proxy (from your chunk 6)
#   - kfs        : KFS output of your main CB-LW model
#   - rho state  : AH[,5] in your current specification
#
# Optionally adds 95% CI with fixed spread measurement sd (e.g., 25 bps)
# ============================================================

suppressPackageStartupMessages({
  library(zoo)
})

# --- checks ---
stopifnot(exists("kfs"))
stopifnot(exists("D"))
stopifnot("r_star" %in% names(D))

AH <- kfs$alphahat
n  <- nrow(AH)

# --- aligned dates ---
dts <- if (exists("dates") && length(dates) >= n) dates[seq_len(n)] else D$date[seq_len(n)]
dts <- as.yearqtr(dts)

# --- components ---
rho_hat   <- as.numeric(AH[, 5])                 # rho_t (pp) from main model (state 5)
rstar_hat <- as.numeric(D$r_star)[seq_len(n)]    # r*_EA(t) (pp)

# --- model-based BiH natural rate ---
rstar_bih_model <- rstar_hat + rho_hat

# ============================================================
# OPTIONAL: 95% CI using a fixed spread measurement error sd
# (recommended when KFAS estimates spread ME ~ 0 and CIs disappear)
# ============================================================

add_ci <- TRUE
spread_me_sd_bps <- 25   # <-- change (e.g., 10, 15, 25); 25 bps = 0.25 pp

if (add_ci) {
  z <- 1.96
  ci_half_pp <- z * (spread_me_sd_bps / 100)   # convert bps -> pp, then 95% half-width
  ci_lo <- rstar_bih_model - ci_half_pp
  ci_hi <- rstar_bih_model + ci_half_pp
}

# --- plot ---
yl <- if (add_ci) range(c(ci_lo, ci_hi, rstar_hat), na.rm = TRUE) else range(c(rstar_bih_model, rstar_hat), na.rm = TRUE)

plot(dts, rstar_bih_model, type="l", lwd=2, ylim=yl,
     main = if (add_ci)
       paste0("Model-based natural real rate r* for BiH (CB): r*_EA + rho (95% CI)\n",
              "Spread ME sd fixed at ", spread_me_sd_bps, " bps")
     else
       "Model-based natural real rate r* for BiH (CB): r*_EA + rho",
     xlab="", ylab="% points")

# EA anchor
lines(dts, rstar_hat, lty=2, lwd=2)

# CI band (clean shading + borders)
if (add_ci) {
  ord <- order(dts)
  polygon(
    x = c(as.numeric(dts[ord]), rev(as.numeric(dts[ord]))),
    y = c(ci_lo[ord], rev(ci_hi[ord])),
    col = "grey85", border = NA
  )
  # redraw main lines on top of shading
  lines(dts, rstar_bih_model, lwd=2)
  lines(dts, rstar_hat, lty=2, lwd=2)
  lines(dts, ci_lo, lty=3)
  lines(dts, ci_hi, lty=3)
}

abline(h=0, lty=3, col="grey60")

legend("bottomleft",
       legend = if (add_ci)
         c("r*_BiH, model (point)",  "95% CI")
       else
         c("r*_BiH, model (point)", "r*_EA (point)"),
       lty = if (add_ci) c(1,2,NA) else c(1,2),
       lwd = if (add_ci) c(2,2,8)  else c(2,2),
       pch = if (add_ci) c(NA,NA,15) else c(NA,NA),
       pt.cex = if (add_ci) c(NA,NA,1.5) else c(NA,NA),
       col = if (add_ci) c("black","black","grey85") else c("black","black"),
       bty="n")
