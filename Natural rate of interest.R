# ============================================================
# Natural real rate for BiH under Currency Board (CB-LW):
#   r*_BiH(t) = r*_EA(t) + rho_dev(t)
#
# Requires from your CB-LW script:
#   - D with columns: date, r_ea (or euribor6m_ea & pi_ea_exp), and r_star (EA neutral real rate)
#     (In your script: r_star was computed from the local-level model on r_ea)
#   - kfs from the CB-LW fit, with state 4 = rho (premium deviation, pp)
#   - dates (yearqtr) aligned to kfs$alphahat rows (your CB-LW dates)
# ============================================================

suppressPackageStartupMessages({
  library(zoo)
})

# --- 0) Grab CB-LW smoothed states safely ---
stopifnot(exists("kfs"))
AH <- kfs$alphahat
n  <- nrow(AH)

# --- 1) Dates aligned to CB-LW output ---
# Prefer the CB-LW dates object you used for plots:
if (exists("dates") && length(dates) >= n) {
  dts <- dates[seq_len(n)]
} else if (exists("D") && "date" %in% names(D) && length(D$date) >= n) {
  dts <- D$date[seq_len(n)]
} else {
  stop("Can't find aligned dates. Use the same 'dates' you used in the CB-LW plots.")
}

# --- 2) Extract rho deviation (state 4 in CB-LW) ---
rho_dev <- as.numeric(AH[, 5])  # "Premium rho (deviation, pp)"

# --- 3) Bring in EA neutral real rate r*_EA aligned to same sample ---
# Your CB-LW script computed:
#   r_ea = euribor6m_ea - lag(pi_ea,1)
#   r_star = local-level smoother on r_ea
# so we need r_star aligned with the same rows used in the CB-LW model.
#
# Best case: you kept r_star in D as a column (recommended).
if (exists("D") && "r_star" %in% names(D)) {
  rstar_ea <- as.numeric(D$r_star)[seq_len(n)]
} else if (exists("r_star")) {
  # If r_star exists as a standalone vector in memory
  rstar_ea <- as.numeric(r_star)[seq_len(n)]
} else {
  stop("EA r_star not found. Ensure your stance step saved r_star into D$r_star or keep r_star vector.")
}

# --- 4) Construct BiH natural real rate ---
rstar_bih <- rstar_ea + rho_dev  # both in percentage points

# --- 5) Plot (clean y-limits; no weird scientific ticks) ---
yl <- range(rstar_bih, rstar_ea, na.rm = TRUE)

plot(dts, rstar_bih, type="l", lwd=2,
     main="Natural real rate r* for BiH (CB): r*_EA + rho",
     xlab="", ylab="% points", ylim=yl)

lines(dts, rstar_ea, lty=2, lwd=2)

abline(h=0, lty=3, col="grey60")

legend("bottomleft",
       legend = c("r*_BiH = r*_EA + rho"),
#                  "r*_EA (neutral EA real rate)"),
       lty = 1, #c(1, 2),
       lwd = 2,
       bty = "n")
