# ============================================================
# HP potential output + HP vs LW output gap (single copy/paste)
# REQUIRE (already in your LW workspace):
#   - dates : yearqtr vector
#   - y     : 100*log(real GDP) numeric (same you used in LW)
#   - LW gap either:
#       (a) x_dev  : your LW "cyclical" component in 100*log units  (you plot x_dev/100), OR
#       (b) y_star : LW potential in 100*log units (so LW gap = y - y_star)
#
# OUTPUT:
#   - y_hp_star : HP potential (100*log units)
#   - gap_hp_pct: HP output gap (%)
#   - gap_lw_pct: LW output gap (%)
#   - Plot: HP vs LW output gap overlay
# ============================================================

suppressPackageStartupMessages({
  library(zoo)
})

# --- 0) Safety + coercion (avoid ts/mts issues) ---
stopifnot(exists("dates"), exists("y"))
dates <- zoo::as.yearqtr(dates)
y     <- as.numeric(y)

# Align to same length (defensive)
n <- min(length(dates), length(y))
dates <- dates[seq_len(n)]
y     <- y[seq_len(n)]

stopifnot(n >= 20)

# --- 1) Compute HP potential output (quarterly lambda = 1600) ---
if (!requireNamespace("mFilter", quietly = TRUE)) {
  install.packages("mFilter")
}
hp <- mFilter::hpfilter(y, freq = 1600, type = "lambda")  # y can be plain numeric
y_hp_star <- as.numeric(hp$trend)                         # HP potential (same units as y)

# HP output gap in percent (since y is 100*log)
gap_hp_pct <- (y - y_hp_star) / 100

# --- 2) Get LW output gap in percent (robust to your object names) ---
if (exists("x_dev")) {
  # This matches your plot(dates, x_dev/100, ...) line
  x_dev <- as.numeric(x_dev)
  n2 <- min(length(x_dev), n)
  gap_lw_pct <- x_dev[seq_len(n2)] / 100
  dates2 <- dates[seq_len(n2)]
  gap_hp_pct2 <- gap_hp_pct[seq_len(n2)]
} else if (exists("y_star")) {
  y_star <- as.numeric(y_star)
  n2 <- min(length(y_star), n)
  gap_lw_pct <- (y[seq_len(n2)] - y_star[seq_len(n2)]) / 100
  dates2 <- dates[seq_len(n2)]
  gap_hp_pct2 <- gap_hp_pct[seq_len(n2)]
} else {
  stop("Need either 'x_dev' (preferred) OR 'y_star' from the LW model to compute LW gap.")
}

# --- 3) Plot: HP vs LW output gap (overlay) ---
yl <- range(c(gap_lw_pct, gap_hp_pct2), na.rm = TRUE)
pad <- 0.10 * diff(yl); if (!is.finite(pad) || pad == 0) pad <- 0.01
yl <- yl + c(-pad, pad)

plot(dates2, gap_lw_pct,
     type="l", lwd=2, ylim=yl,
     main="Output gap: LW vs HP filter",
     xlab="", ylab="Percent")

lines(dates2, gap_hp_pct2, lwd=2, lty=2)
abline(h=0, lty=3, col="grey60")

legend("topright",
       legend=c("LW output gap","HP output gap (lambda=1600)"),
       lty=c(1,2), lwd=2, bty="n")

# --- 4) (Optional) quick numeric check ---
cat("Correlation(LW gap, HP gap) =", round(cor(gap_lw_pct, gap_hp_pct2, use="complete.obs"), 3), "\n")

# Objects you may want later:
# y_hp_star, gap_hp_pct, gap_lw_pct
