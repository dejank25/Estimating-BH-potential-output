# ============================================================
# ROBUSTNESS: Output gaps from 4 methods (single copy/paste)
#   1) Model-based gap (your model): x_dev OR (y - y_star)
#   2) HP filter gap (lambda=1600)
#   3) Band-pass gap (Baxter–King, 6–32 quarters)  [FIXED]
#   4) Production Function (PF) gap (aligned by dates)
#
# REQUIRE:
#   - dates : yearqtr vector
#   - y     : 100*log(real GDP) numeric (same as in the model)
#   - x_dev OR y_star from your model workspace
#   - PF gap optional:
#       out_pf with columns date_yq (yearqtr) and gap_pct (percent)
#       OR vectors gap_pf_pct (percent) and dates_pf (yearqtr)
# ============================================================

suppressPackageStartupMessages({
  library(zoo)
})

if (!requireNamespace("mFilter", quietly = TRUE)) install.packages("mFilter")

# ---------------------------
# 0) Safety + coercion
# ---------------------------
stopifnot(exists("dates"), exists("y"))
dates <- zoo::as.yearqtr(dates)
y     <- as.numeric(y)

# Strip any attributes that can confuse filters
y <- as.numeric(unname(y))
dates <- dates[seq_len(min(length(dates), length(y)))]
y <- y[seq_len(length(dates))]

stopifnot(length(y) >= 20)

# ---------------------------
# 1) MODEL-BASED output gap (percent)
# ---------------------------
if (exists("x_dev")) {
  x_dev <- as.numeric(unname(x_dev))
  n1 <- min(length(x_dev), length(y))
  dates_m  <- dates[seq_len(n1)]
  y_m      <- y[seq_len(n1)]
  gap_model <- x_dev[seq_len(n1)] / 100
} else if (exists("y_star")) {
  y_star <- as.numeric(unname(y_star))
  n1 <- min(length(y_star), length(y))
  dates_m  <- dates[seq_len(n1)]
  y_m      <- y[seq_len(n1)]
  gap_model <- (y_m - y_star[seq_len(n1)]) / 100
} else {
  stop("Need either 'x_dev' OR 'y_star' to compute the model-based gap.")
}

stopifnot(length(y_m) == length(dates_m), length(gap_model) == length(dates_m))

# ---------------------------
# 2) HP filter gap (lambda=1600)
# ---------------------------
hp <- mFilter::hpfilter(y_m, freq = 1600, type = "lambda")
y_hp_star <- as.numeric(hp$trend)
gap_hp <- (y_m - y_hp_star) / 100

# ---------------------------
# 3) Band-pass gap (Baxter–King, 6–32 quarters)  [FIXED]
#    IMPORTANT: do NOT set nfix=0 (causes length mismatch in some versions).
#    Use a "safe" nfix (e.g., 12) OR omit it and accept NA endpoints.
# ---------------------------
bp <- mFilter::bkfilter(as.numeric(y_m), pl = 6, pu = 32, nfix = 6)
gap_bp <- as.numeric(bp$cycle) / 100  # endpoints will be NA; that's normal

# ---------------------------
# 4) PF gap (aligned to dates_m)  [ROBUST]
# ---------------------------
gap_pf <- rep(NA_real_, length(dates_m))

if (exists("out_pf")) {
  stopifnot(all(c("gap_pct") %in% names(out_pf)))
  
  # --- PF dates: accept common column names ---
  if ("date_yq" %in% names(out_pf)) {
    pf_raw <- out_pf$date_yq
  } else if ("date" %in% names(out_pf)) {
    pf_raw <- out_pf$date
  } else if ("Quarter" %in% names(out_pf)) {
    pf_raw <- out_pf$Quarter
  } else {
    stop("out_pf exists but no date column found (expected date_yq/date/Quarter).")
  }
  
  # --- Convert PF dates to yearqtr robustly ---
  pf_dates <- suppressWarnings(zoo::as.yearqtr(pf_raw))
  
  # If conversion failed (NAs), try common string formats:
  if (any(is.na(pf_dates))) {
    s <- tolower(trimws(as.character(pf_raw)))
    s <- gsub("[ -]", "", s)
    
    # handle "2018q1"
    if (all(grepl("^\\d{4}q[1-4]$", s))) {
      yr <- as.integer(substr(s, 1, 4))
      qq <- as.integer(substr(s, 6, 6))
      pf_dates <- zoo::as.yearqtr(paste0(yr, " Q", qq), format = "%Y Q%q")
    } else {
      # handle "2018qtr1" (rare)
      s2 <- gsub("qtr", "q", s)
      pf_dates <- suppressWarnings(zoo::as.yearqtr(s2))
    }
  }
  
  # PF gap (percent -> fraction)
  pf_gap <- as.numeric(out_pf$gap_pct) / 100
  
  # Align by dates
  idx <- match(dates_m, pf_dates)
  gap_pf <- pf_gap[idx]
  
  # Diagnostics
  cat("PF alignment check:\n")
  cat("  PF obs total      =", sum(is.finite(pf_gap)), "\n")
  cat("  PF dates finite   =", sum(!is.na(pf_dates)), "\n")
  cat("  PF matches on sample =", sum(is.finite(gap_pf)), "out of", length(gap_pf), "\n")
  
} else if (exists("gap_pf_pct") && exists("dates_pf")) {
  
  pf_dates <- zoo::as.yearqtr(dates_pf)
  pf_gap   <- as.numeric(gap_pf_pct) / 100
  idx <- match(dates_m, pf_dates)
  gap_pf <- pf_gap[idx]
  
  cat("PF alignment check (vectors): matches =", sum(is.finite(gap_pf)), "\n")
  
} else {
  message("PF gap not found. Create 'out_pf' (with date_yq + gap_pct) or provide gap_pf_pct + dates_pf.")
}

# ============================================================
# ADD PF gap to the overlay (works with out_clean / out from PF script)
# REQUIRE:
#   - dates_m (yearqtr) + gap_model, gap_hp, gap_bp already computed
#   - PF results object exists as: out_clean OR out OR out_pf
#     and contains:
#       * gap_pct  (PF output gap in percent)
#       * date_yq OR date_q OR date (some date representation)
# ============================================================

suppressPackageStartupMessages(library(zoo))

# --- 1) Pick the PF object that exists ---
PF <- NULL
if (exists("out_clean")) PF <- out_clean
if (is.null(PF) && exists("out")) PF <- out
if (is.null(PF) && exists("out_pf")) PF <- out_pf
stopifnot(!is.null(PF))

stopifnot("gap_pct" %in% names(PF))
pf_gap <- as.numeric(PF$gap_pct) / 100  # convert percent -> fraction (same units as other gaps)

# --- 2) Get PF dates robustly and convert to yearqtr ---
pf_dates <- NULL

if ("date_yq" %in% names(PF)) {
  # sometimes already yearqtr, sometimes string; handle both
  pf_dates <- suppressWarnings(as.yearqtr(PF$date_yq))
}

if (is.null(pf_dates) || all(is.na(pf_dates))) {
  if ("date_q" %in% names(PF)) {
    # date_q is Date -> convert to yearqtr
    pf_dates <- as.yearqtr(as.Date(PF$date_q))
  } else if ("date" %in% names(PF)) {
    # could be "2018q1"/"2018 Q1"/Date; try yearqtr conversion
    pf_dates <- suppressWarnings(as.yearqtr(PF$date))
    if (all(is.na(pf_dates))) {
      s <- tolower(trimws(as.character(PF$date)))
      s <- gsub("[ -]", "", s)
      if (all(grepl("^\\d{4}q[1-4]$", s))) {
        yr <- as.integer(substr(s, 1, 4))
        qq <- as.integer(substr(s, 6, 6))
        pf_dates <- as.yearqtr(paste0(yr, " Q", qq), format = "%Y Q%q")
      }
    }
  }
}

stopifnot(!is.null(pf_dates))

# --- 3) Align PF gap to the model sample dates_m ---
idx_pf <- match(dates_m, pf_dates)
gap_pf <- pf_gap[idx_pf]

cat("PF alignment: finite PF points on model sample =",
    sum(is.finite(gap_pf)), "out of", length(gap_pf), "\n")

# --- 4) Re-plot overlay WITH PF ---
stack <- c(gap_model, gap_hp, gap_bp, gap_pf)
yl <- range(stack, na.rm = TRUE)
pad <- 0.10 * diff(yl); if (!is.finite(pad) || pad == 0) pad <- 0.01
yl <- yl + c(-pad, pad)

plot(dates_m, gap_model, type="l", lwd=2, ylim=yl,
     main="Output gap comparison: Model vs HP vs Band-pass vs PF",
     xlab="", ylab="Percent")

lines(dates_m, gap_hp, lwd=2, lty=2)
lines(dates_m, gap_bp, lwd=2, lty=3)

# PF line (draw if at least 3 points matched)
if (sum(is.finite(gap_pf)) >= 3) {
  lines(dates_m, gap_pf, lwd=2, lty=4)
} else {
  cat("PF not drawn because too few matched points. Check PF dates vs dates_m.\n")
}

abline(h=0, lty=5, col="grey60")

leg  <- c("Model-based gap", "HP gap (lambda=1600)", "Band-pass gap (BK 6–32q)")
ltys <- c(1, 2, 3)
if (sum(is.finite(gap_pf)) >= 3) {
  leg  <- c(leg, "PF gap")
  ltys <- c(ltys, 4)
}

legend("bottomright", legend=leg, lty=ltys, lwd=2, bty="n")

# ---------------------------
# 6) Quick correlations (optional)
# ---------------------------
cat("Correlations with Model gap (complete obs):\n")
cat("  corr(Model, HP) =", round(cor(gap_model, gap_hp, use="complete.obs"), 3), "\n")
cat("  corr(Model, BP) =", round(cor(gap_model, gap_bp, use="complete.obs"), 3), "\n")
if (any(is.finite(gap_pf))) {
  cat("  corr(Model, PF) =", round(cor(gap_model, gap_pf, use="complete.obs"), 3), "\n")
}
