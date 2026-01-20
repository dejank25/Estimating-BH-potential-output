# ============================================================
# CB-LW (Currency-board LW-style) model in KFAS
# FULL COPY/PASTE SCRIPT — meaningful potential, gap (%), and time-varying trend growth
#
# What makes this version "meaningful":
# 1) Correct spread scaling (bps -> percentage points) + centering (rho is deviation)
# 2) Gap is anchored via x_gap = x_bar + x_dev (prevents drift)
# 3) Trend growth g is allowed to move (q_g not tiny), so it's not flat
# 4) Output gap is reported in % of potential (x_gap/100)
# 5) Imported inflation handled by adjusting inflation: pi_adj = pi_bih - gamma*pi_ea
#
# INPUTS:
#  - df: columns date, rgdp_sa, cpi_index
#  - cds_croatia.xlsx: date column + CDS/spread column (daily OK)
#
# OUTPUTS (objects in workspace):
#  - results_df with: date, y, y_star, gap_pct, g_ann_pct, rho_pp, pi_tr
#  - plots: y vs y*, gap (%), rho (pp), g annualized (%)
# ============================================================

suppressPackageStartupMessages({
  library(KFAS)
  library(zoo)
  library(dplyr)
  library(eurostat)
  library(rsdmx)
  library(readxl)
})

# ---------------------------
# 0) BiH data
# ---------------------------
stopifnot(exists("df"))
stopifnot(all(c("date","rgdp_sa","cpi_index") %in% names(df)))

df <- df %>%
  mutate(date = as.yearqtr(date)) %>%
  arrange(date) %>%
  mutate(
    y_bih  = 100 * log(as.numeric(rgdp_sa)),
    pi_bih = 100 * (as.numeric(cpi_index) / dplyr::lag(as.numeric(cpi_index), 4) - 1)
  ) %>%
  filter(is.finite(y_bih), is.finite(pi_bih))

# ---------------------------
# 1) EURIBOR 6m quarterly (ECB SDMX)
# ---------------------------
ecb_url <- paste0(
  "https://data-api.ecb.europa.eu/service/data/FM/",
  "Q.U2.EUR.RT.MM.EURIBOR6MD_.HSTA",
  "?format=sdmx-2.1&startPeriod=2000-Q1"
)

eur_sdmx <- readSDMX(ecb_url)
eur_raw  <- as.data.frame(eur_sdmx)

time_col <- if ("obsTime" %in% names(eur_raw)) "obsTime" else names(eur_raw)[grepl("TIME", names(eur_raw), ignore.case=TRUE)][1]
val_col  <- if ("obsValue" %in% names(eur_raw)) "obsValue" else names(eur_raw)[grepl("OBS_VALUE", names(eur_raw), ignore.case=TRUE)][1]

eur_df <- eur_raw %>%
  transmute(
    date_chr = sub("-(Q[1-4])", " \\1", .data[[time_col]]),
    date = as.yearqtr(date_chr, format = "%Y Q%q"),
    euribor6m_ea = as.numeric(.data[[val_col]])
  ) %>%
  filter(is.finite(euribor6m_ea)) %>%
  arrange(date)

# ---------------------------
# 2) EA quarterly real GDP growth (Eurostat)
# ---------------------------
gdp_ea_raw <- get_eurostat("namq_10_gdp", time_format = "date")
geo_pick <- if ("EA20" %in% gdp_ea_raw$geo) "EA20" else if ("EA19" %in% gdp_ea_raw$geo) "EA19" else "EA"

gdp_ea <- gdp_ea_raw %>%
  filter(geo == geo_pick, na_item %in% c("B1GQ","B1G")) %>%
  { if ("s_adj" %in% names(.)) filter(., s_adj %in% c("SCA","SA","CA")) else . } %>%
  { if ("unit" %in% names(.)) filter(., unit %in% c("CLV10_MEUR","CLV05_MEUR","CLV_I15","CLV_I10")) else . } %>%
  arrange(TIME_PERIOD) %>%
  transmute(date = as.yearqtr(TIME_PERIOD), gdp_ea = as.numeric(values)) %>%
  group_by(date) %>% summarise(gdp_ea = mean(gdp_ea, na.rm=TRUE), .groups="drop") %>%
  arrange(date) %>%
  mutate(gdp_ea_gr = 400 * (log(gdp_ea) - log(dplyr::lag(gdp_ea, 1)))) %>%
  filter(is.finite(gdp_ea_gr))

# ---------------------------
# 3) EA HICP monthly -> quarterly y/y inflation (Eurostat)
# ---------------------------
hicp_raw <- get_eurostat("prc_hicp_midx", time_format = "date")

hicp_ea <- hicp_raw %>%
  filter(geo == geo_pick) %>%
  { if ("coicop" %in% names(.)) filter(., coicop %in% c("CP00")) else . } %>%
  arrange(TIME_PERIOD) %>%
  transmute(date_m = as.yearmon(TIME_PERIOD), hicp_ea = as.numeric(values)) %>%
  mutate(date = as.yearqtr(date_m, frac = 1)) %>%
  group_by(date) %>% summarise(hicp_ea = mean(hicp_ea, na.rm=TRUE), .groups="drop") %>%
  arrange(date) %>%
  mutate(pi_ea = 100 * (hicp_ea / dplyr::lag(hicp_ea, 4) - 1)) %>%
  filter(is.finite(pi_ea))

# ---------------------------
# 4) Croatia CDS Excel -> quarterly mean -> spread proxy
# ---------------------------
cds_raw <- read_excel("cds_croatia.xlsx")
nm_l <- tolower(names(cds_raw))

date_col <- names(cds_raw)[which(nm_l %in% c("date","datum","time","day"))[1]]
if (is.na(date_col)) date_col <- names(cds_raw)[1]

spread_candidates <- which(grepl("cds|spread|bps", nm_l))
spread_col <- if (length(spread_candidates) > 0) names(cds_raw)[spread_candidates[1]] else names(cds_raw)[2]

cds_croatia <- cds_raw %>%
  transmute(
    date  = as.Date(.data[[date_col]]),
    spread_raw = as.numeric(.data[[spread_col]])
  ) %>%
  filter(!is.na(date), is.finite(spread_raw)) %>%
  mutate(date = as.yearqtr(date)) %>%
  group_by(date) %>%
  summarise(spread_raw = mean(spread_raw, na.rm=TRUE), .groups="drop") %>%
  arrange(date)

# ---------------------------
# 5) Merge
# ---------------------------
D <- df %>%
  select(date, y_bih, pi_bih) %>%
  inner_join(eur_df, by="date") %>%
  inner_join(gdp_ea %>% select(date, gdp_ea_gr), by="date") %>%
  inner_join(hicp_ea %>% select(date, pi_ea), by="date") %>%
  inner_join(cds_croatia, by="date") %>%
  arrange(date)

stopifnot(nrow(D) >= 12)

# ============================================================
# 6) EA neutral real rate proxy (smoothed EA real short rate)
#    FIX: use low-frequency inflation expectations (MA8),
#    not lagged realized y/y inflation (avoids COVID -7% dip)
# ============================================================

D <- D %>%
  arrange(date) %>%
  mutate(
    # Expected inflation proxy: trailing 8-quarter moving average of y/y HICP
    pi_ea_exp = as.numeric(stats::filter(pi_ea, rep(1/8, 8), sides = 1)),
    
    # EA real short rate proxy (pp): nominal EURIBOR 6m minus expected inflation
    r_ea = euribor6m_ea - pi_ea_exp
  )

cat("UNIT CHECK (EURIBOR 6m): min/mean/max = ",
    min(D$euribor6m_ea, na.rm=TRUE),
    mean(D$euribor6m_ea, na.rm=TRUE),
    max(D$euribor6m_ea, na.rm=TRUE), "\n")

cat("UNIT CHECK (pi_ea_exp MA8): min/mean/max = ",
    min(D$pi_ea_exp, na.rm=TRUE),
    mean(D$pi_ea_exp, na.rm=TRUE),
    max(D$pi_ea_exp, na.rm=TRUE), "\n")

cat("UNIT CHECK (r_ea): min/mean/max = ",
    min(D$r_ea, na.rm=TRUE),
    mean(D$r_ea, na.rm=TRUE),
    max(D$r_ea, na.rm=TRUE), "\n")

# Local-level smoother for r*_EA (EA neutral real rate proxy)
mod_r <- SSModel(D$r_ea ~ SSMtrend(1, Q = matrix(NA)), H = matrix(NA))
fit_r <- fitSSM(
  mod_r,
  inits  = log(c(var(D$r_ea, na.rm=TRUE)/10, var(D$r_ea, na.rm=TRUE)/10)),
  method = "BFGS"
)
kfs_r <- KFS(fit_r$model, smoothing = "state")

D$r_star <- as.numeric(kfs_r$alphahat)


# ---------------------------
# 7) Prepare series + scaling
# ---------------------------
dates  <- D$date
y      <- as.numeric(D$y_bih)
pi     <- as.numeric(D$pi_bih)
pi_for <- as.numeric(D$pi_ea)

# CDS scaling:
# If your CDS is in bps (typical), convert to percentage points:
spr_pp <- as.numeric(D$spread_raw) / 100
# Center so rho is deviation (no permanent drag)
spr_pp_c <- spr_pp - mean(spr_pp, na.rm=TRUE)

# Imported inflation adjustment
gamma <- 0.25
pi_adj <- pi - gamma * pi_for

Tn <- length(y)
cat("Tn =", Tn, "| range:", as.character(min(dates)), "to", as.character(max(dates)), "\n")

Y <- cbind(y, pi_adj, spr_pp_c)
storage.mode(Y) <- "double"
dimnames(Y) <- NULL

# ============================================================
# 8) CB-LW meaningful KFAS model
# States: [ y*, g, x_bar, x_dev, rho, pi_tr ]
# x_gap = x_bar + x_dev  (in 100*log points)
# gap_pct = x_gap / 100
# ============================================================

# Coefficients (moderate, stable)
phi1 <- 0.85     # gap persistence
lam  <- 0.10     # premium -> gap
beta <- 0.20     # gap -> inflation

m <- 6; p <- 3; k <- 6

# Measurement matrix
Zmat <- matrix(0, p, m)
# y = y* + x_bar + x_dev + e_y
Zmat[1,1] <- 1
Zmat[1,3] <- 1
Zmat[1,4] <- 1
# pi_adj = pi_tr + beta*(x_bar + x_dev) + e_pi
Zmat[2,3] <- beta
Zmat[2,4] <- beta
Zmat[2,6] <- 1
# spr_pp_c = rho + e_s
Zmat[3,5] <- 1

# Transition matrix
Tmat <- matrix(0, m, m)
# y* = y* + g + eta_y*
Tmat[1,1] <- 1; Tmat[1,2] <- 1
# g RW (time-varying trend growth)
Tmat[2,2] <- 1
# x_bar RW (slow-moving anchor level)
Tmat[3,3] <- 1
# x_dev AR(1) + premium effect
Tmat[4,4] <- phi1
Tmat[4,5] <- -lam
# rho RW
Tmat[5,5] <- 1
# pi_tr RW
Tmat[6,6] <- 1

Rmat <- diag(1, m)

# Variances to estimate (IMPORTANT: q_g not tiny -> g moves)
par0 <- log(c(
  q_y    = 1e-2,
  q_g    = 5e-3,   # <-- key change: allow time-varying growth
  q_xbar = 1e-4,
  q_xdev = 1e-1,
  q_rho  = 5e-3,
  q_pi   = 5e-3,
  h_y    = 1e-2,
  h_pi   = 1e-1,
  h_s    = 5e-2
))

Q0 <- diag(as.numeric(exp(par0[1:6])), k, k); dimnames(Q0) <- NULL
H0 <- diag(as.numeric(exp(par0[7:9])), p, p); dimnames(H0) <- NULL
storage.mode(Q0) <- "double"; storage.mode(H0) <- "double"

# Diffuse init: y*, g, x_bar, rho, pi_tr diffuse; x_dev finite
a1 <- rep(0, m)

P1inf <- diag(0, m)
P1inf[1,1] <- 1  # y*
P1inf[2,2] <- 1  # g
P1inf[3,3] <- 1  # x_bar
P1inf[5,5] <- 1  # rho
P1inf[6,6] <- 1  # pi_tr

P1 <- diag(0, m)
P1[4,4] <- 1     # x_dev finite

# Build model
model0 <- SSModel(
  Y ~ -1 + SSMcustom(
    Z = unname(Zmat),
    T = unname(Tmat),
    R = unname(Rmat),
    Q = Q0,
    a1 = a1,
    P1 = unname(P1),
    P1inf = unname(P1inf)
  ),
  H = H0
)

ok <- is.SSModel(model0, na.check = TRUE, return.logical = TRUE)
cat("is.SSModel:", ok, "\n")
stopifnot(ok)

updatefn <- function(par, model) {
  qv <- as.numeric(exp(par[1:6]))
  hv <- as.numeric(exp(par[7:9]))
  model$Q[,,1] <- diag(qv, 6, 6)
  model$H[,,1] <- diag(hv, 3, 3)
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
a <- kfs$alphahat

# Extract states
y_star  <- a[,1]
g_trend <- a[,2]
x_bar   <- a[,3]
x_dev   <- a[,4]
rho_pp  <- a[,5]       # premium deviation in percentage points
pi_tr   <- a[,6]

x_gap   <- x_bar + x_dev          # in 100*log points
gap_pct <- x_gap / 100            # approx percent deviation
g_ann_pct <- 4 * g_trend          # annualized %, correct for y=100*log

# Results DF
results_df <- data.frame(
  date      = dates,
  y         = y,
  y_star    = y_star,
  gap_pct   = gap_pct,
  rho_pp    = rho_pp,
  g_ann_pct = g_ann_pct,
  pi_adj    = pi_adj,
  pi_tr     = pi_tr
)

# ---------------------------
# 9) Plots (meaningful units)
# ---------------------------
par(
  mfrow = c(2, 2),
  mar   = c(3.2, 3.4, 2.8, 1.2),  # bottom, left, top, right
  oma   = c(0, 0, 0, 0)
)


# Actual vs Potential GDP — deviation from base period

# Ensure plain vectors (avoid ts / method issues)
dates  <- zoo::as.yearqtr(dates)
y      <- as.numeric(y)
y_star <- as.numeric(y_star)

# Choose base period (pre-COVID benchmark)
base_date <- zoo::as.yearqtr("2019 Q4")
base_idx  <- which(dates == base_date)
if (length(base_idx) != 1) stop("Base date not found uniquely.")

# Deviations from base (percent, since y is 100*log)
y_dev_act <- y      - y[base_idx]
y_dev_pot <- y_star - y_star[base_idx]

# y-axis limits with padding
ylim <- range(c(y_dev_act, y_dev_pot), na.rm = TRUE)
pad  <- 0.05 * diff(ylim)
ylim <- ylim + c(-pad, pad)

# Plot
plot(dates, y_dev_act,
     type = "l",
     lwd  = 2,
     ylim = ylim,
     main = paste(
       "Actual vs Potential GDP",
       "(deviation from 2019 Q4)",
       sep = "\n"
     ),
     font.main = 2,        # <-- makes BOTH lines bold
     ylab = "Deviation (%)",
     xlab = ""
)



lines(dates, y_dev_pot, lwd = 2, lty = 2)
abline(h = 0, col = "grey60", lty = 3)

legend("topleft",
       legend = c("Actual GDP", "Potential GDP"),
       lty    = c(1, 2),
       lwd    = 2,
       bty    = "n")

plot(dates, x_dev/100, type="l", lwd=2, main="Output gap (cyclical, %)", xlab="", ylab="%")
abline(h=0,lty=2)

plot(dates, rho_pp, type="l", lwd=2, main="Premium rho (deviation, pp)", xlab="", ylab="pp")
abline(h=0,lty=2)

plot(dates, g_ann_pct, type="l", lwd=2, main="Trend growth g (annualized, %)", xlab="", ylab="%")
abline(h=0,lty=2)

cat("\nCB-LW meaningful model estimated.\n")
cat("Converged:", fit$optim.out$convergence == 0, "\n")
cat("exp(par):\n")
print(exp(fit$optim.out$par))
