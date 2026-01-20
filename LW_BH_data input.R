# ============================================================
# BiH LW input preparation
# - CPI is an index (no SA needed)
# - Real GDP is raw (needs seasonal adjustment -> SA -> log)
# Output: vectors y, pi, i (and optionally r) aligned quarterly
# ============================================================

suppressPackageStartupMessages({
  library(zoo)
  library(dplyr)
  library(seasonal)   # install.packages("seasonal")
  library(x13binary)  # install.packages("x13binary")
})

# ---------------------------
# 1) LOAD DATA (edit paths)
# ---------------------------
rgdp_raw <- read.csv("rgdp.csv", stringsAsFactors = FALSE)     # date, rgdp_raw
cpi_idx  <- read.csv("cpi.csv", stringsAsFactors = FALSE)         # date, cpi_index
eur      <- read.csv("euribor6m.csv", stringsAsFactors = FALSE)       # date, euribor6m

# Dates must be like "2000 Q1" or "2000Q1" or "2000-01-01" etc.
# Best: "YYYY-Qx". If you have another format, adjust parsing here.
rgdp_raw$date <- as.yearqtr(rgdp_raw$date)
cpi_idx$date  <- as.yearqtr(cpi_idx$date)
eur$date      <- as.yearqtr(eur$date)

# Ensure numeric
rgdp_raw$rgdp_raw   <- as.numeric(rgdp_raw$rgdp)
cpi_idx$cpi_index   <- as.numeric(cpi_idx$cpi)
eur$euribor6m       <- as.numeric(eur$euribor6m)

# ---------------------------
# 2) SEASONALLY ADJUST RGDP (X-13)
# ---------------------------
# X-13 likes "ts" with frequency=4 (quarterly).
# Build ts object spanning the available sample.
start_y <- as.integer(format(min(rgdp_raw$date), "%Y"))
start_q <- as.integer(cycle(as.ts(rgdp_raw$date))[1])  # robust-ish, but see below

# More robust start quarter:
start_q <- as.integer(round((as.numeric(min(rgdp_raw$date)) - floor(as.numeric(min(rgdp_raw$date)))) * 4 + 1))

rgdp_ts <- ts(rgdp_raw$rgdp_raw, start = c(start_y, start_q), frequency = 4)

# Seasonal adjustment (auto spec). You can tweak: transform.function="log" if needed.
sa_fit <- seas(rgdp_ts)

rgdp_sa_ts <- final(sa_fit)  # seasonally adjusted series (same freq/length)

# Convert back to a data frame aligned with original dates
rgdp_sa <- data.frame(
  date    = rgdp_raw$date,
  rgdp_sa = as.numeric(rgdp_sa_ts)
)

# ---------------------------
# 3) MERGE QUARTERLY DATA
# ---------------------------
df <- rgdp_sa %>%
  inner_join(cpi_idx, by = "date") %>%
  inner_join(eur, by = "date") %>%
  arrange(date)

# ---------------------------
# 4) TRANSFORMATIONS FOR LW INPUT
# ---------------------------

# (A) Output: y = log(seasonally adjusted real GDP level)
df$y <- 100*log(df$rgdp_sa)

# (B) Inflation from CPI index:
# Choose ONE definition and keep it throughout the model.
# Recommended for stability: y/y (q/q-4)
df$pi <- 100 * (df$cpi_index / dplyr::lag(df$cpi_index, 4) - 1)

# (C) Nominal rate: EURIBOR 6m (already in percent)
df$i <- df$euribor6m

# (D) Expected inflation proxy (if you need it): lagged inflation
df$pi_exp <- dplyr::lag(df$pi, 1)

# (E) Real rate proxy (used in IS): r = i - expected inflation
df$r <- df$i - df$pi_exp

# ---------------------------
# 5) CLEAN SAMPLE (drop NA created by lags)
# ---------------------------
df <- df %>%
  filter(!is.na(y), !is.na(pi), !is.na(i), !is.na(r))

# Final vectors for the LW script
y     <- df$y
pi    <- df$pi
i     <- df$i
r     <- df$r
dates <- df$date


# Limit to 2017Q1 on wards
DF<-df
df_2017 <- df %>%
  mutate(date = as.yearqtr(date)) %>%   # safe even if already yearqtr
  filter(date >= as.yearqtr("2017 Q1")) %>%
  arrange(date)
df<-df_2017
# ---------------------------
# 6) QUICK DIAGNOSTICS
# ---------------------------
print("BiH LW inputs ready: y (log SA rgdp), pi (CPI y/y), i (euribor), r (real proxy).")
print(head(df, 8))

plot(dates, y,  type="l", main="log real GDP (SA)", xlab="", ylab="")
plot(dates, pi, type="l", main="Inflation (CPI y/y)", xlab="", ylab="")
plot(dates, r,  type="l", main="Real rate proxy: i - pi_exp", xlab="", ylab="")
