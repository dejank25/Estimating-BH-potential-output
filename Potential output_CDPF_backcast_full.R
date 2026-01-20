## =========================
## Fully runnable PF + Seasonal Adjustment (GDP) code
## =========================

## 0) Packages
pkgs <- c("dplyr", "zoo", "ggplot2", "seasonal", "mFilter", "stringr")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(zoo)
library(ggplot2)
library(seasonal)
library(mFilter)
library(stringr)

## 1) Helper: parse "2008q1" -> yearqtr  (VECTORISED correctly)
to_yearqtr_safe <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[ -]", "", x)
  
  ok <- grepl("^\\d{4}q[1-4]$", x)
  if (any(!ok, na.rm = TRUE)) {
    bad <- unique(x[!ok])
    stop(paste0(
      "date must look like '2008q1'. Bad values: ",
      paste(head(bad, 10), collapse = ", "),
      if (length(bad) > 10) " ..." else ""
    ))
  }
  
  yr <- as.integer(substr(x, 1, 4))
  qq <- as.integer(substr(x, 6, 6))
  as.yearqtr(paste0(yr, " Q", qq), format = "%Y Q%q")
}

## 2) Helper: parse numbers with commas
num_parse <- function(x) {
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x))
  x[x == ""] <- NA
  x <- gsub(",", "", x)
  suppressWarnings(as.numeric(x))
}

## 3) Seasonal adjustment of GDP using X-13 (with fallback if X-13 fails)
seasonal_adjust_gdp <- function(date_yq, y, y_is_log = TRUE) {
  stopifnot(length(date_yq) == length(y))
  
  if (y_is_log) {
    y_level <- exp(y)
  } else {
    y_level <- y
  }
  
  # Build quarterly ts
  start_year <- as.integer(format(min(date_yq), "%Y"))
  start_q    <- as.integer(format(min(date_yq), "%q"))
  y_ts <- ts(as.numeric(y_level), start = c(start_year, start_q), frequency = 4)
  
  # Try X-13; if it fails, fallback to STL decomposition
  y_sa_level <- tryCatch({
    fit <- seas(y_ts)
    as.numeric(final(fit))
  }, error = function(e) {
    message("X-13 failed (seasonal::seas). Falling back to STL: ", e$message)
    stl_fit <- stl(y_ts, s.window = "periodic", robust = TRUE)
    as.numeric(seasadj(stl_fit))
  })
  
  y_sa <- if (y_is_log) log(y_sa_level) else y_sa_level
  list(y_sa = y_sa)
}

## 4) Backcast labor force wf to cover full sample (FIXED)
backcast_wf <- function(date_yq, wf) {
  t_all <- seq_along(date_yq)
  idx <- which(!is.na(wf))
  
  if (length(idx) < 6) stop("Too few wf observations to backcast reliably (need at least ~6 quarters).")
  
  dat <- data.frame(wf = wf[idx], t = t_all[idx])
  fit <- lm(wf ~ t, data = dat)
  
  wf_hat <- as.numeric(predict(fit, newdata = data.frame(t = t_all)))
  
  wf_back <- wf
  wf_back[is.na(wf_back)] <- wf_hat[is.na(wf_back)]
  wf_back
}

## 5) Capital stock via PIM
build_capital_pim <- function(i_real, y_log_sa = NULL, delta_annual = 0.05) {
  if (all(is.na(i_real))) stop("Investment series i is all NA.")
  if (any(i_real <= 0, na.rm = TRUE)) stop("Investment i must be positive (real).")
  
  delta_q <- 1 - (1 - delta_annual)^(1/4)
  
  if (!is.null(y_log_sa)) {
    y_level_sa <- exp(y_log_sa)
    ky_guess <- 2.5
    k0 <- ky_guess * y_level_sa[which(!is.na(y_level_sa))[1]]
  } else {
    i0 <- i_real[which(!is.na(i_real))[1]]
    k0 <- i0 / delta_q
  }
  
  k <- rep(NA_real_, length(i_real))
  k[1] <- k0
  for (t in 2:length(i_real)) {
    it <- i_real[t]
    if (is.na(it)) it <- i_real[t - 1]
    k[t] <- (1 - delta_q) * k[t - 1] + it
  }
  k
}

## 6) PF potential output estimator
estimate_potential_pf <- function(df,
                                  alpha = 0.35,
                                  hp_lambda = 1600,
                                  delta_annual = 0.05,
                                  y_is_log = TRUE,
                                  u_in_percent = FALSE) {
  
  df <- df %>%
    mutate(
      date_yq = to_yearqtr_safe(date),
      y  = num_parse(y),
      u  = num_parse(u),
      i  = num_parse(i),
      wf = num_parse(wf)
    ) %>%
    arrange(date_yq)
  
  if (u_in_percent) df$u <- df$u / 100
  
  ## Seasonal adjustment
  sa <- seasonal_adjust_gdp(df$date_yq, df$y, y_is_log = y_is_log)
  df$y_sa <- sa$y_sa
  
  ## Backcast labor force
  df$wf_backcast <- backcast_wf(df$date_yq, df$wf)
  
  ## Employment proxy
  df$emp <- df$wf_backcast * (1 - df$u)
  
  ## NAWRU via HP on unemployment
  u_hp <- hpfilter(ts(df$u, frequency = 4), freq = hp_lambda)$trend
  df$nawru <- pmin(pmax(as.numeric(u_hp), 0), 1)
  
  df$emp_pot <- df$wf_backcast * (1 - df$nawru)
  
  ## Capital stock
  df$k <- build_capital_pim(df$i, y_log_sa = df$y_sa, delta_annual = delta_annual)
  
  ## Logs
  df$k_log      <- log(df$k)
  df$l_log      <- log(df$emp)
  df$l_pot_log  <- log(df$emp_pot)
  
  ## TFP + trend
  df$tfp <- df$y_sa - alpha * df$k_log - (1 - alpha) * df$l_log
  df$tfp_trend <- as.numeric(hpfilter(ts(df$tfp, frequency = 4), freq = hp_lambda)$trend)
  
  ## Potential output
  df$y_pot <- df$tfp_trend + alpha * df$k_log + (1 - alpha) * df$l_pot_log
  
  ## Output gap
  df$gap_log <- df$y_sa - df$y_pot
  df$gap_pct <- 100 * df$gap_log
  
  ## Date for plotting
  df$date_q <- as.Date(df$date_yq)
  
  df
}

## =========================
## 7) RUN (edit only the file path if needed)
## =========================
df <- read.csv("Data_macro_input1.csv", stringsAsFactors = FALSE)

out <- estimate_potential_pf(
  df,
  alpha = 0.35,
  hp_lambda = 1600,
  delta_annual = 0.05,
  y_is_log = TRUE,
  u_in_percent = FALSE
)

out_clean <- out %>%
  filter(date_yq >= as.yearqtr("2011 Q1"))

## =========================
## 8) PLOTS
## =========================

# Log GDP (SA) vs log potential
ggplot(out_clean, aes(x = date_q)) +
  geom_line(aes(y = y_sa)) +
  geom_line(aes(y = y_pot), linetype = "dashed") +
  labs(
    title = "log Real GDP (seasonally adjusted) vs log Potential Output (PF)",
    x = "Quarter",
    y = "log level"
  ) +
  theme_minimal()

# Output gap
ggplot(out_clean, aes(x = date_q, y = gap_pct)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Output Gap (PF, using seasonally adjusted GDP)",
    x = "Quarter",
    y = "Percent (approx.)"
  ) +
  theme_minimal()

## =========================
## 9) CHECK backcast quickly
## =========================
head(out %>% select(date, wf, wf_backcast), 20)
