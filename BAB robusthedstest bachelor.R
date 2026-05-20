
# ==============================================================================
# STEP 5.1: ROBUSTHEDSTEST — KORRELATIONS- OG VOLATILITETSVINDUER
# ==============================================================================
# Formål:
#   Tester om BAB-strategiens resultater er robuste over for centrale valg i
#   betaestimationen.
#
# Baseline i hovedanalysen:
#   - Korrelationsvindue: 260 uger
#   - Volatilitetsvindue: 52 uger
#
# Robusthedstest:
#   - Korrelationsvinduer: 52, 156 og 260 uger
#   - Volatilitetsvinduer : 26, 52 og 104 uger
#
# Vigtigt:
#   Performance beregnes på månedlige, kompoundede BAB-afkast, så Sharpe-ratioer
#   og alphaer er direkte sammenlignelige med hovedtabellen i STEP 4.7.
#
# Output:
#   - robusthed_corr_vol_grid.csv
#   - robusthed_corr_vol_grid_summary.rds
#   - robusthed_corr_vol_sharpe_matrix.csv
#   - robusthed_corr_vol_heatmap.png
#   - robusthed_corr_vol_sharpe_bar.png
#   - bab_robust_corrXX_volYY.rds for hver specifikation
# ==============================================================================

cat("\n=== STEP 5.1: Robusthedstest — korrelations- og volatilitetsvinduer ===\n")

library(data.table)
library(sandwich)
library(lmtest)
library(ggplot2)

# ------------------------------------------------------------------------------
# 0) Sikkerhed: find korrekt inputdata
# ------------------------------------------------------------------------------

# Robusthedstesten skal bruge data_weekly_clean, dvs. ugentlig data uden
# allerede beregnede betaer fra baseline.
if (!exists("data_weekly_clean")) {
  if (file.exists("data_weekly_clean.rds")) {
    data_weekly_clean <- readRDS("data_weekly_clean.rds")
    cat("Indlæst data_weekly_clean.rds\n")
  } else if (exists("data_weekly")) {
    data_weekly_clean <- copy(data_weekly)
    cat("Bruger data_weekly som fallback\n")
  } else {
    stop("Kan ikke finde data_weekly_clean eller data_weekly.")
  }
}

setDT(data_weekly_clean)

required_cols <- c(
  "id", "iso_year", "iso_week", "date",
  "ret_exc_wk", "mkt_exc_wk", "me"
)

missing_cols <- setdiff(required_cols, names(data_weekly_clean))

if (length(missing_cols) > 0) {
  stop(
    "data_weekly_clean mangler følgende nødvendige kolonner: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 1) Hjælpefunktion: månedlig performance-statistik
# ------------------------------------------------------------------------------

calc_bab_performance_monthly <- function(bab_dt,
                                         label,
                                         corr_window_weeks,
                                         vol_window_weeks) {
  
  bab_dt <- copy(bab_dt)
  
  bab_dt <- bab_dt[
    !is.na(r_bab) &
      !is.na(r_mkt) &
      is.finite(r_bab) &
      is.finite(r_mkt)
  ]
  
  if (nrow(bab_dt) < 52) {
    warning("Færre end 52 ugentlige observationer i ", label)
    
    return(data.table(
      test = label,
      corr_window_weeks = corr_window_weeks,
      vol_window_weeks = vol_window_weeks,
      corr_years = round(corr_window_weeks / 52, 2),
      vol_years = round(vol_window_weeks / 52, 2),
      n_weeks = nrow(bab_dt),
      n_months = NA_integer_,
      start_date = as.IDate(NA),
      end_date = as.IDate(NA),
      mean_monthly_pct = NA_real_,
      nw_t_mean = NA_real_,
      nw_p_mean = NA_real_,
      ann_return_pct = NA_real_,
      ann_vol_pct = NA_real_,
      sharpe = NA_real_,
      max_drawdown_pct = NA_real_,
      capm_alpha_m_pct = NA_real_,
      capm_alpha_t = NA_real_,
      capm_alpha_p = NA_real_,
      capm_beta = NA_real_,
      capm_r2 = NA_real_,
      beta_low_avg = NA_real_,
      beta_high_avg = NA_real_,
      beta_spread_avg = NA_real_,
      n_low_avg = NA_real_,
      n_high_avg = NA_real_
    ))
  }
  
  # --------------------------------------------------------------------------
  # Aggreger ugentlige BAB-afkast til månedlige, kompoundede afkast
  # --------------------------------------------------------------------------
  
  bab_dt[, `:=`(
    yr  = as.integer(format(date, "%Y")),
    mth = as.integer(format(date, "%m"))
  )]
  
  bab_mth <- bab_dt[, .(
    r_bab_m = prod(1 + r_bab, na.rm = TRUE) - 1,
    r_mkt_m = prod(1 + r_mkt, na.rm = TRUE) - 1,
    
    # Deskriptive porteføljemål
    beta_low  = mean(beta_low,  na.rm = TRUE),
    beta_high = mean(beta_high, na.rm = TRUE),
    n_low     = mean(n_low,     na.rm = TRUE),
    n_high    = mean(n_high,    na.rm = TRUE),
    
    date = max(date, na.rm = TRUE)
  ), by = .(yr, mth)]
  
  setorder(bab_mth, yr, mth)
  
  bab_mth <- bab_mth[
    !is.na(r_bab_m) &
      !is.na(r_mkt_m) &
      is.finite(r_bab_m) &
      is.finite(r_mkt_m)
  ]
  
  if (nrow(bab_mth) < 24) {
    warning("Færre end 24 månedlige observationer i ", label)
  }
  
  # --------------------------------------------------------------------------
  # Mean return-test med Newey-West
  # --------------------------------------------------------------------------
  
  mod_int <- lm(r_bab_m ~ 1, data = bab_mth)
  
  nw_int <- coeftest(
    mod_int,
    vcov = NeweyWest(mod_int, lag = 3, prewhite = FALSE)
  )
  
  # --------------------------------------------------------------------------
  # CAPM-alpha med Newey-West
  # --------------------------------------------------------------------------
  
  mod_capm <- lm(r_bab_m ~ r_mkt_m, data = bab_mth)
  
  nw_capm <- coeftest(
    mod_capm,
    vcov = NeweyWest(mod_capm, lag = 3, prewhite = FALSE)
  )
  
  # --------------------------------------------------------------------------
  # Annualiserede nøgletal baseret på månedlige afkast
  # --------------------------------------------------------------------------
  
  r <- bab_mth$r_bab_m
  
  mu    <- mean(r, na.rm = TRUE)
  sigma <- sd(r, na.rm = TRUE)
  
  ann_ret    <- mu * 12
  ann_vol    <- sigma * sqrt(12)
  ann_sharpe <- ann_ret / ann_vol
  
  cum_ret  <- cumprod(1 + r)
  peak     <- cummax(cum_ret)
  drawdown <- (cum_ret - peak) / peak
  max_dd   <- min(drawdown, na.rm = TRUE)
  
  data.table(
    test = label,
    corr_window_weeks = corr_window_weeks,
    vol_window_weeks = vol_window_weeks,
    corr_years = round(corr_window_weeks / 52, 2),
    vol_years = round(vol_window_weeks / 52, 2),
    
    n_weeks = nrow(bab_dt),
    n_months = nrow(bab_mth),
    start_date = min(bab_mth$date, na.rm = TRUE),
    end_date = max(bab_mth$date, na.rm = TRUE),
    
    mean_monthly_pct = mu * 100,
    nw_t_mean = nw_int[1, 3],
    nw_p_mean = nw_int[1, 4],
    
    ann_return_pct = ann_ret * 100,
    ann_vol_pct = ann_vol * 100,
    sharpe = ann_sharpe,
    max_drawdown_pct = max_dd * 100,
    
    capm_alpha_m_pct = nw_capm[1, 1] * 100,
    capm_alpha_t = nw_capm[1, 3],
    capm_alpha_p = nw_capm[1, 4],
    capm_beta = nw_capm[2, 1],
    capm_r2 = summary(mod_capm)$r.squared,
    
    beta_low_avg = mean(bab_mth$beta_low, na.rm = TRUE),
    beta_high_avg = mean(bab_mth$beta_high, na.rm = TRUE),
    beta_spread_avg = mean(bab_mth$beta_high - bab_mth$beta_low, na.rm = TRUE),
    n_low_avg = mean(bab_mth$n_low, na.rm = TRUE),
    n_high_avg = mean(bab_mth$n_high, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# 2) Hjælpefunktion: konstruér BAB med valgt correlation og volatility window
# ------------------------------------------------------------------------------

run_bab_corr_vol_test <- function(data_weekly_input,
                                  corr_window_weeks,
                                  vol_window_weeks,
                                  min_obs_corr,
                                  min_obs_vol,
                                  shrinkage_weight = 0.6,
                                  label = NULL,
                                  save_bab = TRUE) {
  
  if (is.null(label)) {
    label <- paste0(
      "Corr ", corr_window_weeks,
      " / Vol ", vol_window_weeks
    )
  }
  
  cat("\n------------------------------------------------------------\n")
  cat("Robusthedstest:", label, "\n")
  cat("Korrelationsvindue:", corr_window_weeks, "uger",
      "| Minimum obs.:", min_obs_corr, "\n")
  cat("Volatilitetsvindue :", vol_window_weeks, "uger",
      "| Minimum obs.:", min_obs_vol, "\n")
  cat("Shrinkage-vægt     :", shrinkage_weight, "\n")
  cat("------------------------------------------------------------\n")
  
  dt <- copy(data_weekly_input)
  
  dt <- dt[
    !is.na(id) &
      !is.na(iso_year) &
      !is.na(iso_week) &
      !is.na(date) &
      !is.na(ret_exc_wk) &
      !is.na(mkt_exc_wk) &
      !is.na(me) &
      me > 0 &
      is.finite(ret_exc_wk) &
      is.finite(mkt_exc_wk) &
      is.finite(me)
  ]
  
  setorder(dt, id, iso_year, iso_week)
  
  # --------------------------------------------------------------------------
  # 2.1 Rolling korrelation mellem aktie og marked
  # --------------------------------------------------------------------------
  
  cat("Beregner rolling korrelation...\n")
  
  dt[, rho_alt := frollapply(
    seq_len(.N),
    n = corr_window_weeks,
    FUN = function(idx) {
      xi <- ret_exc_wk[idx]
      xm <- mkt_exc_wk[idx]
      
      ok <- !is.na(xi) & !is.na(xm) & is.finite(xi) & is.finite(xm)
      
      if (sum(ok) < min_obs_corr) return(NA_real_)
      
      out <- suppressWarnings(cor(xi[ok], xm[ok]))
      
      if (!is.finite(out)) return(NA_real_)
      
      out
    },
    fill = NA_real_
  ), by = id]
  
  cat(sprintf(
    "Gyldig rho_alt: %.1f%%\n",
    100 * mean(!is.na(dt$rho_alt))
  ))
  
  # --------------------------------------------------------------------------
  # 2.2 Rolling volatilitet for aktier
  # --------------------------------------------------------------------------
  
  cat("Beregner rolling volatilitet for aktier...\n")
  
  dt[, sigma_i_alt := frollapply(
    ret_exc_wk,
    n = vol_window_weeks,
    FUN = function(x) {
      ok <- !is.na(x) & is.finite(x)
      if (sum(ok) < min_obs_vol) return(NA_real_)
      out <- sd(x[ok])
      if (!is.finite(out) || out <= 0) return(NA_real_)
      out
    },
    fill = NA_real_
  ), by = id]
  
  cat(sprintf(
    "Gyldig sigma_i_alt: %.1f%%\n",
    100 * mean(!is.na(dt$sigma_i_alt))
  ))
  
  # --------------------------------------------------------------------------
  # 2.3 Rolling volatilitet for markedet
  # --------------------------------------------------------------------------
  
  cat("Beregner rolling volatilitet for markedet...\n")
  
  mkt_sigma_alt <- dt[, .(
    mkt_exc_wk = first(mkt_exc_wk)
  ), by = .(iso_year, iso_week)]
  
  setorder(mkt_sigma_alt, iso_year, iso_week)
  
  mkt_sigma_alt[, sigma_m_alt := frollapply(
    mkt_exc_wk,
    n = vol_window_weeks,
    FUN = function(x) {
      ok <- !is.na(x) & is.finite(x)
      if (sum(ok) < min_obs_vol) return(NA_real_)
      out <- sd(x[ok])
      if (!is.finite(out) || out <= 0) return(NA_real_)
      out
    },
    fill = NA_real_
  )]
  
  dt <- mkt_sigma_alt[, .(
    iso_year,
    iso_week,
    sigma_m_alt
  )][
    dt,
    on = .(iso_year, iso_week)
  ]
  
  setorder(dt, id, iso_year, iso_week)
  
  # --------------------------------------------------------------------------
  # 2.4 Beregn raw beta, shrinkage og lag
  # --------------------------------------------------------------------------
  
  cat("Beregner beta, shrinkage og lag...\n")
  
  dt[, beta_raw_alt := rho_alt * (sigma_i_alt / sigma_m_alt)]
  
  dt[
    !is.finite(beta_raw_alt),
    beta_raw_alt := NA_real_
  ]
  
  dt[, beta_alt := shrinkage_weight * beta_raw_alt +
       (1 - shrinkage_weight) * 1]
  
  dt[, beta_lag_alt := shift(
    beta_alt,
    n = 1,
    type = "lag"
  ), by = id]
  
  cat(sprintf(
    "Gyldig beta_lag_alt: %.1f%%\n",
    100 * mean(!is.na(dt$beta_lag_alt))
  ))
  
  # --------------------------------------------------------------------------
  # 2.5 Konstruér BAB-portefølje med rank-vægte
  # --------------------------------------------------------------------------
  
  data_bab_alt <- dt[
    !is.na(beta_lag_alt) &
      !is.na(ret_exc_wk) &
      !is.na(mkt_exc_wk) &
      !is.na(me) &
      is.finite(beta_lag_alt) &
      is.finite(ret_exc_wk) &
      is.finite(mkt_exc_wk) &
      is.finite(me) &
      me > 0
  ]
  
  cat(sprintf(
    "Gyldige aktie-uge-observationer: %d\n",
    nrow(data_bab_alt)
  ))
  
  cat(sprintf(
    "Antal uger: %d\n",
    uniqueN(data_bab_alt[, paste(iso_year, iso_week)])
  ))
  
  setorder(data_bab_alt, iso_year, iso_week, beta_lag_alt)
  
  data_bab_alt[, rnk := rank(
    beta_lag_alt,
    ties.method = "average"
  ), by = .(iso_year, iso_week)]
  
  data_bab_alt[, z := rnk - mean(rnk), by = .(iso_year, iso_week)]
  
  data_bab_alt[, beta_med := median(
    beta_lag_alt,
    na.rm = TRUE
  ), by = .(iso_year, iso_week)]
  
  data_bab_alt[, portfolio := fifelse(
    beta_lag_alt <= beta_med,
    "low",
    "high"
  )]
  
  data_bab_alt[portfolio == "low",
               w := -z / sum(-z, na.rm = TRUE),
               by = .(iso_year, iso_week)]
  
  data_bab_alt[portfolio == "high",
               w := z / sum(z, na.rm = TRUE),
               by = .(iso_year, iso_week)]
  
  data_bab_alt <- data_bab_alt[
    !is.na(w) &
      is.finite(w)
  ]
  
  low_alt <- data_bab_alt[portfolio == "low", .(
    r_low    = sum(w * ret_exc_wk,   na.rm = TRUE),
    beta_low = sum(w * beta_lag_alt, na.rm = TRUE),
    n_low    = .N
  ), by = .(iso_year, iso_week)]
  
  high_alt <- data_bab_alt[portfolio == "high", .(
    r_high    = sum(w * ret_exc_wk,   na.rm = TRUE),
    beta_high = sum(w * beta_lag_alt, na.rm = TRUE),
    n_high    = .N
  ), by = .(iso_year, iso_week)]
  
  mkt_alt <- data_bab_alt[, .(
    r_mkt = first(mkt_exc_wk),
    date  = max(date, na.rm = TRUE)
  ), by = .(iso_year, iso_week)]
  
  bab_alt <- merge(
    low_alt,
    high_alt,
    by = c("iso_year", "iso_week")
  )
  
  bab_alt <- merge(
    bab_alt,
    mkt_alt,
    by = c("iso_year", "iso_week")
  )
  
  setorder(bab_alt, iso_year, iso_week)
  
  # Undgå division med nul eller ugyldige porteføljebetaer
  bab_alt <- bab_alt[
    !is.na(beta_low) &
      !is.na(beta_high) &
      is.finite(beta_low) &
      is.finite(beta_high) &
      abs(beta_low) > 1e-6 &
      abs(beta_high) > 1e-6
  ]
  
  bab_alt[, r_bab := (1 / beta_low) * r_low -
            (1 / beta_high) * r_high]
  
  bab_alt <- bab_alt[
    !is.na(r_bab) &
      !is.na(r_mkt) &
      is.finite(r_bab) &
      is.finite(r_mkt)
  ]
  
  cat(sprintf(
    "BAB-portefølje konstrueret: %d uger\n",
    nrow(bab_alt)
  ))
  
  # --------------------------------------------------------------------------
  # 2.6 Månedlig performance
  # --------------------------------------------------------------------------
  
  out_summary <- calc_bab_performance_monthly(
    bab_dt = bab_alt,
    label = label,
    corr_window_weeks = corr_window_weeks,
    vol_window_weeks = vol_window_weeks
  )
  
  cat(sprintf("Månedligt gns. afkast : %+.4f%%\n",
              out_summary$mean_monthly_pct))
  cat(sprintf("Ann. afkast           : %+.2f%%\n",
              out_summary$ann_return_pct))
  cat(sprintf("Ann. vol.             :  %.2f%%\n",
              out_summary$ann_vol_pct))
  cat(sprintf("Sharpe-ratio          :  %.3f\n",
              out_summary$sharpe))
  cat(sprintf("CAPM-alpha, måned     : %+.4f%%, t = %.3f\n",
              out_summary$capm_alpha_m_pct,
              out_summary$capm_alpha_t))
  
  # --------------------------------------------------------------------------
  # 2.7 Gem individuel BAB-serie
  # --------------------------------------------------------------------------
  
  if (isTRUE(save_bab)) {
    
    file_stub <- paste0(
      "bab_robust_corr",
      corr_window_weeks,
      "_vol",
      vol_window_weeks,
      ".rds"
    )
    
    saveRDS(bab_alt, file_stub)
    
    cat("Gemt:", file_stub, "\n")
  }
  
  list(
    summary = out_summary,
    bab     = bab_alt,
    data    = data_bab_alt
  )
}

# ------------------------------------------------------------------------------
# 3) Definér grid og minimumskrav
# ------------------------------------------------------------------------------

corr_windows <- c(52, 156, 260)
vol_windows  <- c(26, 52, 104)

# Minimum observationer:
# - Korrelation: matcher jeres tidligere logik og baseline:
#   52 uger  -> 26 obs.
#   156 uger -> 78 obs.
#   260 uger -> 156 obs.
# - Volatilitet:
#   26 uger  -> 13 obs.
#   52 uger  -> 26 obs.
#   104 uger -> 52 obs.

min_obs_corr_map <- data.table(
  corr_window_weeks = c(52, 156, 260),
  min_obs_corr      = c(26, 78, 156)
)

min_obs_vol_map <- data.table(
  vol_window_weeks = c(26, 52, 104),
  min_obs_vol      = c(13, 26, 52)
)

# ------------------------------------------------------------------------------
# 4) Kør grid: 3 korrelationsvinduer x 3 volatilitetsvinduer
# ------------------------------------------------------------------------------

robust_grid_results <- list()
robust_grid_bab     <- list()

counter <- 0L

for (cw in corr_windows) {
  for (vw in vol_windows) {
    
    counter <- counter + 1L
    
    label_i <- paste0("Corr ", cw, " / Vol ", vw)
    
    min_corr_i <- min_obs_corr_map[
      corr_window_weeks == cw,
      min_obs_corr
    ]
    
    min_vol_i <- min_obs_vol_map[
      vol_window_weeks == vw,
      min_obs_vol
    ]
    
    res_i <- run_bab_corr_vol_test(
      data_weekly_input = data_weekly_clean,
      corr_window_weeks = cw,
      vol_window_weeks  = vw,
      min_obs_corr      = min_corr_i,
      min_obs_vol       = min_vol_i,
      shrinkage_weight  = 0.6,
      label             = label_i,
      save_bab          = TRUE
    )
    
    robust_grid_results[[counter]] <- res_i$summary
    robust_grid_bab[[label_i]]     <- res_i$bab
  }
}

robust_corr_vol_grid <- rbindlist(
  robust_grid_results,
  fill = TRUE
)

setorder(
  robust_corr_vol_grid,
  corr_window_weeks,
  vol_window_weeks
)

cat("\n=== Robusthedsmatrix: korrelations- og volatilitetsvinduer ===\n")
print(robust_corr_vol_grid)

# ------------------------------------------------------------------------------
# 5) Gem samlet output
# ------------------------------------------------------------------------------

fwrite(
  robust_corr_vol_grid,
  "robusthed_corr_vol_grid.csv"
)

saveRDS(
  robust_corr_vol_grid,
  "robusthed_corr_vol_grid_summary.rds"
)

saveRDS(
  robust_grid_bab,
  "robusthed_corr_vol_grid_bab_series.rds"
)

cat("\nGemt:\n")
cat("  -> robusthed_corr_vol_grid.csv\n")
cat("  -> robusthed_corr_vol_grid_summary.rds\n")
cat("  -> robusthed_corr_vol_grid_bab_series.rds\n")

# ------------------------------------------------------------------------------
# 6) Lav Sharpe-matrix til tabel
# ------------------------------------------------------------------------------

sharpe_matrix <- dcast(
  robust_corr_vol_grid,
  corr_window_weeks ~ vol_window_weeks,
  value.var = "sharpe"
)

setnames(
  sharpe_matrix,
  old = as.character(vol_windows),
  new = paste0("Vol_", vol_windows)
)

cat("\nSharpe-matrix, månedlige afkast:\n")
print(sharpe_matrix)

fwrite(
  sharpe_matrix,
  "robusthed_corr_vol_sharpe_matrix.csv"
)

cat("Gemt: robusthed_corr_vol_sharpe_matrix.csv\n")

# ------------------------------------------------------------------------------
# 7) Heatmap: Sharpe-ratio
# ------------------------------------------------------------------------------

plot_grid <- copy(robust_corr_vol_grid)

plot_grid[, corr_label := paste0(corr_window_weeks, " uger")]
plot_grid[, vol_label  := paste0(vol_window_weeks,  " uger")]

plot_grid[, corr_label := factor(
  corr_label,
  levels = paste0(corr_windows, " uger")
)]

plot_grid[, vol_label := factor(
  vol_label,
  levels = paste0(vol_windows, " uger")
)]

p_heat <- ggplot(
  plot_grid,
  aes(
    x = vol_label,
    y = corr_label,
    fill = sharpe
  )
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(
    aes(label = sprintf("%.3f", sharpe)),
    size = 4
  ) +
  labs(
    title = "Robusthedstest: Sharpe-ratio ved alternative betaestimationsvinduer",
    subtitle = "Sharpe-ratio beregnes på månedlige, kompoundede BAB-afkast og annualiseres med √12",
    x = "Volatilitetsvindue",
    y = "Korrelationsvindue",
    fill = "Sharpe-ratio"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    panel.grid = element_blank()
  )

ggsave(
  filename = "robusthed_corr_vol_heatmap.png",
  plot = p_heat,
  width = 9,
  height = 6,
  dpi = 300
)

cat("Gemt: robusthed_corr_vol_heatmap.png\n")

# ------------------------------------------------------------------------------
# 8) Barplot: Sharpe-ratio på tværs af specifikationer
# ------------------------------------------------------------------------------

plot_grid[, test_label := paste0(
  "Corr ", corr_window_weeks,
  " / Vol ", vol_window_weeks
)]

plot_grid[, test_label := factor(
  test_label,
  levels = plot_grid[
    order(corr_window_weeks, vol_window_weeks),
    test_label
  ]
)]

p_bar <- ggplot(
  plot_grid,
  aes(
    x = test_label,
    y = sharpe
  )
) +
  geom_col(width = 0.65, fill = "steelblue") +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  labs(
    title = "Robusthedstest: BAB Sharpe-ratio ved alternative betaestimationsvinduer",
    subtitle = "Sharpe-ratio beregnes på månedlige, kompoundede afkast",
    x = NULL,
    y = "Annualiseret Sharpe-ratio"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(
  filename = "robusthed_corr_vol_sharpe_bar.png",
  plot = p_bar,
  width = 11,
  height = 6,
  dpi = 300
)

cat("Gemt: robusthed_corr_vol_sharpe_bar.png\n")

# ------------------------------------------------------------------------------
# 9) Automatisk opsummering
# ------------------------------------------------------------------------------

best_spec <- robust_corr_vol_grid[which.max(sharpe)]
worst_spec <- robust_corr_vol_grid[which.min(sharpe)]

cat("\n=== Automatisk opsummering ===\n")

cat(sprintf(
  "Højeste Sharpe-ratio: %.3f ved korrelationsvindue %d uger og volatilitetsvindue %d uger.\n",
  best_spec$sharpe,
  best_spec$corr_window_weeks,
  best_spec$vol_window_weeks
))

cat(sprintf(
  "Laveste Sharpe-ratio: %.3f ved korrelationsvindue %d uger og volatilitetsvindue %d uger.\n",
  worst_spec$sharpe,
  worst_spec$corr_window_weeks,
  worst_spec$vol_window_weeks
))

cat(sprintf(
  "Baseline-specifikationen Corr 260 / Vol 52 har Sharpe-ratio: %.3f.\n",
  robust_corr_vol_grid[
    corr_window_weeks == 260 & vol_window_weeks == 52,
    sharpe
  ]
))

cat("\n=== STEP 5.1 FÆRDIG — korrelations- og volatilitetsvinduer testet ===\n")

# ==============================================================================
# STEP 5.2: ROBUSTHEDSTEST — SHRINKAGE = 0.9
# ==============================================================================
# Formål:
#   Test om BAB-resultaterne er robuste over for ændring i beta-shrinkage.
#
# Baseline:
#   beta = 0.6 * beta_raw + 0.4
#
# Robusthedstest:
#   beta = 0.9 * beta_raw + 0.1
#
# Fortolkning:
#   En højere shrinkage-vægt på beta_raw betyder, at den estimerede beta
#   i mindre grad trækkes mod 1. Testen undersøger derfor, om BAB-resultatet
#   afhænger af, hvor kraftigt beta-estimaterne glattes.
# ==============================================================================

cat("\n=== STEP 5.2: Robusthedstest — shrinkage = 0.9 ===\n")

library(data.table)
library(sandwich)
library(lmtest)
library(ggplot2)

# ------------------------------------------------------------------------------
# Hjælpefunktion: konstruér BAB med alternativ shrinkage
# ------------------------------------------------------------------------------

run_bab_shrinkage_test <- function(data_weekly_input,
                                   shrinkage_weight = 0.9,
                                   label = "Shrinkage 0.9") {
  
  cat("\n------------------------------------------------------------\n")
  cat("Robusthedstest:", label, "\n")
  cat("Shrinkage-vægt på beta_raw:", shrinkage_weight, "\n")
  cat("Shrinkage mod 1:", 1 - shrinkage_weight, "\n")
  cat("------------------------------------------------------------\n")
  
  dt <- copy(data_weekly_input)
  setorder(dt, id, iso_year, iso_week)
  
  # --------------------------------------------------------------------------
  # 1) Brug eksisterende beta_raw fra hovedanalysen, hvis den findes
  # --------------------------------------------------------------------------
  # data_weekly fra STEP 2 indeholder allerede:
  #   beta_raw = rho * (sigma_i / sigma_m)
  #
  # Vi ændrer derfor kun shrinkage-formlen og holder alt andet konstant.
  # --------------------------------------------------------------------------
  
  if (!"beta_raw" %in% names(dt)) {
    stop("data_weekly_input skal indeholde beta_raw. Brug data_weekly fra STEP 2.")
  }
  
  dt[, beta_shrink_alt := shrinkage_weight * beta_raw + (1 - shrinkage_weight) * 1]
  
  # Lag én uge for at undgå look-ahead bias
  dt[, beta_lag_alt := shift(beta_shrink_alt, n = 1, type = "lag"), by = id]
  
  # --------------------------------------------------------------------------
  # 2) Konstruér BAB-portefølje med alternativ shrinkage-beta
  # --------------------------------------------------------------------------
  
  data_bab_alt <- dt[
    !is.na(beta_lag_alt) &
      !is.na(ret_exc_wk) &
      !is.na(me) &
      is.finite(beta_lag_alt)
  ]
  
  cat(sprintf("Gyldige aktie-uge-observationer: %d\n", nrow(data_bab_alt)))
  cat(sprintf("Antal uger: %d\n", uniqueN(data_bab_alt[, paste(iso_year, iso_week)])))
  
  # Rangér aktier efter alternativ lagget beta
  data_bab_alt[, rnk := rank(beta_lag_alt, ties.method = "average"),
               by = .(iso_year, iso_week)]
  
  data_bab_alt[, z := rnk - mean(rnk),
               by = .(iso_year, iso_week)]
  
  data_bab_alt[, beta_med := median(beta_lag_alt, na.rm = TRUE),
               by = .(iso_year, iso_week)]
  
  data_bab_alt[, portfolio := fifelse(beta_lag_alt <= beta_med, "low", "high")]
  
  # Rank-vægte som i hovedanalysen
  data_bab_alt[portfolio == "low",
               w := -z / sum(-z, na.rm = TRUE),
               by = .(iso_year, iso_week)]
  
  data_bab_alt[portfolio == "high",
               w :=  z / sum( z, na.rm = TRUE),
               by = .(iso_year, iso_week)]
  
  # Lav-beta ben
  low_alt <- data_bab_alt[portfolio == "low", .(
    r_low    = sum(w * ret_exc_wk,   na.rm = TRUE),
    beta_low = sum(w * beta_lag_alt, na.rm = TRUE),
    n_low    = .N
  ), by = .(iso_year, iso_week)]
  
  # Høj-beta ben
  high_alt <- data_bab_alt[portfolio == "high", .(
    r_high    = sum(w * ret_exc_wk,   na.rm = TRUE),
    beta_high = sum(w * beta_lag_alt, na.rm = TRUE),
    n_high    = .N
  ), by = .(iso_year, iso_week)]
  
  # Markedsafkast
  mkt_alt <- data_bab_alt[, .(
    r_mkt = first(mkt_exc_wk)
  ), by = .(iso_year, iso_week)]
  
  # Sammensæt BAB
  bab_alt <- merge(low_alt, high_alt, by = c("iso_year", "iso_week"))
  bab_alt <- merge(bab_alt, mkt_alt,  by = c("iso_year", "iso_week"))
  
  # Tilføj dato
  dates_alt <- data_bab_alt[, .(
    date = max(date)
  ), by = .(iso_year, iso_week)]
  
  bab_alt <- dates_alt[bab_alt, on = .(iso_year, iso_week)]
  setorder(bab_alt, iso_year, iso_week)
  
  # BAB-afkast
  bab_alt[, r_bab := (1 / beta_low) * r_low - (1 / beta_high) * r_high]
  
  # Fjern ugyldige observationer
  bab_alt <- bab_alt[
    !is.na(r_bab) &
      is.finite(r_bab) &
      !is.na(r_mkt) &
      is.finite(r_mkt)
  ]
  
  # --------------------------------------------------------------------------
  # 3) Performance-statistik
  # --------------------------------------------------------------------------
  
  mod_int <- lm(r_bab ~ 1, data = bab_alt)
  nw_int  <- coeftest(
    mod_int,
    vcov = NeweyWest(mod_int, lag = 4, prewhite = FALSE)
  )
  
  mod_capm <- lm(r_bab ~ r_mkt, data = bab_alt)
  nw_capm  <- coeftest(
    mod_capm,
    vcov = NeweyWest(mod_capm, lag = 4, prewhite = FALSE)
  )
  
  mu    <- mean(bab_alt$r_bab, na.rm = TRUE)
  sigma <- sd(bab_alt$r_bab, na.rm = TRUE)
  
  ann_ret    <- mu * 52
  ann_vol    <- sigma * sqrt(52)
  ann_sharpe <- ann_ret / ann_vol
  
  cum_ret  <- cumprod(1 + bab_alt$r_bab)
  peak     <- cummax(cum_ret)
  drawdown <- (cum_ret - peak) / peak
  max_dd   <- min(drawdown, na.rm = TRUE)
  
  out <- data.table(
    test              = label,
    shrinkage_weight  = shrinkage_weight,
    shrinkage_to_one  = 1 - shrinkage_weight,
    n_weeks           = nrow(bab_alt),
    start_date        = min(bab_alt$date),
    end_date          = max(bab_alt$date),
    
    mean_weekly_pct   = mu * 100,
    nw_t_mean         = nw_int[1, 3],
    nw_p_mean         = nw_int[1, 4],
    
    ann_return_pct    = ann_ret * 100,
    ann_vol_pct       = ann_vol * 100,
    sharpe            = ann_sharpe,
    max_drawdown_pct  = max_dd * 100,
    
    capm_alpha_w_pct  = nw_capm[1, 1] * 100,
    capm_alpha_t      = nw_capm[1, 3],
    capm_alpha_p      = nw_capm[1, 4],
    capm_beta         = nw_capm[2, 1],
    capm_r2           = summary(mod_capm)$r.squared,
    
    beta_low_avg      = mean(bab_alt$beta_low,  na.rm = TRUE),
    beta_high_avg     = mean(bab_alt$beta_high, na.rm = TRUE),
    n_low_avg         = mean(bab_alt$n_low,     na.rm = TRUE),
    n_high_avg        = mean(bab_alt$n_high,    na.rm = TRUE)
  )
  
  cat(sprintf("Ann. afkast      : %+.2f%%\n", out$ann_return_pct))
  cat(sprintf("Sharpe ratio     : %.3f\n",   out$sharpe))
  cat(sprintf("CAPM alfa, uge   : %+.4f%%, t = %.3f\n",
              out$capm_alpha_w_pct, out$capm_alpha_t))
  cat(sprintf("CAPM beta        : %.3f\n", out$capm_beta))
  cat(sprintf("Max drawdown     : %+.2f%%\n", out$max_drawdown_pct))
  
  return(list(
    summary = out,
    bab     = bab_alt,
    data    = data_bab_alt
  ))
}

# ------------------------------------------------------------------------------
# Kør shrinkage-test
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Kør shrinkage-test for 0.9 og 0.3
# ------------------------------------------------------------------------------

robust_shrink_09 <- run_bab_shrinkage_test(
  data_weekly_input = data_weekly,
  shrinkage_weight  = 0.9,
  label             = "Shrinkage 0.9"
)

robust_shrink_03 <- run_bab_shrinkage_test(
  data_weekly_input = data_weekly,
  shrinkage_weight  = 0.3,
  label             = "Shrinkage 0.3"
)

# ------------------------------------------------------------------------------
# Sammenlign med baseline
# ------------------------------------------------------------------------------

cat("\n=== STEP 5.2: Sammenligning med baseline ===\n")

baseline_mod_int <- lm(r_bab ~ 1, data = bab_clean)
baseline_nw_int  <- coeftest(
  baseline_mod_int,
  vcov = NeweyWest(baseline_mod_int, lag = 4, prewhite = FALSE)
)

baseline_mod_capm <- lm(r_bab ~ r_mkt, data = bab_clean)
baseline_nw_capm  <- coeftest(
  baseline_mod_capm,
  vcov = NeweyWest(baseline_mod_capm, lag = 4, prewhite = FALSE)
)

baseline_mu    <- mean(bab_clean$r_bab, na.rm = TRUE)
baseline_sigma <- sd(bab_clean$r_bab, na.rm = TRUE)

baseline_cum <- cumprod(1 + bab_clean$r_bab)
baseline_dd  <- (baseline_cum - cummax(baseline_cum)) / cummax(baseline_cum)

baseline_shrink_summary <- data.table(
  test              = "Baseline shrinkage 0.6",
  shrinkage_weight  = 0.6,
  shrinkage_to_one  = 0.4,
  n_weeks           = nrow(bab_clean),
  start_date        = min(bab_clean$date),
  end_date          = max(bab_clean$date),
  
  mean_weekly_pct   = baseline_mu * 100,
  nw_t_mean         = baseline_nw_int[1, 3],
  nw_p_mean         = baseline_nw_int[1, 4],
  
  ann_return_pct    = baseline_mu * 52 * 100,
  ann_vol_pct       = baseline_sigma * sqrt(52) * 100,
  sharpe            = (baseline_mu / baseline_sigma) * sqrt(52),
  max_drawdown_pct  = min(baseline_dd, na.rm = TRUE) * 100,
  
  capm_alpha_w_pct  = baseline_nw_capm[1, 1] * 100,
  capm_alpha_t      = baseline_nw_capm[1, 3],
  capm_alpha_p      = baseline_nw_capm[1, 4],
  capm_beta         = baseline_nw_capm[2, 1],
  capm_r2           = summary(baseline_mod_capm)$r.squared,
  
  beta_low_avg      = mean(bab_clean$beta_low,  na.rm = TRUE),
  beta_high_avg     = mean(bab_clean$beta_high, na.rm = TRUE),
  n_low_avg         = mean(bab_clean$n_low,     na.rm = TRUE),
  n_high_avg        = mean(bab_clean$n_high,    na.rm = TRUE)
)

robust_shrink_summary <- rbindlist(list(
  baseline_shrink_summary,
  robust_shrink_09$summary,
  robust_shrink_03$summary
), fill = TRUE)

print(robust_shrink_summary)

# Gem tabel
fwrite(robust_shrink_summary, "robusthed_shrinkage_03_09.csv")
cat("\nGemt: robusthed_shrinkage_03_09.csv\n")

# Gem RDS-filer
saveRDS(robust_shrink_summary, "robusthed_shrinkage_03_09_summary.rds")
saveRDS(robust_shrink_09$bab,  "bab_robust_shrinkage_09.rds")
saveRDS(robust_shrink_03$bab,  "bab_robust_shrinkage_03.rds")

cat("Gemt:\n")
cat("  -> robusthed_shrinkage_03_09_summary.rds\n")
cat("  -> bab_robust_shrinkage_09.rds\n")
cat("  -> bab_robust_shrinkage_03.rds\n")

# ------------------------------------------------------------------------------
# Plot: kumuleret afkast baseline vs shrinkage 0.9 og 0.3
# ------------------------------------------------------------------------------

bab_base_plot <- copy(bab_clean[, .(date, r_bab)])
bab_base_plot[, test := "Baseline shrinkage 0.6"]

bab_shrink_09_plot <- copy(robust_shrink_09$bab[, .(date, r_bab)])
bab_shrink_09_plot[, test := "Shrinkage 0.9"]

bab_shrink_03_plot <- copy(robust_shrink_03$bab[, .(date, r_bab)])
bab_shrink_03_plot[, test := "Shrinkage 0.3"]

shrink_plot_data <- rbindlist(list(
  bab_base_plot,
  bab_shrink_09_plot,
  bab_shrink_03_plot
), fill = TRUE)

setorder(shrink_plot_data, test, date)

shrink_plot_data[, cum_return := cumprod(1 + r_bab) - 1, by = test]

p_shrink <- ggplot(
  shrink_plot_data,
  aes(x = date, y = cum_return * 100, color = test, linetype = test)
) +
  geom_line(linewidth = 0.75) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title    = "Robusthedstest: BAB-afkast ved ændret shrinkage",
    subtitle = "Sammenligning af baseline shrinkage 0,6 samt alternative shrinkage-vægte 0,9 og 0,3",
    x        = NULL,
    y        = "Kumuleret afkast (%)",
    color    = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey30")
  )

ggsave(
  filename = "robusthed_shrinkage_03_09_kumuleret_afkast.png",
  plot     = p_shrink,
  width    = 12,
  height   = 6,
  dpi      = 300
)

cat("Gemt: robusthed_shrinkage_03_09_kumuleret_afkast.png\n")

# ------------------------------------------------------------------------------
# Plot: Sharpe-ratio baseline vs shrinkage 0.9 og 0.3
# ------------------------------------------------------------------------------

p_shrink_bar <- ggplot(
  robust_shrink_summary,
  aes(x = test, y = sharpe)
) +
  geom_col(width = 0.60, fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title    = "Robusthedstest: Sharpe-ratio ved ændret shrinkage",
    subtitle = "Annualiseret Sharpe-ratio for BAB-strategien",
    x        = NULL,
    y        = "Sharpe-ratio"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(
  filename = "robusthed_shrinkage_03_09_sharpe.png",
  plot     = p_shrink_bar,
  width    = 8,
  height   = 5,
  dpi      = 300
)

cat("Gemt: robusthed_shrinkage_03_09_sharpe.png\n")

cat("\n=== STEP 5.2 FÆRDIG — robusthedstest for shrinkage = 0.9 og 0.3 er gennemført ===\n")

# ==============================================================================
# STEP 5.3: ROBUSTHEDSTEST — ANTAL PORTEFØLJER I PORTEFØLJEKONSTRUKTIONEN
# ==============================================================================
# Formål:
#   Test om BAB-resultaterne er robuste over for antallet af beta-sorterede
#   porteføljer, der anvendes i konstruktionen.
#
# Baseline:
#   - Low/high-porteføljer dannes omkring medianen med rank-vægte.
#
# Robusthedstest:
#   - Kvintiler: long P1 og short P5
#   - Deciler:   long P1 og short P10
#
# Fortolkning:
#   Testen undersøger, om BAB-effekten primært drives af den brede low/high-
#   konstruktion, eller om effekten også findes, når strategien kun bruger
#   yderporteføljerne i beta-fordelingen.
# ==============================================================================

cat("\n=== STEP 5.3: Robusthedstest — kvintiler vs. deciler ===\n")

library(data.table)
library(sandwich)
library(lmtest)
library(ggplot2)

# ------------------------------------------------------------------------------
# Hjælpefunktion: BAB baseret på yderporteføljer
# ------------------------------------------------------------------------------

run_bab_nport_test <- function(data_bab_input,
                               n_portfolios,
                               label,
                               min_stocks_leg = 5) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Robusthedstest:", label, "\n")
  cat("Antal beta-sorterede porteføljer:", n_portfolios, "\n")
  cat("Long-ben: P1 | Short-ben: P", n_portfolios, "\n", sep = "")
  cat("------------------------------------------------------------\n")
  
  dt <- copy(data_bab_input)
  
  if (!all(c("id", "date", "iso_year", "iso_week",
             "beta_lag", "ret_exc_wk", "mkt_exc_wk") %in% names(dt))) {
    stop("data_bab_input skal indeholde id, date, iso_year, iso_week, beta_lag, ret_exc_wk og mkt_exc_wk.")
  }
  
  dt <- dt[
    !is.na(beta_lag) &
      !is.na(ret_exc_wk) &
      !is.na(mkt_exc_wk) &
      is.finite(beta_lag) &
      is.finite(ret_exc_wk) &
      is.finite(mkt_exc_wk)
  ]
  
  setorder(dt, iso_year, iso_week, beta_lag)
  
  # --------------------------------------------------------------------------
  # 1) Tildel aktier til beta-sorterede porteføljer hver uge
  # --------------------------------------------------------------------------
  # Metoden nedenfor svarer til en ntile-opdeling:
  #   portfolio_num = 1          -> laveste beta
  #   portfolio_num = n_portfolios -> højeste beta
  # --------------------------------------------------------------------------
  
  dt[, n_week := .N, by = .(iso_year, iso_week)]
  
  dt[, beta_rank := frank(beta_lag, ties.method = "first"),
     by = .(iso_year, iso_week)]
  
  dt[, portfolio_num := pmin(
    n_portfolios,
    floor((beta_rank - 1) * n_portfolios / n_week) + 1
  )]
  
  # Behold kun yderporteføljerne
  dt_extreme <- dt[portfolio_num %in% c(1, n_portfolios)]
  
  dt_extreme[, side := fifelse(
    portfolio_num == 1,
    "low",
    "high"
  )]
  
  # --------------------------------------------------------------------------
  # 2) Beregn porteføljeafkast og porteføljebeta
  # --------------------------------------------------------------------------
  # Her anvendes equal-weighting inden for hver yderportefølje.
  # Det gør testen enkel og direkte sammenlignelig med klassiske
  # sorteringsporteføljer.
  # --------------------------------------------------------------------------
  
  low_port <- dt_extreme[side == "low", .(
    r_low    = mean(ret_exc_wk, na.rm = TRUE),
    beta_low = mean(beta_lag,   na.rm = TRUE),
    n_low    = .N
  ), by = .(iso_year, iso_week)]
  
  high_port <- dt_extreme[side == "high", .(
    r_high    = mean(ret_exc_wk, na.rm = TRUE),
    beta_high = mean(beta_lag,   na.rm = TRUE),
    n_high    = .N
  ), by = .(iso_year, iso_week)]
  
  mkt_weekly <- dt[, .(
    r_mkt = first(mkt_exc_wk),
    date  = max(date)
  ), by = .(iso_year, iso_week)]
  
  bab_alt <- merge(low_port, high_port, by = c("iso_year", "iso_week"))
  bab_alt <- merge(bab_alt, mkt_weekly, by = c("iso_year", "iso_week"))
  
  setorder(bab_alt, iso_year, iso_week)
  
  # Krav om minimum antal aktier i begge ben
  bab_alt <- bab_alt[n_low >= min_stocks_leg & n_high >= min_stocks_leg]
  
  # BAB-afkast: long lav-beta skaleret til beta = 1,
  # short høj-beta skaleret til beta = 1
  bab_alt[, r_bab := (1 / beta_low) * r_low - (1 / beta_high) * r_high]
  
  bab_alt <- bab_alt[
    !is.na(r_bab) &
      is.finite(r_bab) &
      !is.na(r_mkt) &
      is.finite(r_mkt)
  ]
  
  # --------------------------------------------------------------------------
  # 3) Performance-statistik
  # --------------------------------------------------------------------------
  
  mod_int <- lm(r_bab ~ 1, data = bab_alt)
  nw_int  <- coeftest(
    mod_int,
    vcov = NeweyWest(mod_int, lag = 4, prewhite = FALSE)
  )
  
  mod_capm <- lm(r_bab ~ r_mkt, data = bab_alt)
  nw_capm  <- coeftest(
    mod_capm,
    vcov = NeweyWest(mod_capm, lag = 4, prewhite = FALSE)
  )
  
  mu    <- mean(bab_alt$r_bab, na.rm = TRUE)
  sigma <- sd(bab_alt$r_bab, na.rm = TRUE)
  
  ann_ret    <- mu * 52
  ann_vol    <- sigma * sqrt(52)
  ann_sharpe <- ann_ret / ann_vol
  
  cum_ret  <- cumprod(1 + bab_alt$r_bab)
  peak     <- cummax(cum_ret)
  drawdown <- (cum_ret - peak) / peak
  max_dd   <- min(drawdown, na.rm = TRUE)
  
  out <- data.table(
    test              = label,
    n_portfolios      = n_portfolios,
    long_portfolio    = "P1",
    short_portfolio   = paste0("P", n_portfolios),
    weighting         = "Equal-weighted yderporteføljer",
    n_weeks           = nrow(bab_alt),
    start_date        = min(bab_alt$date),
    end_date          = max(bab_alt$date),
    
    mean_weekly_pct   = mu * 100,
    nw_t_mean         = nw_int[1, 3],
    nw_p_mean         = nw_int[1, 4],
    
    ann_return_pct    = ann_ret * 100,
    ann_vol_pct       = ann_vol * 100,
    sharpe            = ann_sharpe,
    max_drawdown_pct  = max_dd * 100,
    
    capm_alpha_w_pct  = nw_capm[1, 1] * 100,
    capm_alpha_t      = nw_capm[1, 3],
    capm_alpha_p      = nw_capm[1, 4],
    capm_beta         = nw_capm[2, 1],
    capm_r2           = summary(mod_capm)$r.squared,
    
    beta_low_avg      = mean(bab_alt$beta_low,  na.rm = TRUE),
    beta_high_avg     = mean(bab_alt$beta_high, na.rm = TRUE),
    n_low_avg         = mean(bab_alt$n_low,     na.rm = TRUE),
    n_high_avg        = mean(bab_alt$n_high,    na.rm = TRUE)
  )
  
  cat(sprintf("Ann. afkast      : %+.2f%%\n", out$ann_return_pct))
  cat(sprintf("Sharpe ratio     : %.3f\n",   out$sharpe))
  cat(sprintf("CAPM alfa, uge   : %+.4f%%, t = %.3f\n",
              out$capm_alpha_w_pct, out$capm_alpha_t))
  cat(sprintf("CAPM beta        : %.3f\n", out$capm_beta))
  cat(sprintf("Max drawdown     : %+.2f%%\n", out$max_drawdown_pct))
  cat(sprintf("Gns. antal aktier i low-ben : %.0f\n", out$n_low_avg))
  cat(sprintf("Gns. antal aktier i high-ben: %.0f\n", out$n_high_avg))
  
  return(list(
    summary = out,
    bab     = bab_alt,
    data    = dt_extreme
  ))
}

# ------------------------------------------------------------------------------
# Kør robusthedstest: kvintiler og deciler
# ------------------------------------------------------------------------------

robust_quintile <- run_bab_nport_test(
  data_bab_input = data_bab,
  n_portfolios   = 5,
  label          = "Kvintiler: P1 vs. P5"
)

robust_decile <- run_bab_nport_test(
  data_bab_input = data_bab,
  n_portfolios   = 10,
  label          = "Deciler: P1 vs. P10"
)

# ------------------------------------------------------------------------------
# Baseline-statistik til sammenligning
# ------------------------------------------------------------------------------

cat("\n=== STEP 5.3: Sammenligning med baseline ===\n")

baseline_mod_int <- lm(r_bab ~ 1, data = bab_clean)
baseline_nw_int  <- coeftest(
  baseline_mod_int,
  vcov = NeweyWest(baseline_mod_int, lag = 4, prewhite = FALSE)
)

baseline_mod_capm <- lm(r_bab ~ r_mkt, data = bab_clean)
baseline_nw_capm  <- coeftest(
  baseline_mod_capm,
  vcov = NeweyWest(baseline_mod_capm, lag = 4, prewhite = FALSE)
)

baseline_mu    <- mean(bab_clean$r_bab, na.rm = TRUE)
baseline_sigma <- sd(bab_clean$r_bab, na.rm = TRUE)

baseline_cum <- cumprod(1 + bab_clean$r_bab)
baseline_dd  <- (baseline_cum - cummax(baseline_cum)) / cummax(baseline_cum)

baseline_nport_summary <- data.table(
  test              = "Baseline median-split",
  n_portfolios      = NA_integer_,
  long_portfolio    = "Lav beta",
  short_portfolio   = "Høj beta",
  weighting         = "Rank-vægtet median-split",
  n_weeks           = nrow(bab_clean),
  start_date        = min(bab_clean$date),
  end_date          = max(bab_clean$date),
  
  mean_weekly_pct   = baseline_mu * 100,
  nw_t_mean         = baseline_nw_int[1, 3],
  nw_p_mean         = baseline_nw_int[1, 4],
  
  ann_return_pct    = baseline_mu * 52 * 100,
  ann_vol_pct       = baseline_sigma * sqrt(52) * 100,
  sharpe            = (baseline_mu / baseline_sigma) * sqrt(52),
  max_drawdown_pct  = min(baseline_dd, na.rm = TRUE) * 100,
  
  capm_alpha_w_pct  = baseline_nw_capm[1, 1] * 100,
  capm_alpha_t      = baseline_nw_capm[1, 3],
  capm_alpha_p      = baseline_nw_capm[1, 4],
  capm_beta         = baseline_nw_capm[2, 1],
  capm_r2           = summary(baseline_mod_capm)$r.squared,
  
  beta_low_avg      = mean(bab_clean$beta_low,  na.rm = TRUE),
  beta_high_avg     = mean(bab_clean$beta_high, na.rm = TRUE),
  n_low_avg         = mean(bab_clean$n_low,     na.rm = TRUE),
  n_high_avg        = mean(bab_clean$n_high,    na.rm = TRUE)
)

robust_nport_summary <- rbindlist(list(
  baseline_nport_summary,
  robust_quintile$summary,
  robust_decile$summary
), fill = TRUE)

print(robust_nport_summary)

# Gem tabel
fwrite(robust_nport_summary, "robusthed_portefoljer_kvintiler_deciler.csv")
cat("\nGemt: robusthed_portefoljer_kvintiler_deciler.csv\n")

# Gem RDS-filer
saveRDS(robust_nport_summary, "robusthed_portefoljer_kvintiler_deciler_summary.rds")
saveRDS(robust_quintile$bab,  "bab_robust_kvintiler.rds")
saveRDS(robust_decile$bab,    "bab_robust_deciler.rds")

cat("Gemt:\n")
cat("  -> robusthed_portefoljer_kvintiler_deciler_summary.rds\n")
cat("  -> bab_robust_kvintiler.rds\n")
cat("  -> bab_robust_deciler.rds\n")

# ------------------------------------------------------------------------------
# Plot: kumuleret afkast
# ------------------------------------------------------------------------------

bab_base_plot <- copy(bab_clean[, .(date, r_bab)])
bab_base_plot[, test := "Baseline median-split"]

bab_quintile_plot <- copy(robust_quintile$bab[, .(date, r_bab)])
bab_quintile_plot[, test := "Kvintiler: P1 vs. P5"]

bab_decile_plot <- copy(robust_decile$bab[, .(date, r_bab)])
bab_decile_plot[, test := "Deciler: P1 vs. P10"]

nport_plot_data <- rbindlist(list(
  bab_base_plot,
  bab_quintile_plot,
  bab_decile_plot
), fill = TRUE)

setorder(nport_plot_data, test, date)

nport_plot_data[, cum_return := cumprod(1 + r_bab) - 1, by = test]

p_nport <- ggplot(
  nport_plot_data,
  aes(x = date, y = cum_return * 100, color = test, linetype = test)
) +
  geom_line(linewidth = 0.75) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title    = "Robusthedstest: BAB-afkast ved alternative porteføljeinddelinger",
    subtitle = "Baseline median-split sammenlignet med kvintiler og deciler",
    x        = NULL,
    y        = "Kumuleret afkast (%)",
    color    = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey30")
  )

ggsave(
  filename = "robusthed_portefoljer_kvintiler_deciler_kumuleret_afkast.png",
  plot     = p_nport,
  width    = 12,
  height   = 6,
  dpi      = 300
)

cat("Gemt: robusthed_portefoljer_kvintiler_deciler_kumuleret_afkast.png\n")

# ------------------------------------------------------------------------------
# Plot: Sharpe-ratio
# ------------------------------------------------------------------------------

p_nport_bar <- ggplot(
  robust_nport_summary,
  aes(x = test, y = sharpe)
) +
  geom_col(width = 0.60, fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title    = "Robusthedstest: Sharpe-ratio ved alternative porteføljeinddelinger",
    subtitle = "Annualiseret Sharpe-ratio for BAB-strategien",
    x        = NULL,
    y        = "Sharpe-ratio"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    axis.text.x   = element_text(angle = 15, hjust = 1)
  )

ggsave(
  filename = "robusthed_portefoljer_kvintiler_deciler_sharpe.png",
  plot     = p_nport_bar,
  width    = 9,
  height   = 5,
  dpi      = 300
)

cat("Gemt: robusthed_portefoljer_kvintiler_deciler_sharpe.png\n")

cat("\n=== STEP 5.3 FÆRDIG — robusthedstest for kvintiler og deciler er gennemført ===\n")

# ==============================================================================
# STEP 5.4: ROBUSTHEDSTEST — VÆGTNINGSMETODE
# ==============================================================================
# Formål:
#   Test om BAB-resultaterne er robuste over for vægtningen inden for
#   low- og high-beta-porteføljerne.
#
# Baseline:
#   - Rank-vægtet median-split
#
# Robusthedstest:
#   - Equal-weighted median-split
#   - Value-weighted median-split
#
# Fortolkning:
#   Equal-weighted tester, om resultaterne drives af den generelle beta-effekt
#   på tværs af aktier.
#
#   Value-weighted tester, om effekten fortsat findes, når større selskaber
#   får større vægt i porteføljerne.
# ==============================================================================

cat("\n=== STEP 5.4: Robusthedstest — equal-weighted vs. value-weighted ===\n")

library(data.table)
library(sandwich)
library(lmtest)
library(ggplot2)

# ------------------------------------------------------------------------------
# Hjælpefunktion: BAB med alternativ vægtning
# ------------------------------------------------------------------------------

run_bab_weighting_test <- function(data_bab_input,
                                   weighting_method = c("equal", "value"),
                                   label,
                                   min_stocks_leg = 5) {
  
  weighting_method <- match.arg(weighting_method)
  
  cat("\n------------------------------------------------------------\n")
  cat("Robusthedstest:", label, "\n")
  cat("Vægtning:", weighting_method, "\n")
  cat("------------------------------------------------------------\n")
  
  dt <- copy(data_bab_input)
  
  required_cols <- c(
    "id", "date", "iso_year", "iso_week",
    "beta_lag", "ret_exc_wk", "mkt_exc_wk", "me"
  )
  
  if (!all(required_cols %in% names(dt))) {
    stop("data_bab_input mangler en eller flere nødvendige kolonner: ",
         paste(setdiff(required_cols, names(dt)), collapse = ", "))
  }
  
  dt <- dt[
    !is.na(beta_lag) &
      !is.na(ret_exc_wk) &
      !is.na(mkt_exc_wk) &
      !is.na(me) &
      me > 0 &
      is.finite(beta_lag) &
      is.finite(ret_exc_wk) &
      is.finite(mkt_exc_wk) &
      is.finite(me)
  ]
  
  setorder(dt, iso_year, iso_week, beta_lag)
  
  # --------------------------------------------------------------------------
  # 1) Median-split som i baseline
  # --------------------------------------------------------------------------
  
  dt[, beta_med := median(beta_lag, na.rm = TRUE),
     by = .(iso_year, iso_week)]
  
  dt[, portfolio := fifelse(beta_lag <= beta_med, "low", "high")]
  
  # --------------------------------------------------------------------------
  # 2) Beregn vægte inden for hvert ben
  # --------------------------------------------------------------------------
  
  if (weighting_method == "equal") {
    
    # Equal-weighted:
    # Alle aktier inden for hvert ben får samme vægt.
    dt[, w_alt := 1 / .N, by = .(iso_year, iso_week, portfolio)]
    
  } else if (weighting_method == "value") {
    
    # Value-weighted:
    # Aktier vægtes efter markedsværdi inden for hvert ben.
    dt[, w_alt := me / sum(me, na.rm = TRUE),
       by = .(iso_year, iso_week, portfolio)]
  }
  
  # Sikkerhedstjek: vægte summerer til 1 i hvert ben
  weight_check <- dt[, .(
    w_sum = sum(w_alt, na.rm = TRUE)
  ), by = .(iso_year, iso_week, portfolio)]
  
  if (any(abs(weight_check$w_sum - 1) > 1e-6, na.rm = TRUE)) {
    warning("Nogle porteføljevægte summerer ikke præcist til 1.")
  }
  
  # --------------------------------------------------------------------------
  # 3) Beregn low- og high-beta-ben
  # --------------------------------------------------------------------------
  
  low_port <- dt[portfolio == "low", .(
    r_low    = sum(w_alt * ret_exc_wk, na.rm = TRUE),
    beta_low = sum(w_alt * beta_lag,   na.rm = TRUE),
    n_low    = .N
  ), by = .(iso_year, iso_week)]
  
  high_port <- dt[portfolio == "high", .(
    r_high    = sum(w_alt * ret_exc_wk, na.rm = TRUE),
    beta_high = sum(w_alt * beta_lag,   na.rm = TRUE),
    n_high    = .N
  ), by = .(iso_year, iso_week)]
  
  mkt_weekly <- dt[, .(
    r_mkt = first(mkt_exc_wk),
    date  = max(date)
  ), by = .(iso_year, iso_week)]
  
  # --------------------------------------------------------------------------
  # 4) Konstruér BAB
  # --------------------------------------------------------------------------
  
  bab_alt <- merge(low_port, high_port, by = c("iso_year", "iso_week"))
  bab_alt <- merge(bab_alt, mkt_weekly, by = c("iso_year", "iso_week"))
  
  setorder(bab_alt, iso_year, iso_week)
  
  # Minimum antal aktier i begge ben
  bab_alt <- bab_alt[n_low >= min_stocks_leg & n_high >= min_stocks_leg]
  
  # BAB-afkast:
  # long lav-beta skaleret til beta = 1
  # short høj-beta skaleret til beta = 1
  bab_alt[, r_bab := (1 / beta_low) * r_low - (1 / beta_high) * r_high]
  
  bab_alt <- bab_alt[
    !is.na(r_bab) &
      is.finite(r_bab) &
      !is.na(r_mkt) &
      is.finite(r_mkt)
  ]
  
  # --------------------------------------------------------------------------
  # 5) Performance-statistik
  # --------------------------------------------------------------------------
  
  mod_int <- lm(r_bab ~ 1, data = bab_alt)
  nw_int  <- coeftest(
    mod_int,
    vcov = NeweyWest(mod_int, lag = 4, prewhite = FALSE)
  )
  
  mod_capm <- lm(r_bab ~ r_mkt, data = bab_alt)
  nw_capm  <- coeftest(
    mod_capm,
    vcov = NeweyWest(mod_capm, lag = 4, prewhite = FALSE)
  )
  
  mu    <- mean(bab_alt$r_bab, na.rm = TRUE)
  sigma <- sd(bab_alt$r_bab, na.rm = TRUE)
  
  ann_ret    <- mu * 52
  ann_vol    <- sigma * sqrt(52)
  ann_sharpe <- ann_ret / ann_vol
  
  cum_ret  <- cumprod(1 + bab_alt$r_bab)
  peak     <- cummax(cum_ret)
  drawdown <- (cum_ret - peak) / peak
  max_dd   <- min(drawdown, na.rm = TRUE)
  
  out <- data.table(
    test              = label,
    weighting         = weighting_method,
    n_weeks           = nrow(bab_alt),
    start_date        = min(bab_alt$date),
    end_date          = max(bab_alt$date),
    
    mean_weekly_pct   = mu * 100,
    nw_t_mean         = nw_int[1, 3],
    nw_p_mean         = nw_int[1, 4],
    
    ann_return_pct    = ann_ret * 100,
    ann_vol_pct       = ann_vol * 100,
    sharpe            = ann_sharpe,
    max_drawdown_pct  = max_dd * 100,
    
    capm_alpha_w_pct  = nw_capm[1, 1] * 100,
    capm_alpha_t      = nw_capm[1, 3],
    capm_alpha_p      = nw_capm[1, 4],
    capm_beta         = nw_capm[2, 1],
    capm_r2           = summary(mod_capm)$r.squared,
    
    beta_low_avg      = mean(bab_alt$beta_low,  na.rm = TRUE),
    beta_high_avg     = mean(bab_alt$beta_high, na.rm = TRUE),
    n_low_avg         = mean(bab_alt$n_low,     na.rm = TRUE),
    n_high_avg        = mean(bab_alt$n_high,    na.rm = TRUE)
  )
  
  cat(sprintf("Ann. afkast      : %+.2f%%\n", out$ann_return_pct))
  cat(sprintf("Sharpe ratio     : %.3f\n",   out$sharpe))
  cat(sprintf("CAPM alfa, uge   : %+.4f%%, t = %.3f\n",
              out$capm_alpha_w_pct, out$capm_alpha_t))
  cat(sprintf("CAPM beta        : %.3f\n", out$capm_beta))
  cat(sprintf("Max drawdown     : %+.2f%%\n", out$max_drawdown_pct))
  cat(sprintf("Gns. beta low    : %.3f\n", out$beta_low_avg))
  cat(sprintf("Gns. beta high   : %.3f\n", out$beta_high_avg))
  
  return(list(
    summary = out,
    bab     = bab_alt,
    data    = dt
  ))
}

# ------------------------------------------------------------------------------
# Kør robusthedstest: equal-weighted og value-weighted
# ------------------------------------------------------------------------------

robust_equal_weight <- run_bab_weighting_test(
  data_bab_input   = data_bab,
  weighting_method = "equal",
  label            = "Equal-weighted"
)

robust_value_weight <- run_bab_weighting_test(
  data_bab_input   = data_bab,
  weighting_method = "value",
  label            = "Value-weighted"
)

# ------------------------------------------------------------------------------
# Baseline-statistik til sammenligning
# ------------------------------------------------------------------------------

cat("\n=== STEP 5.4: Sammenligning med baseline ===\n")

baseline_mod_int <- lm(r_bab ~ 1, data = bab_clean)
baseline_nw_int  <- coeftest(
  baseline_mod_int,
  vcov = NeweyWest(baseline_mod_int, lag = 4, prewhite = FALSE)
)

baseline_mod_capm <- lm(r_bab ~ r_mkt, data = bab_clean)
baseline_nw_capm  <- coeftest(
  baseline_mod_capm,
  vcov = NeweyWest(baseline_mod_capm, lag = 4, prewhite = FALSE)
)

baseline_mu    <- mean(bab_clean$r_bab, na.rm = TRUE)
baseline_sigma <- sd(bab_clean$r_bab, na.rm = TRUE)

baseline_cum <- cumprod(1 + bab_clean$r_bab)
baseline_dd  <- (baseline_cum - cummax(baseline_cum)) / cummax(baseline_cum)

baseline_weight_summary <- data.table(
  test              = "Baseline rank-weighted",
  weighting         = "rank",
  n_weeks           = nrow(bab_clean),
  start_date        = min(bab_clean$date),
  end_date          = max(bab_clean$date),
  
  mean_weekly_pct   = baseline_mu * 100,
  nw_t_mean         = baseline_nw_int[1, 3],
  nw_p_mean         = baseline_nw_int[1, 4],
  
  ann_return_pct    = baseline_mu * 52 * 100,
  ann_vol_pct       = baseline_sigma * sqrt(52) * 100,
  sharpe            = (baseline_mu / baseline_sigma) * sqrt(52),
  max_drawdown_pct  = min(baseline_dd, na.rm = TRUE) * 100,
  
  capm_alpha_w_pct  = baseline_nw_capm[1, 1] * 100,
  capm_alpha_t      = baseline_nw_capm[1, 3],
  capm_alpha_p      = baseline_nw_capm[1, 4],
  capm_beta         = baseline_nw_capm[2, 1],
  capm_r2           = summary(baseline_mod_capm)$r.squared,
  
  beta_low_avg      = mean(bab_clean$beta_low,  na.rm = TRUE),
  beta_high_avg     = mean(bab_clean$beta_high, na.rm = TRUE),
  n_low_avg         = mean(bab_clean$n_low,     na.rm = TRUE),
  n_high_avg        = mean(bab_clean$n_high,    na.rm = TRUE)
)

robust_weight_summary <- rbindlist(list(
  baseline_weight_summary,
  robust_equal_weight$summary,
  robust_value_weight$summary
), fill = TRUE)

print(robust_weight_summary)

# Gem tabel
fwrite(robust_weight_summary, "robusthed_vaegtning_equal_value.csv")
cat("\nGemt: robusthed_vaegtning_equal_value.csv\n")

# Gem RDS-filer
saveRDS(robust_weight_summary, "robusthed_vaegtning_equal_value_summary.rds")
saveRDS(robust_equal_weight$bab, "bab_robust_equal_weighted.rds")
saveRDS(robust_value_weight$bab, "bab_robust_value_weighted.rds")

cat("Gemt:\n")
cat("  -> robusthed_vaegtning_equal_value_summary.rds\n")
cat("  -> bab_robust_equal_weighted.rds\n")
cat("  -> bab_robust_value_weighted.rds\n")

# ------------------------------------------------------------------------------
# Plot: kumuleret afkast
# ------------------------------------------------------------------------------

bab_base_plot <- copy(bab_clean[, .(date, r_bab)])
bab_base_plot[, test := "Baseline rank-weighted"]

bab_equal_plot <- copy(robust_equal_weight$bab[, .(date, r_bab)])
bab_equal_plot[, test := "Equal-weighted"]

bab_value_plot <- copy(robust_value_weight$bab[, .(date, r_bab)])
bab_value_plot[, test := "Value-weighted"]

weight_plot_data <- rbindlist(list(
  bab_base_plot,
  bab_equal_plot,
  bab_value_plot
), fill = TRUE)

setorder(weight_plot_data, test, date)

weight_plot_data[, cum_return := cumprod(1 + r_bab) - 1, by = test]

p_weight <- ggplot(
  weight_plot_data,
  aes(x = date, y = cum_return * 100, color = test, linetype = test)
) +
  geom_line(linewidth = 0.75) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title    = "Robusthedstest: BAB-afkast ved alternative vægtningsmetoder",
    subtitle = "Baseline rank-weighted sammenlignet med equal-weighted og value-weighted",
    x        = NULL,
    y        = "Kumuleret afkast (%)",
    color    = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey30")
  )

ggsave(
  filename = "robusthed_vaegtning_equal_value_kumuleret_afkast.png",
  plot     = p_weight,
  width    = 12,
  height   = 6,
  dpi      = 300
)

cat("Gemt: robusthed_vaegtning_equal_value_kumuleret_afkast.png\n")

# ------------------------------------------------------------------------------
# Plot: Sharpe-ratio
# ------------------------------------------------------------------------------

p_weight_bar <- ggplot(
  robust_weight_summary,
  aes(x = test, y = sharpe)
) +
  geom_col(width = 0.60, fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title    = "Robusthedstest: Sharpe-ratio ved alternative vægtningsmetoder",
    subtitle = "Annualiseret Sharpe-ratio for BAB-strategien",
    x        = NULL,
    y        = "Sharpe-ratio"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    axis.text.x   = element_text(angle = 15, hjust = 1)
  )

ggsave(
  filename = "robusthed_vaegtning_equal_value_sharpe.png",
  plot     = p_weight_bar,
  width    = 9,
  height   = 5,
  dpi      = 300
)

cat("Gemt: robusthed_vaegtning_equal_value_sharpe.png\n")

cat("\n=== STEP 5.4 FÆRDIG — robusthedstest for vægtning er gennemført ===\n")

# ==============================================================================
# STEP 5.5: ROBUSTHEDSTEST — VALUE-WEIGHTED MED EXPONENTIAL DECAY
# ==============================================================================
# Formål:
#   Test om BAB-resultaterne er robuste over for en porteføljekonstruktion,
#   hvor aktier vægtes efter markedsværdi, men hvor vægtene samtidig justeres
#   for tidsafhængighed ved hjælp af en eksponentiel decay-funktion.
#
# Vægt:
#   w_i ∝ ME_i * lambda^(T - t_i)
#
# hvor:
#   ME_i  = aktiens markedsværdi
#   T     = seneste dato i den pågældende uge
#   t_i   = aktiens observationsdato
#   lambda ∈ (0,1)
#
# Fortolkning:
#   Når lambda < 1, får nyere observationer højere vægt end ældre observationer.
#   Jo lavere lambda er, desto stærkere nedvægtes ældre observationer.
# ==============================================================================

cat("\n=== STEP 5.5: Robusthedstest — value-weighted med exponential decay ===\n")

library(data.table)
library(sandwich)
library(lmtest)
library(ggplot2)

# ------------------------------------------------------------------------------
# Hjælpefunktion: BAB med value-weighting og exponential decay
# ------------------------------------------------------------------------------

run_bab_decay_weight_test <- function(data_bab_input,
                                      lambda = 0.95,
                                      label = "Value-weighted decay",
                                      min_stocks_leg = 5) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Robusthedstest:", label, "\n")
  cat("Lambda:", lambda, "\n")
  cat("------------------------------------------------------------\n")
  
  if (lambda <= 0 || lambda > 1) {
    stop("lambda skal være i intervallet (0, 1].")
  }
  
  dt <- copy(data_bab_input)
  
  required_cols <- c(
    "id", "date", "iso_year", "iso_week",
    "beta_lag", "ret_exc_wk", "mkt_exc_wk", "me"
  )
  
  if (!all(required_cols %in% names(dt))) {
    stop("data_bab_input mangler en eller flere nødvendige kolonner: ",
         paste(setdiff(required_cols, names(dt)), collapse = ", "))
  }
  
  dt <- dt[
    !is.na(beta_lag) &
      !is.na(ret_exc_wk) &
      !is.na(mkt_exc_wk) &
      !is.na(me) &
      !is.na(date) &
      me > 0 &
      is.finite(beta_lag) &
      is.finite(ret_exc_wk) &
      is.finite(mkt_exc_wk) &
      is.finite(me)
  ]
  
  setorder(dt, iso_year, iso_week, beta_lag)
  
  # --------------------------------------------------------------------------
  # 1) Median-split som i baseline
  # --------------------------------------------------------------------------
  
  dt[, beta_med := median(beta_lag, na.rm = TRUE),
     by = .(iso_year, iso_week)]
  
  dt[, portfolio := fifelse(beta_lag <= beta_med, "low", "high")]
  
  # --------------------------------------------------------------------------
  # 2) Beregn decay-faktor
  # --------------------------------------------------------------------------
  # T er seneste observationsdato i den pågældende uge.
  # age_days = T - t_i.
  #
  # Hvis alle aktier i en uge har samme dato, bliver decay-faktoren 1 for alle,
  # og testen svarer i praksis til value-weighting.
  # --------------------------------------------------------------------------
  
  dt[, week_T := max(date, na.rm = TRUE), by = .(iso_year, iso_week)]
  
  dt[, age_days := as.numeric(week_T - date)]
  
  dt[, decay_factor := lambda ^ age_days]
  
  # Kombiner markedsværdi og decay
  dt[, me_decay := me * decay_factor]
  
  # --------------------------------------------------------------------------
  # 3) Beregn porteføljevægte inden for hvert ben
  # --------------------------------------------------------------------------
  
  dt[, w_alt := me_decay / sum(me_decay, na.rm = TRUE),
     by = .(iso_year, iso_week, portfolio)]
  
  weight_check <- dt[, .(
    w_sum = sum(w_alt, na.rm = TRUE)
  ), by = .(iso_year, iso_week, portfolio)]
  
  if (any(abs(weight_check$w_sum - 1) > 1e-6, na.rm = TRUE)) {
    warning("Nogle decay-vægte summerer ikke præcist til 1.")
  }
  
  # --------------------------------------------------------------------------
  # 4) Beregn low- og high-beta-ben
  # --------------------------------------------------------------------------
  
  low_port <- dt[portfolio == "low", .(
    r_low    = sum(w_alt * ret_exc_wk, na.rm = TRUE),
    beta_low = sum(w_alt * beta_lag,   na.rm = TRUE),
    n_low    = .N
  ), by = .(iso_year, iso_week)]
  
  high_port <- dt[portfolio == "high", .(
    r_high    = sum(w_alt * ret_exc_wk, na.rm = TRUE),
    beta_high = sum(w_alt * beta_lag,   na.rm = TRUE),
    n_high    = .N
  ), by = .(iso_year, iso_week)]
  
  mkt_weekly <- dt[, .(
    r_mkt = first(mkt_exc_wk),
    date  = max(date)
  ), by = .(iso_year, iso_week)]
  
  # --------------------------------------------------------------------------
  # 5) Konstruér BAB
  # --------------------------------------------------------------------------
  
  bab_alt <- merge(low_port, high_port, by = c("iso_year", "iso_week"))
  bab_alt <- merge(bab_alt, mkt_weekly, by = c("iso_year", "iso_week"))
  
  setorder(bab_alt, iso_year, iso_week)
  
  bab_alt <- bab_alt[n_low >= min_stocks_leg & n_high >= min_stocks_leg]
  
  bab_alt[, r_bab := (1 / beta_low) * r_low - (1 / beta_high) * r_high]
  
  bab_alt <- bab_alt[
    !is.na(r_bab) &
      is.finite(r_bab) &
      !is.na(r_mkt) &
      is.finite(r_mkt)
  ]
  
  # --------------------------------------------------------------------------
  # 6) Performance-statistik
  # --------------------------------------------------------------------------
  
  mod_int <- lm(r_bab ~ 1, data = bab_alt)
  nw_int  <- coeftest(
    mod_int,
    vcov = NeweyWest(mod_int, lag = 4, prewhite = FALSE)
  )
  
  mod_capm <- lm(r_bab ~ r_mkt, data = bab_alt)
  nw_capm  <- coeftest(
    mod_capm,
    vcov = NeweyWest(mod_capm, lag = 4, prewhite = FALSE)
  )
  
  mu    <- mean(bab_alt$r_bab, na.rm = TRUE)
  sigma <- sd(bab_alt$r_bab, na.rm = TRUE)
  
  ann_ret    <- mu * 52
  ann_vol    <- sigma * sqrt(52)
  ann_sharpe <- ann_ret / ann_vol
  
  cum_ret  <- cumprod(1 + bab_alt$r_bab)
  peak     <- cummax(cum_ret)
  drawdown <- (cum_ret - peak) / peak
  max_dd   <- min(drawdown, na.rm = TRUE)
  
  out <- data.table(
    test              = label,
    weighting         = "value-weighted decay",
    lambda            = lambda,
    n_weeks           = nrow(bab_alt),
    start_date        = min(bab_alt$date),
    end_date          = max(bab_alt$date),
    
    mean_weekly_pct   = mu * 100,
    nw_t_mean         = nw_int[1, 3],
    nw_p_mean         = nw_int[1, 4],
    
    ann_return_pct    = ann_ret * 100,
    ann_vol_pct       = ann_vol * 100,
    sharpe            = ann_sharpe,
    max_drawdown_pct  = max_dd * 100,
    
    capm_alpha_w_pct  = nw_capm[1, 1] * 100,
    capm_alpha_t      = nw_capm[1, 3],
    capm_alpha_p      = nw_capm[1, 4],
    capm_beta         = nw_capm[2, 1],
    capm_r2           = summary(mod_capm)$r.squared,
    
    beta_low_avg      = mean(bab_alt$beta_low,  na.rm = TRUE),
    beta_high_avg     = mean(bab_alt$beta_high, na.rm = TRUE),
    n_low_avg         = mean(bab_alt$n_low,     na.rm = TRUE),
    n_high_avg        = mean(bab_alt$n_high,    na.rm = TRUE),
    
    avg_age_days      = mean(dt$age_days, na.rm = TRUE),
    max_age_days      = max(dt$age_days,  na.rm = TRUE)
  )
  
  cat(sprintf("Ann. afkast      : %+.2f%%\n", out$ann_return_pct))
  cat(sprintf("Sharpe ratio     : %.3f\n",   out$sharpe))
  cat(sprintf("CAPM alfa, uge   : %+.4f%%, t = %.3f\n",
              out$capm_alpha_w_pct, out$capm_alpha_t))
  cat(sprintf("CAPM beta        : %.3f\n", out$capm_beta))
  cat(sprintf("Max drawdown     : %+.2f%%\n", out$max_drawdown_pct))
  cat(sprintf("Gns. observationsalder: %.2f dage\n", out$avg_age_days))
  cat(sprintf("Maks. observationsalder: %.0f dage\n", out$max_age_days))
  
  return(list(
    summary = out,
    bab     = bab_alt,
    data    = dt
  ))
}

# ------------------------------------------------------------------------------
# Kør robusthedstest med flere lambda-værdier
# ------------------------------------------------------------------------------
# lambda = 1.00 svarer til almindelig value-weighting uden decay.
# lambda = 0.95 giver moderat decay.
# lambda = 0.90 giver stærkere decay.
# ------------------------------------------------------------------------------

robust_decay_100 <- run_bab_decay_weight_test(
  data_bab_input = data_bab,
  lambda         = 1.00,
  label          = "Value-weighted uden decay"
)

robust_decay_095 <- run_bab_decay_weight_test(
  data_bab_input = data_bab,
  lambda         = 0.8,
  label          = "Value-weighted decay lambda 0.8"
)

robust_decay_090 <- run_bab_decay_weight_test(
  data_bab_input = data_bab,
  lambda         = 0.6,
  label          = "Value-weighted decay lambda 0.60"
)

# ------------------------------------------------------------------------------
# Baseline-statistik til sammenligning
# ------------------------------------------------------------------------------

cat("\n=== STEP 5.5: Sammenligning med baseline ===\n")

baseline_mod_int <- lm(r_bab ~ 1, data = bab_clean)
baseline_nw_int  <- coeftest(
  baseline_mod_int,
  vcov = NeweyWest(baseline_mod_int, lag = 4, prewhite = FALSE)
)

baseline_mod_capm <- lm(r_bab ~ r_mkt, data = bab_clean)
baseline_nw_capm  <- coeftest(
  baseline_mod_capm,
  vcov = NeweyWest(baseline_mod_capm, lag = 4, prewhite = FALSE)
)

baseline_mu    <- mean(bab_clean$r_bab, na.rm = TRUE)
baseline_sigma <- sd(bab_clean$r_bab, na.rm = TRUE)

baseline_cum <- cumprod(1 + bab_clean$r_bab)
baseline_dd  <- (baseline_cum - cummax(baseline_cum)) / cummax(baseline_cum)

baseline_decay_summary <- data.table(
  test              = "Baseline rank-weighted",
  weighting         = "rank",
  lambda            = NA_real_,
  n_weeks           = nrow(bab_clean),
  start_date        = min(bab_clean$date),
  end_date          = max(bab_clean$date),
  
  mean_weekly_pct   = baseline_mu * 100,
  nw_t_mean         = baseline_nw_int[1, 3],
  nw_p_mean         = baseline_nw_int[1, 4],
  
  ann_return_pct    = baseline_mu * 52 * 100,
  ann_vol_pct       = baseline_sigma * sqrt(52) * 100,
  sharpe            = (baseline_mu / baseline_sigma) * sqrt(52),
  max_drawdown_pct  = min(baseline_dd, na.rm = TRUE) * 100,
  
  capm_alpha_w_pct  = baseline_nw_capm[1, 1] * 100,
  capm_alpha_t      = baseline_nw_capm[1, 3],
  capm_alpha_p      = baseline_nw_capm[1, 4],
  capm_beta         = baseline_nw_capm[2, 1],
  capm_r2           = summary(baseline_mod_capm)$r.squared,
  
  beta_low_avg      = mean(bab_clean$beta_low,  na.rm = TRUE),
  beta_high_avg     = mean(bab_clean$beta_high, na.rm = TRUE),
  n_low_avg         = mean(bab_clean$n_low,     na.rm = TRUE),
  n_high_avg        = mean(bab_clean$n_high,    na.rm = TRUE),
  
  avg_age_days      = NA_real_,
  max_age_days      = NA_real_
)

robust_decay_summary <- rbindlist(list(
  baseline_decay_summary,
  robust_decay_100$summary,
  robust_decay_095$summary,
  robust_decay_090$summary
), fill = TRUE)

print(robust_decay_summary)

# Gem tabel
fwrite(robust_decay_summary, "robusthed_value_decay.csv")
cat("\nGemt: robusthed_value_decay.csv\n")

# Gem RDS-filer
saveRDS(robust_decay_summary, "robusthed_value_decay_summary.rds")
saveRDS(robust_decay_100$bab, "bab_robust_value_no_decay.rds")
saveRDS(robust_decay_095$bab, "bab_robust_value_decay_095.rds")
saveRDS(robust_decay_090$bab, "bab_robust_value_decay_090.rds")

cat("Gemt:\n")
cat("  -> robusthed_value_decay_summary.rds\n")
cat("  -> bab_robust_value_no_decay.rds\n")
cat("  -> bab_robust_value_decay_095.rds\n")
cat("  -> bab_robust_value_decay_090.rds\n")

# ------------------------------------------------------------------------------
# Plot: kumuleret afkast
# ------------------------------------------------------------------------------

bab_base_plot <- copy(bab_clean[, .(date, r_bab)])
bab_base_plot[, test := "Baseline rank-weighted"]

bab_decay_100_plot <- copy(robust_decay_100$bab[, .(date, r_bab)])
bab_decay_100_plot[, test := "Value-weighted uden decay"]

bab_decay_095_plot <- copy(robust_decay_095$bab[, .(date, r_bab)])
bab_decay_095_plot[, test := "Decay lambda 0.8"]

bab_decay_090_plot <- copy(robust_decay_090$bab[, .(date, r_bab)])
bab_decay_090_plot[, test := "Decay lambda 0.60"]

decay_plot_data <- rbindlist(list(
  bab_base_plot,
  bab_decay_100_plot,
  bab_decay_095_plot,
  bab_decay_090_plot
), fill = TRUE)

setorder(decay_plot_data, test, date)

decay_plot_data[, cum_return := cumprod(1 + r_bab) - 1, by = test]

p_decay <- ggplot(
  decay_plot_data,
  aes(x = date, y = cum_return * 100, color = test, linetype = test)
) +
  geom_line(linewidth = 0.75) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title    = "Robusthedstest: BAB-afkast med value-weighting og exponential decay",
    subtitle = "Baseline sammenlignet med value-weighted porteføljer med forskellige decay-parametre",
    x        = NULL,
    y        = "Kumuleret afkast (%)",
    color    = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey30")
  )

ggsave(
  filename = "robusthed_value_decay_kumuleret_afkast.png",
  plot     = p_decay,
  width    = 12,
  height   = 6,
  dpi      = 300
)

cat("Gemt: robusthed_value_decay_kumuleret_afkast.png\n")

# ------------------------------------------------------------------------------
# Plot: Sharpe-ratio
# ------------------------------------------------------------------------------

p_decay_bar <- ggplot(
  robust_decay_summary,
  aes(x = test, y = sharpe)
) +
  geom_col(width = 0.60, fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title    = "Robusthedstest: Sharpe-ratio med value-weighting og decay",
    subtitle = "Annualiseret Sharpe-ratio for BAB-strategien",
    x        = NULL,
    y        = "Sharpe-ratio"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    axis.text.x   = element_text(angle = 15, hjust = 1)
  )

ggsave(
  filename = "robusthed_value_decay_sharpe.png",
  plot     = p_decay_bar,
  width    = 10,
  height   = 5,
  dpi      = 300
)

cat("Gemt: robusthed_value_decay_sharpe.png\n")

cat("\n=== STEP 5.5 FÆRDIG — robusthedstest for value-weighted decay er gennemført ===\n")



# ==============================================================================
# SAMLET OVERSIGT: SHARPE-RATIOER FRA ALLE ROBUSTHEDSTEST
# ==============================================================================

cat("\n=== Samlet oversigt over Sharpe-ratioer fra robusthedstest ===\n")

sr_robusthed <- rbindlist(list(
  
  robust_summary[, .(
    robusthedstest = "5.1 Korrelationsvinduer",
    test,
    sharpe
  )],
  
  robust_shrink_summary[, .(
    robusthedstest = "5.2 Shrinkage",
    test,
    sharpe
  )],
  
  robust_nport_summary[, .(
    robusthedstest = "5.3 Porteføljeinddeling",
    test,
    sharpe
  )],
  
  robust_weight_summary[, .(
    robusthedstest = "5.4 Vægtningsmetode",
    test,
    sharpe
  )],
  
  robust_decay_summary[, .(
    robusthedstest = "5.5 Value-weighted decay",
    test,
    sharpe
  )]
  
), fill = TRUE)

# Afrund for pænere output
sr_robusthed[, sharpe := round(sharpe, 3)]

print(sr_robusthed)

# Gem som CSV
fwrite(sr_robusthed, "samlet_sharpe_robusthedstest.csv")

cat("\nGemt: samlet_sharpe_robusthedstest.csv\n")


# ==============================================================================
# STEP 6: EMPIRISK TEST AF HANDELSOMKOSTNINGER
# ==============================================================================
# Formål:
#   Gør diskussionen af implementerbarhed og handelsomkostninger empirisk.
#
# Metode:
#   1. Rekonstruerer BAB-porteføljens aktievægte:
#        signed_weight_i,t =
#          + w_i,t / beta_low,t   for low-beta-ben
#          - w_i,t / beta_high,t  for high-beta-ben
#
#   2. Beregner ugentlig turnover som:
#        gross_turnover_t = sum_i |signed_weight_i,t - signed_weight_i,t-1|
#
#      Dette er den samlede handlede porteføljeværdi relativt til kapital.
#      Der rapporteres også conventionel turnover:
#        conventional_turnover_t = 0.5 * gross_turnover_t
#
#   3. Beregner nettoafkast:
#        r_BAB,net,t = r_BAB,t - transaction_cost_t
#        transaction_cost_t = cost_bps / 10000 * gross_turnover_t
#
#   4. Aggregér ugentlige nettoafkast til månedlige afkast, så resultaterne
#      kan sammenlignes direkte med hovedtabellen.
#
# Output:
#   - handelsomkostninger_turnover_weekly.csv
#   - handelsomkostninger_turnover_summary.csv
#   - handelsomkostninger_performance.csv
#   - handelsomkostninger_break_even.csv
#   - handelsomkostninger_nettoafkast.png
#   - handelsomkostninger_sharpe.png
# ==============================================================================

cat("\n=== STEP 6: Empirisk test af handelsomkostninger ===\n")

library(data.table)
library(ggplot2)
library(sandwich)
library(lmtest)

# ------------------------------------------------------------------------------
# 0) Sikkerhed: objekter og nødvendige kolonner
# ------------------------------------------------------------------------------

if (!exists("data_bab")) {
  if (file.exists("data_bab.rds")) {
    data_bab <- readRDS("data_bab.rds")
    cat("Indlæst data_bab.rds\n")
  } else {
    stop("Objektet data_bab findes ikke, og data_bab.rds kunne ikke findes.")
  }
}

if (!exists("bab_clean")) {
  if (exists("bab")) {
    bab_clean <- bab[!is.na(r_bab)]
    cat("Bruger bab som bab_clean\n")
  } else if (file.exists("bab_portfolios.rds")) {
    bab_clean <- readRDS("bab_portfolios.rds")
    bab_clean <- bab_clean[!is.na(r_bab)]
    cat("Indlæst bab_portfolios.rds som bab_clean\n")
  } else {
    stop("Objektet bab_clean/bab findes ikke, og bab_portfolios.rds kunne ikke findes.")
  }
}

setDT(data_bab)
setDT(bab_clean)

required_data_bab <- c(
  "id", "iso_year", "iso_week", "date",
  "portfolio", "w", "ret_exc_wk"
)

required_bab <- c(
  "iso_year", "iso_week", "date",
  "r_bab", "r_mkt", "beta_low", "beta_high",
  "n_low", "n_high"
)

missing_data_bab <- setdiff(required_data_bab, names(data_bab))
missing_bab      <- setdiff(required_bab,      names(bab_clean))

if (length(missing_data_bab) > 0) {
  stop(
    "data_bab mangler følgende nødvendige kolonner: ",
    paste(missing_data_bab, collapse = ", ")
  )
}

if (length(missing_bab) > 0) {
  stop(
    "bab_clean mangler følgende nødvendige kolonner: ",
    paste(missing_bab, collapse = ", ")
  )
}

# Arbejd på kopier
dt_bab  <- copy(data_bab)
bab_tc  <- copy(bab_clean)

# ------------------------------------------------------------------------------
# 1) Konstruér ugentlig tidsindeks
# ------------------------------------------------------------------------------

weeks <- bab_tc[, .(
  iso_year,
  iso_week,
  date = max(date, na.rm = TRUE)
)]

setorder(weeks, iso_year, iso_week)
weeks[, week_id := .I]

bab_tc <- merge(
  bab_tc,
  weeks[, .(iso_year, iso_week, week_id)],
  by = c("iso_year", "iso_week"),
  all.x = TRUE
)

setorder(bab_tc, week_id)

# ------------------------------------------------------------------------------
# 2) Rekonstruér signed BAB-vægte på aktieniveau
# ------------------------------------------------------------------------------

# Merge porteføljebetaer ind på aktieniveau
weights_bab <- merge(
  dt_bab[
    !is.na(w) &
      !is.na(portfolio) &
      portfolio %in% c("low", "high") &
      is.finite(w),
    .(
      id,
      iso_year,
      iso_week,
      date,
      portfolio,
      w,
      ret_exc_wk
    )
  ],
  bab_tc[, .(
    iso_year,
    iso_week,
    week_id,
    beta_low,
    beta_high
  )],
  by = c("iso_year", "iso_week"),
  all.x = TRUE
)

weights_bab <- weights_bab[
  !is.na(week_id) &
    !is.na(beta_low) &
    !is.na(beta_high) &
    is.finite(beta_low) &
    is.finite(beta_high) &
    abs(beta_low) > 1e-6 &
    abs(beta_high) > 1e-6
]

# Signed, beta-skalerede vægte
weights_bab[, signed_weight := fifelse(
  portfolio == "low",
  w / beta_low,
  -w / beta_high
)]

weights_bab <- weights_bab[
  !is.na(signed_weight) &
    is.finite(signed_weight)
]

# ------------------------------------------------------------------------------
# 3) Valider at vægtene rekonstruerer BAB-afkastet
# ------------------------------------------------------------------------------

recon_bab <- weights_bab[, .(
  r_bab_reconstructed = sum(signed_weight * ret_exc_wk, na.rm = TRUE),
  gross_exposure      = sum(abs(signed_weight), na.rm = TRUE),
  net_exposure        = sum(signed_weight, na.rm = TRUE),
  n_stocks            = .N
), by = week_id]

recon_check <- merge(
  bab_tc[, .(week_id, date, r_bab, r_mkt, beta_low, beta_high, n_low, n_high)],
  recon_bab,
  by = "week_id",
  all.x = TRUE
)

recon_check[, reconstruction_error := r_bab_reconstructed - r_bab]

cat("\nValidering af rekonstruerede BAB-vægte:\n")
cat(sprintf(
  "Maksimal absolut rekonstruktionsfejl: %.10f\n",
  max(abs(recon_check$reconstruction_error), na.rm = TRUE)
))

cat(sprintf(
  "Gennemsnitlig gross exposure: %.3f\n",
  mean(recon_check$gross_exposure, na.rm = TRUE)
))

cat(sprintf(
  "Gennemsnitlig net exposure: %.3f\n",
  mean(recon_check$net_exposure, na.rm = TRUE)
))

# Gem validering
fwrite(recon_check, "handelsomkostninger_reconstruction_check.csv")

# ------------------------------------------------------------------------------
# 4) Beregn ugentlig turnover
# ------------------------------------------------------------------------------

# Aktuelle vægte
current_w <- weights_bab[, .(
  id,
  week_id,
  signed_weight
)]

# Forrige uges vægte flyttes én uge frem, så de kan sammenlignes med aktuelle vægte
previous_w <- current_w[, .(
  id,
  week_id = week_id + 1L,
  signed_weight_lag = signed_weight
)]

# Full outer merge sikrer, at både nye aktier og aktier, der forlader porteføljen,
# indgår i turnover.
turnover_long <- merge(
  current_w,
  previous_w,
  by = c("id", "week_id"),
  all = TRUE
)

# Behold kun faktiske BAB-uger
turnover_long <- turnover_long[
  week_id %in% weeks$week_id
]

turnover_long[is.na(signed_weight),     signed_weight := 0]
turnover_long[is.na(signed_weight_lag), signed_weight_lag := 0]

turnover_long[, trade_abs := abs(signed_weight - signed_weight_lag)]

turnover_weekly <- turnover_long[, .(
  gross_turnover        = sum(trade_abs, na.rm = TRUE),
  conventional_turnover = 0.5 * sum(trade_abs, na.rm = TRUE),
  n_names_traded        = sum(trade_abs > 0, na.rm = TRUE)
), by = week_id]

turnover_weekly <- merge(
  weeks,
  turnover_weekly,
  by = "week_id",
  all.x = TRUE
)

setorder(turnover_weekly, week_id)

# Første uge er initial portfolio formation. Den kan rapporteres, men ekskluderes
# fra gennemsnitlig løbende turnover.
turnover_weekly[, is_initial_week := week_id == min(week_id, na.rm = TRUE)]

turnover_summary <- turnover_weekly[is_initial_week == FALSE, .(
  avg_weekly_gross_turnover        = mean(gross_turnover, na.rm = TRUE),
  median_weekly_gross_turnover     = median(gross_turnover, na.rm = TRUE),
  avg_weekly_conventional_turnover = mean(conventional_turnover, na.rm = TRUE),
  annual_gross_turnover            = mean(gross_turnover, na.rm = TRUE) * 52,
  annual_conventional_turnover     = mean(conventional_turnover, na.rm = TRUE) * 52,
  avg_names_traded                 = mean(n_names_traded, na.rm = TRUE),
  median_names_traded              = median(n_names_traded, na.rm = TRUE)
)]

cat("\nTurnover-summary:\n")
print(turnover_summary)

fwrite(turnover_weekly,  "handelsomkostninger_turnover_weekly.csv")
fwrite(turnover_summary, "handelsomkostninger_turnover_summary.csv")

cat("\nGemt:\n")
cat("  -> handelsomkostninger_turnover_weekly.csv\n")
cat("  -> handelsomkostninger_turnover_summary.csv\n")

# ------------------------------------------------------------------------------
# 5) Break-even handelsomkostning
# ------------------------------------------------------------------------------

# Break-even beregnes som den handelsomkostning i bps, der reducerer det
# gennemsnitlige ugentlige BAB-afkast til nul:
#
# mean(r_BAB - c * gross_turnover) = 0
# c = mean(r_BAB) / mean(gross_turnover)
# bps = c * 10000

bab_turnover <- merge(
  bab_tc[, .(
    week_id,
    date,
    r_bab,
    r_mkt,
    beta_low,
    beta_high,
    n_low,
    n_high
  )],
  turnover_weekly[, .(
    week_id,
    gross_turnover,
    conventional_turnover,
    is_initial_week
  )],
  by = "week_id",
  all.x = TRUE
)

setDT(bab_turnover)

# Sikkerhed: sørg for at is_initial_week er logisk og uden NA
bab_turnover[, is_initial_week := as.logical(is_initial_week)]
bab_turnover[is.na(is_initial_week), is_initial_week := FALSE]

# Initial porteføljeopbygning sættes til 0 i omkostningsberegningen,
# så første uges etablering ikke dominerer analysen
bab_turnover[, gross_turnover_for_cost := fifelse(
  is_initial_week == TRUE,
  0,
  gross_turnover
)]

# Break-even beregnes ekskl. initial uge
break_even_bps <- bab_turnover[
  is_initial_week == FALSE &
    !is.na(r_bab) &
    !is.na(gross_turnover) &
    is.finite(r_bab) &
    is.finite(gross_turnover) &
    gross_turnover > 0,
  .(
    break_even_bps_mean_return =
      mean(r_bab, na.rm = TRUE) / mean(gross_turnover, na.rm = TRUE) * 10000,
    
    mean_weekly_bab_return_pct =
      mean(r_bab, na.rm = TRUE) * 100,
    
    mean_weekly_gross_turnover =
      mean(gross_turnover, na.rm = TRUE),
    
    mean_weekly_conventional_turnover =
      mean(conventional_turnover, na.rm = TRUE),
    
    n_weeks_used =
      .N
  )
]

cat("\nBreak-even handelsomkostning:\n")
print(break_even_bps)

fwrite(break_even_bps, "handelsomkostninger_break_even.csv")
cat("Gemt: handelsomkostninger_break_even.csv\n")

# ------------------------------------------------------------------------------
# 6) Performancefunktion for nettoafkast
# ------------------------------------------------------------------------------

calc_tc_performance <- function(dt, cost_bps) {
  
  dt <- copy(dt)
  
  # Nettoafkast efter handelsomkostninger
  dt[, r_bab_net := r_bab - (cost_bps / 10000) * gross_turnover_for_cost]
  
  dt <- dt[
    !is.na(r_bab_net) &
      !is.na(r_mkt) &
      is.finite(r_bab_net) &
      is.finite(r_mkt)
  ]
  
  # Aggreger til månedlige afkast, så tallene matcher hovedtabellen
  dt[, `:=`(
    yr  = as.integer(format(date, "%Y")),
    mth = as.integer(format(date, "%m"))
  )]
  
  dt_mth <- dt[, .(
    r_bab_net_m = prod(1 + r_bab_net, na.rm = TRUE) - 1,
    r_mkt_m     = prod(1 + r_mkt,     na.rm = TRUE) - 1,
    tc_paid_m   = sum((cost_bps / 10000) * gross_turnover_for_cost, na.rm = TRUE),
    turnover_m  = sum(gross_turnover_for_cost, na.rm = TRUE),
    date        = max(date, na.rm = TRUE)
  ), by = .(yr, mth)]
  
  dt_mth <- dt_mth[
    !is.na(r_bab_net_m) &
      !is.na(r_mkt_m) &
      is.finite(r_bab_net_m) &
      is.finite(r_mkt_m)
  ]
  
  # Mean return-test
  mod_int <- lm(r_bab_net_m ~ 1, data = dt_mth)
  
  nw_int <- coeftest(
    mod_int,
    vcov = NeweyWest(mod_int, lag = 3, prewhite = FALSE)
  )
  
  # CAPM-alpha
  mod_capm <- lm(r_bab_net_m ~ r_mkt_m, data = dt_mth)
  
  nw_capm <- coeftest(
    mod_capm,
    vcov = NeweyWest(mod_capm, lag = 3, prewhite = FALSE)
  )
  
  r <- dt_mth$r_bab_net_m
  
  mu    <- mean(r, na.rm = TRUE)
  sigma <- sd(r, na.rm = TRUE)
  
  ann_ret    <- mu * 12
  ann_vol    <- sigma * sqrt(12)
  ann_sharpe <- ann_ret / ann_vol
  
  cum_ret  <- cumprod(1 + r)
  peak     <- cummax(cum_ret)
  drawdown <- (cum_ret - peak) / peak
  max_dd   <- min(drawdown, na.rm = TRUE)
  
  data.table(
    cost_bps = cost_bps,
    n_months = nrow(dt_mth),
    start_date = min(dt_mth$date, na.rm = TRUE),
    end_date   = max(dt_mth$date, na.rm = TRUE),
    
    mean_monthly_return_pct = mu * 100,
    mean_return_t           = nw_int[1, 3],
    mean_return_p           = nw_int[1, 4],
    
    ann_return_pct = ann_ret * 100,
    ann_vol_pct    = ann_vol * 100,
    sharpe         = ann_sharpe,
    max_drawdown_pct = max_dd * 100,
    
    capm_alpha_m_pct = nw_capm[1, 1] * 100,
    capm_alpha_t     = nw_capm[1, 3],
    capm_alpha_p     = nw_capm[1, 4],
    capm_beta        = nw_capm[2, 1],
    capm_r2          = summary(mod_capm)$r.squared,
    
    avg_monthly_tc_paid_pct = mean(dt_mth$tc_paid_m, na.rm = TRUE) * 100,
    avg_monthly_turnover    = mean(dt_mth$turnover_m, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# 7) Kør handelsomkostningsscenarier
# ------------------------------------------------------------------------------

# cost_bps er én-vejs handelsomkostning pr. handlet værdi.
# Eksempel:
#   10 bps betyder 0,10% af den handlede porteføljeværdi.
cost_levels_bps <- c(0, 5, 10, 25, 50, 75, 100)

tc_performance <- rbindlist(lapply(
  cost_levels_bps,
  function(c_bps) {
    calc_tc_performance(
      dt = bab_turnover,
      cost_bps = c_bps
    )
  }
), fill = TRUE)

cat("\nPerformance efter handelsomkostninger:\n")
print(tc_performance)

fwrite(tc_performance, "handelsomkostninger_performance.csv")
saveRDS(tc_performance, "handelsomkostninger_performance.rds")

cat("\nGemt:\n")
cat("  -> handelsomkostninger_performance.csv\n")
cat("  -> handelsomkostninger_performance.rds\n")

# ------------------------------------------------------------------------------
# 8) Konstruér månedlige nettoafkastserier til plot
# ------------------------------------------------------------------------------

net_return_series <- rbindlist(lapply(cost_levels_bps, function(c_bps) {
  
  dt <- copy(bab_turnover)
  
  dt[, r_bab_net := r_bab - (c_bps / 10000) * gross_turnover_for_cost]
  
  dt[, `:=`(
    yr  = as.integer(format(date, "%Y")),
    mth = as.integer(format(date, "%m"))
  )]
  
  out <- dt[, .(
    r_bab_net_m = prod(1 + r_bab_net, na.rm = TRUE) - 1,
    date        = max(date, na.rm = TRUE)
  ), by = .(yr, mth)]
  
  out[, cost_bps := c_bps]
  out[, cost_label := paste0(c_bps, " bps")]
  
  out
  
}), fill = TRUE)

setorder(net_return_series, cost_bps, date)

net_return_series[, cum_net_return := cumprod(1 + r_bab_net_m) - 1,
                  by = cost_bps]

fwrite(net_return_series, "handelsomkostninger_nettoafkast_series.csv")
cat("Gemt: handelsomkostninger_nettoafkast_series.csv\n")

# ------------------------------------------------------------------------------
# 9) Plot: Kumulativt nettoafkast efter handelsomkostninger
# ------------------------------------------------------------------------------

p_tc_cum <- ggplot(
  net_return_series,
  aes(
    x = date,
    y = cum_net_return * 100,
    color = cost_label,
    linetype = cost_label
  )
) +
  geom_line(linewidth = 0.75) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    color = "grey40"
  ) +
  labs(
    title = "BAB-strategiens nettoafkast efter handelsomkostninger",
    subtitle = "Månedlige, kompoundede nettoafkast ved alternative én-vejs handelsomkostninger",
    x = NULL,
    y = "Kumulativt nettoafkast (%)",
    color = "Handelsomkostning",
    linetype = "Handelsomkostning"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(
  filename = "handelsomkostninger_nettoafkast.png",
  plot = p_tc_cum,
  width = 12,
  height = 6,
  dpi = 300
)

cat("Gemt: handelsomkostninger_nettoafkast.png\n")

# ------------------------------------------------------------------------------
# 10) Plot: Sharpe-ratio efter handelsomkostninger
# ------------------------------------------------------------------------------

p_tc_sharpe <- ggplot(
  tc_performance,
  aes(
    x = factor(cost_bps),
    y = sharpe
  )
) +
  geom_col(width = 0.65, fill = "steelblue") +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  labs(
    title = "BAB Sharpe-ratio efter handelsomkostninger",
    subtitle = "Sharpe-ratio beregnet på månedlige nettoafkast og annualiseret med √12",
    x = "Én-vejs handelsomkostning, bps",
    y = "Annualiseret Sharpe-ratio"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(
  filename = "handelsomkostninger_sharpe.png",
  plot = p_tc_sharpe,
  width = 9,
  height = 5,
  dpi = 300
)

cat("Gemt: handelsomkostninger_sharpe.png\n")

# ------------------------------------------------------------------------------
# 11) Plot: Annualiseret nettoafkast og CAPM-alpha
# ------------------------------------------------------------------------------

tc_plot_long <- melt(
  tc_performance[, .(
    cost_bps,
    ann_return_pct,
    capm_alpha_m_pct
  )],
  id.vars = "cost_bps",
  variable.name = "metric",
  value.name = "value"
)

tc_plot_long[, metric := fcase(
  metric == "ann_return_pct",  "Annualiseret nettoafkast",
  metric == "capm_alpha_m_pct", "CAPM-alpha, månedlig"
)]

p_tc_return_alpha <- ggplot(
  tc_plot_long,
  aes(
    x = factor(cost_bps),
    y = value,
    fill = metric
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  labs(
    title = "Afkast og CAPM-alpha efter handelsomkostninger",
    subtitle = "Nettoafkast beregnes efter ugentlig turnover og én-vejs handelsomkostninger",
    x = "Én-vejs handelsomkostning, bps",
    y = "Procent",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(
  filename = "handelsomkostninger_return_alpha.png",
  plot = p_tc_return_alpha,
  width = 10,
  height = 5.5,
  dpi = 300
)

cat("Gemt: handelsomkostninger_return_alpha.png\n")

# ------------------------------------------------------------------------------
# 12) Automatisk opsummering
# ------------------------------------------------------------------------------

base_row <- tc_performance[cost_bps == 0]
cost_25  <- tc_performance[cost_bps == 25]
cost_50  <- tc_performance[cost_bps == 50]

cat("\n=== Automatisk opsummering: handelsomkostninger ===\n")

cat(sprintf(
  "Baseline uden handelsomkostninger: Sharpe = %.3f, ann. afkast = %.2f%%.\n",
  base_row$sharpe,
  base_row$ann_return_pct
))

if (nrow(cost_25) == 1) {
  cat(sprintf(
    "Ved 25 bps handelsomkostning: Sharpe = %.3f, ann. afkast = %.2f%%.\n",
    cost_25$sharpe,
    cost_25$ann_return_pct
  ))
}

if (nrow(cost_50) == 1) {
  cat(sprintf(
    "Ved 50 bps handelsomkostning: Sharpe = %.3f, ann. afkast = %.2f%%.\n",
    cost_50$sharpe,
    cost_50$ann_return_pct
  ))
}

cat(sprintf(
  "Break-even handelsomkostning baseret på gennemsnitligt ugentligt afkast: %.2f bps.\n",
  break_even_bps$break_even_bps_mean_return
))

cat("\n=== STEP 7 FÆRDIG — handelsomkostninger testet empirisk ===\n")


# ==============================================================================
# STEP 7: EMPIRISK LEVERAGE-ANALYSE
# ==============================================================================
# Formål:
#   Kvantificerer BAB-strategiens implicitte leverage.
#
#   BAB konstrueres som:
#      r_BAB = (1 / beta_low) * r_low - (1 / beta_high) * r_high
#
#   Derfor måler:
#      1 / beta_low  = implicit leverage i long low-beta-benet
#      1 / beta_high = beta-skalering i short high-beta-benet
#
#   Analysen anvendes i diskussionen af funding constraints og praktisk
#   implementerbarhed.
#
# Output:
#   - bab_leverage_weekly.csv
#   - bab_leverage_summary.csv
#   - bab_leverage_by_decade.csv
#   - bab_leverage_timeseries.png
#   - bab_gross_exposure_histogram.png
# ==============================================================================

cat("\n=== STEP 7: Empirisk leverage-analyse ===\n")

library(data.table)
library(ggplot2)

# ------------------------------------------------------------------------------
# 0) Sikkerhed: hent BAB-serien
# ------------------------------------------------------------------------------

if (!exists("bab_clean")) {
  if (exists("bab")) {
    bab_clean <- bab[!is.na(r_bab)]
    cat("Bruger bab som bab_clean\n")
  } else if (file.exists("bab_portfolios.rds")) {
    bab_clean <- readRDS("bab_portfolios.rds")
    bab_clean <- bab_clean[!is.na(r_bab)]
    cat("Indlæst bab_portfolios.rds som bab_clean\n")
  } else {
    stop("Kan ikke finde bab_clean, bab eller bab_portfolios.rds.")
  }
}

setDT(bab_clean)

required_cols <- c(
  "date", "r_bab", "r_low", "r_high",
  "beta_low", "beta_high",
  "n_low", "n_high"
)

missing_cols <- setdiff(required_cols, names(bab_clean))

if (length(missing_cols) > 0) {
  stop(
    "bab_clean mangler følgende nødvendige kolonner: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 1) Beregn ugentlig leverage og eksponering
# ------------------------------------------------------------------------------

lev_weekly <- copy(bab_clean)

lev_weekly <- lev_weekly[
  !is.na(beta_low) &
    !is.na(beta_high) &
    is.finite(beta_low) &
    is.finite(beta_high) &
    abs(beta_low) > 1e-6 &
    abs(beta_high) > 1e-6
]

lev_weekly[, `:=`(
  leverage_low     = 1 / beta_low,
  leverage_high    = 1 / beta_high,
  gross_exposure   = (1 / beta_low) + (1 / beta_high),
  net_exposure     = (1 / beta_low) - (1 / beta_high),
  beta_spread      = beta_high - beta_low,
  total_n_stocks   = n_low + n_high,
  year             = as.integer(format(date, "%Y"))
)]

# Ekstra indikatorer
lev_weekly[, high_leverage_3x := gross_exposure > 3]
lev_weekly[, high_leverage_4x := gross_exposure > 4]
lev_weekly[, high_leverage_5x := gross_exposure > 5]

# ------------------------------------------------------------------------------
# 2) Samlet leverage-summary
# ------------------------------------------------------------------------------

lev_summary <- lev_weekly[, .(
  n_weeks = .N,
  start_date = min(date, na.rm = TRUE),
  end_date   = max(date, na.rm = TRUE),
  
  avg_beta_low  = mean(beta_low, na.rm = TRUE),
  avg_beta_high = mean(beta_high, na.rm = TRUE),
  avg_beta_spread = mean(beta_spread, na.rm = TRUE),
  
  avg_leverage_low  = mean(leverage_low, na.rm = TRUE),
  median_leverage_low = median(leverage_low, na.rm = TRUE),
  p95_leverage_low = quantile(leverage_low, 0.95, na.rm = TRUE),
  max_leverage_low = max(leverage_low, na.rm = TRUE),
  
  avg_leverage_high = mean(leverage_high, na.rm = TRUE),
  median_leverage_high = median(leverage_high, na.rm = TRUE),
  
  avg_gross_exposure = mean(gross_exposure, na.rm = TRUE),
  median_gross_exposure = median(gross_exposure, na.rm = TRUE),
  p95_gross_exposure = quantile(gross_exposure, 0.95, na.rm = TRUE),
  max_gross_exposure = max(gross_exposure, na.rm = TRUE),
  
  avg_net_exposure = mean(net_exposure, na.rm = TRUE),
  median_net_exposure = median(net_exposure, na.rm = TRUE),
  
  pct_weeks_gross_above_3x = mean(high_leverage_3x, na.rm = TRUE) * 100,
  pct_weeks_gross_above_4x = mean(high_leverage_4x, na.rm = TRUE) * 100,
  pct_weeks_gross_above_5x = mean(high_leverage_5x, na.rm = TRUE) * 100,
  
  avg_total_n_stocks = mean(total_n_stocks, na.rm = TRUE)
)]

cat("\nSamlet leverage-summary:\n")
print(lev_summary)

fwrite(lev_weekly, "bab_leverage_weekly.csv")
fwrite(lev_summary, "bab_leverage_summary.csv")

cat("\nGemt:\n")
cat("  -> bab_leverage_weekly.csv\n")
cat("  -> bab_leverage_summary.csv\n")

# ------------------------------------------------------------------------------
# 3) Leverage opdelt på delperioder
# ------------------------------------------------------------------------------

lev_weekly[, period := fifelse(
  year <= 1999, "1991-1999",
  fifelse(year <= 2009, "2000-2009",
          fifelse(year <= 2019, "2010-2019", "2020-2022"))
)]

lev_by_period <- lev_weekly[, .(
  n_weeks = .N,
  
  avg_beta_low = mean(beta_low, na.rm = TRUE),
  avg_beta_high = mean(beta_high, na.rm = TRUE),
  
  avg_leverage_low = mean(leverage_low, na.rm = TRUE),
  avg_leverage_high = mean(leverage_high, na.rm = TRUE),
  avg_gross_exposure = mean(gross_exposure, na.rm = TRUE),
  p95_gross_exposure = quantile(gross_exposure, 0.95, na.rm = TRUE),
  max_gross_exposure = max(gross_exposure, na.rm = TRUE),
  
  avg_net_exposure = mean(net_exposure, na.rm = TRUE),
  
  pct_weeks_gross_above_3x = mean(gross_exposure > 3, na.rm = TRUE) * 100,
  pct_weeks_gross_above_4x = mean(gross_exposure > 4, na.rm = TRUE) * 100,
  
  avg_r_bab_weekly_pct = mean(r_bab, na.rm = TRUE) * 100
), by = period][order(period)]

cat("\nLeverage opdelt på delperioder:\n")
print(lev_by_period)

fwrite(lev_by_period, "bab_leverage_by_period.csv")
cat("Gemt: bab_leverage_by_period.csv\n")

# ------------------------------------------------------------------------------
# 4) Identificér uger med højest gross exposure
# ------------------------------------------------------------------------------

top_leverage_weeks <- lev_weekly[order(-gross_exposure)][1:20, .(
  date,
  beta_low,
  beta_high,
  leverage_low,
  leverage_high,
  gross_exposure,
  net_exposure,
  r_bab,
  n_low,
  n_high
)]

cat("\nTop 20 uger med højest gross exposure:\n")
print(top_leverage_weeks)

fwrite(top_leverage_weeks, "bab_top_leverage_weeks.csv")
cat("Gemt: bab_top_leverage_weeks.csv\n")

# ------------------------------------------------------------------------------
# 5) Plot: leverage over tid
# ------------------------------------------------------------------------------

lev_plot <- copy(lev_weekly)

p_lev_time <- ggplot(
  lev_plot,
  aes(x = date, y = gross_exposure)
) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = mean(lev_plot$gross_exposure, na.rm = TRUE),
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 3,
             linetype = "dotted", color = "grey40") +
  labs(
    title = "BAB-strategiens implicitte gross exposure over tid",
    subtitle = "Gross exposure = 1/beta_low + 1/beta_high. Den stiplede linje viser gennemsnittet.",
    x = NULL,
    y = "Gross exposure"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(
  filename = "bab_leverage_timeseries.png",
  plot = p_lev_time,
  width = 12,
  height = 6,
  dpi = 300
)

cat("Gemt: bab_leverage_timeseries.png\n")

# ------------------------------------------------------------------------------
# 6) Plot: histogram over gross exposure
# ------------------------------------------------------------------------------

p_lev_hist <- ggplot(
  lev_plot,
  aes(x = gross_exposure)
) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = mean(lev_plot$gross_exposure, na.rm = TRUE),
             linetype = "dashed", color = "grey30") +
  labs(
    title = "Fordeling af BAB-strategiens implicitte gross exposure",
    subtitle = "Den stiplede linje viser gennemsnitlig gross exposure",
    x = "Gross exposure",
    y = "Antal uger"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(
  filename = "bab_gross_exposure_histogram.png",
  plot = p_lev_hist,
  width = 9,
  height = 5,
  dpi = 300
)

cat("Gemt: bab_gross_exposure_histogram.png\n")

# ------------------------------------------------------------------------------
# 7) Sammenhæng mellem leverage og BAB-afkast
# ------------------------------------------------------------------------------

# Spørgsmål:
#   Er høje leverage-perioder forbundet med lavere eller mere volatile BAB-afkast?

lev_weekly[, gross_exposure_lag := shift(gross_exposure, n = 1, type = "lag")]

lev_reg_data <- lev_weekly[
  !is.na(r_bab) &
    !is.na(gross_exposure_lag) &
    is.finite(r_bab) &
    is.finite(gross_exposure_lag)
]

mod_lev <- lm(r_bab ~ gross_exposure_lag, data = lev_reg_data)

nw_lev <- lmtest::coeftest(
  mod_lev,
  vcov = sandwich::NeweyWest(mod_lev, lag = 4, prewhite = FALSE)
)

cat("\nRegression: BAB-afkast på lagget gross exposure\n")
print(nw_lev)

lev_reg_summary <- data.table(
  alpha = nw_lev[1, 1],
  alpha_t = nw_lev[1, 3],
  beta_gross_exposure_lag = nw_lev[2, 1],
  beta_t = nw_lev[2, 3],
  beta_p = nw_lev[2, 4],
  r_squared = summary(mod_lev)$r.squared,
  n_obs = nobs(mod_lev)
)

fwrite(lev_reg_summary, "bab_leverage_return_regression.csv")
cat("Gemt: bab_leverage_return_regression.csv\n")

# ------------------------------------------------------------------------------
# 8) Automatisk opsummering
# ------------------------------------------------------------------------------

cat("\n=== Automatisk opsummering: leverage ===\n")

cat(sprintf(
  "Gennemsnitlig leverage i low-beta-benet: %.2fx.\n",
  lev_summary$avg_leverage_low
))

cat(sprintf(
  "Gennemsnitlig beta-skalering i high-beta-benet: %.2fx.\n",
  lev_summary$avg_leverage_high
))

cat(sprintf(
  "Gennemsnitlig gross exposure: %.2fx.\n",
  lev_summary$avg_gross_exposure
))

cat(sprintf(
  "95-percentil for gross exposure: %.2fx.\n",
  lev_summary$p95_gross_exposure
))

cat(sprintf(
  "Andel uger med gross exposure over 3x: %.1f%%.\n",
  lev_summary$pct_weeks_gross_above_3x
))

cat("\n=== STEP 7.X FÆRDIG — leverage-statistik beregnet ===\n")

# ==============================================================================
# STEP 7.X: LIKVIDITETS- OG IMPLEMENTERBARHEDSTEST
# ==============================================================================
# Formål:
#   Tester om BAB-strategien kan implementeres i et mere likvidt large-cap-univers.
#
#   Da datasættet ikke indeholder direkte likviditetsmål som volumen eller
#   bid-ask spread, anvendes market equity (me) som proxy for likviditet.
#
# Test:
#   - Full universe
#   - Top 1000 aktier efter lagget market equity
#   - Top 500 aktier efter lagget market equity
#   - Top 250 aktier efter lagget market equity
#
# Metode:
#   - Universe selection: ugentligt
#   - Top-univers bestemmes ud fra me_lag = sidste uges market equity
#   - Beta er allerede estimeret i hovedanalysen
#   - Porteføljekonstruktion: median-split, rank-weighted, beta-skaleret
#   - Performance rapporteres på månedlige, kompoundede afkast
#
# Output:
#   - likviditet_universe_summary.csv
#   - likviditet_universe_summary_formatted.csv
#   - likviditet_universe_table.png
#   - bab_liquidity_full.rds
#   - bab_liquidity_top1000.rds
#   - bab_liquidity_top500.rds
#   - bab_liquidity_top250.rds
# ==============================================================================

cat("\n=== STEP 7.X: Likviditets- og implementerbarhedstest ===\n")

library(data.table)
library(ggplot2)
library(sandwich)
library(lmtest)

# ------------------------------------------------------------------------------
# 0) Hent nødvendige data
# ------------------------------------------------------------------------------

if (!exists("data_weekly")) {
  if (file.exists("data_weekly_beta.rds")) {
    data_weekly <- readRDS("data_weekly_beta.rds")
    cat("Indlæst data_weekly_beta.rds\n")
  } else {
    stop("Kan ikke finde data_weekly eller data_weekly_beta.rds.")
  }
}

setDT(data_weekly)

required_cols <- c(
  "id", "iso_year", "iso_week", "date",
  "me", "beta_lag", "ret_exc_wk", "mkt_exc_wk"
)

missing_cols <- setdiff(required_cols, names(data_weekly))

if (length(missing_cols) > 0) {
  stop(
    "data_weekly mangler følgende kolonner: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 1) Hjælpefunktion: månedlig performance
# ------------------------------------------------------------------------------

calc_monthly_performance <- function(bab_dt, universe_label, universe_n_label) {
  
  bab_dt <- copy(bab_dt)
  
  bab_dt <- bab_dt[
    !is.na(r_bab) &
      !is.na(r_mkt) &
      is.finite(r_bab) &
      is.finite(r_mkt)
  ]
  
  bab_dt[, `:=`(
    yr  = as.integer(format(date, "%Y")),
    mth = as.integer(format(date, "%m"))
  )]
  
  bab_mth <- bab_dt[, .(
    r_bab_m = prod(1 + r_bab, na.rm = TRUE) - 1,
    r_mkt_m = prod(1 + r_mkt, na.rm = TRUE) - 1,
    beta_low  = mean(beta_low,  na.rm = TRUE),
    beta_high = mean(beta_high, na.rm = TRUE),
    n_low     = mean(n_low,     na.rm = TRUE),
    n_high    = mean(n_high,    na.rm = TRUE),
    gross_exposure = mean((1 / beta_low) + (1 / beta_high), na.rm = TRUE),
    date = max(date, na.rm = TRUE)
  ), by = .(yr, mth)]
  
  setorder(bab_mth, yr, mth)
  
  bab_mth <- bab_mth[
    !is.na(r_bab_m) &
      !is.na(r_mkt_m) &
      is.finite(r_bab_m) &
      is.finite(r_mkt_m)
  ]
  
  # Mean return-test
  mod_int <- lm(r_bab_m ~ 1, data = bab_mth)
  
  nw_int <- coeftest(
    mod_int,
    vcov = NeweyWest(mod_int, lag = 3, prewhite = FALSE)
  )
  
  # CAPM-alpha
  mod_capm <- lm(r_bab_m ~ r_mkt_m, data = bab_mth)
  
  nw_capm <- coeftest(
    mod_capm,
    vcov = NeweyWest(mod_capm, lag = 3, prewhite = FALSE)
  )
  
  r <- bab_mth$r_bab_m
  
  mu    <- mean(r, na.rm = TRUE)
  sigma <- sd(r, na.rm = TRUE)
  
  ann_return <- mu * 12
  ann_vol    <- sigma * sqrt(12)
  sharpe     <- ann_return / ann_vol
  
  cum_ret  <- cumprod(1 + r)
  peak     <- cummax(cum_ret)
  drawdown <- (cum_ret - peak) / peak
  max_dd   <- min(drawdown, na.rm = TRUE)
  
  data.table(
    universe = universe_label,
    antal_aktier = universe_n_label,
    
    n_weeks = nrow(bab_dt),
    n_months = nrow(bab_mth),
    start_date = min(bab_mth$date, na.rm = TRUE),
    end_date   = max(bab_mth$date, na.rm = TRUE),
    
    avg_n_low  = mean(bab_dt$n_low,  na.rm = TRUE),
    avg_n_high = mean(bab_dt$n_high, na.rm = TRUE),
    avg_n_total = mean(bab_dt$n_low + bab_dt$n_high, na.rm = TRUE),
    
    mean_monthly_return_pct = mu * 100,
    mean_return_t = nw_int[1, 3],
    mean_return_p = nw_int[1, 4],
    
    ann_return_pct = ann_return * 100,
    ann_vol_pct = ann_vol * 100,
    sharpe = sharpe,
    max_drawdown_pct = max_dd * 100,
    
    capm_alpha_m_pct = nw_capm[1, 1] * 100,
    capm_alpha_t = nw_capm[1, 3],
    capm_alpha_p = nw_capm[1, 4],
    capm_beta = nw_capm[2, 1],
    capm_r2 = summary(mod_capm)$r.squared,
    
    avg_beta_low = mean(bab_dt$beta_low, na.rm = TRUE),
    avg_beta_high = mean(bab_dt$beta_high, na.rm = TRUE),
    avg_gross_exposure = mean(
      (1 / bab_dt$beta_low) + (1 / bab_dt$beta_high),
      na.rm = TRUE
    )
  )
}

# ------------------------------------------------------------------------------
# 2) Hjælpefunktion: konstruér BAB for et givent univers
# ------------------------------------------------------------------------------

run_bab_liquidity_universe <- function(data_weekly_input,
                                       top_n = Inf,
                                       label = "Full universe",
                                       n_label = "Alle",
                                       save_file = NULL) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Kører univers:", label, "\n")
  cat("Antal aktier i univers:", n_label, "\n")
  cat("------------------------------------------------------------\n")
  
  dt <- copy(data_weekly_input)
  
  setorder(dt, id, iso_year, iso_week)
  
  # Lagget market equity for at undgå look-ahead bias
  dt[, me_lag := shift(me, n = 1, type = "lag"), by = id]
  
  dt <- dt[
    !is.na(beta_lag) &
      !is.na(ret_exc_wk) &
      !is.na(mkt_exc_wk) &
      !is.na(me) &
      is.finite(beta_lag) &
      is.finite(ret_exc_wk) &
      is.finite(mkt_exc_wk) &
      is.finite(me)
  ]
  
  # For Top N-univers kræves lagget market equity
  if (is.finite(top_n)) {
    
    dt <- dt[
      !is.na(me_lag) &
        is.finite(me_lag) &
        me_lag > 0
    ]
    
    # Rank aktier efter lagget market equity hver uge
    # 1 = størst
    dt[, me_rank := frank(
      -me_lag,
      ties.method = "first"
    ), by = .(iso_year, iso_week)]
    
    dt <- dt[me_rank <= top_n]
  }
  
  cat(sprintf(
    "Aktie-uge-observationer efter universfilter: %d\n",
    nrow(dt)
  ))
  
  cat(sprintf(
    "Antal uger: %d\n",
    uniqueN(dt[, paste(iso_year, iso_week)])
  ))
  
  # Minimumskrav: der skal være nok aktier til at danne low/high-ben
  dt[, n_universe_week := .N, by = .(iso_year, iso_week)]
  dt <- dt[n_universe_week >= 20]
  
  # --------------------------------------------------------------------------
  # Rank aktier efter beta_lag og lav median-split
  # --------------------------------------------------------------------------
  
  setorder(dt, iso_year, iso_week, beta_lag)
  
  dt[, rnk := rank(
    beta_lag,
    ties.method = "average"
  ), by = .(iso_year, iso_week)]
  
  dt[, z := rnk - mean(rnk), by = .(iso_year, iso_week)]
  
  dt[, beta_med := median(
    beta_lag,
    na.rm = TRUE
  ), by = .(iso_year, iso_week)]
  
  dt[, portfolio := fifelse(
    beta_lag <= beta_med,
    "low",
    "high"
  )]
  
  # Rank-weighting som i hovedanalysen
  dt[portfolio == "low",
     w := -z / sum(-z, na.rm = TRUE),
     by = .(iso_year, iso_week)]
  
  dt[portfolio == "high",
     w := z / sum(z, na.rm = TRUE),
     by = .(iso_year, iso_week)]
  
  dt <- dt[
    !is.na(w) &
      is.finite(w)
  ]
  
  # --------------------------------------------------------------------------
  # Porteføljeafkast og porteføljebeta
  # --------------------------------------------------------------------------
  
  low_port <- dt[portfolio == "low", .(
    r_low    = sum(w * ret_exc_wk, na.rm = TRUE),
    beta_low = sum(w * beta_lag,   na.rm = TRUE),
    n_low    = .N
  ), by = .(iso_year, iso_week)]
  
  high_port <- dt[portfolio == "high", .(
    r_high    = sum(w * ret_exc_wk, na.rm = TRUE),
    beta_high = sum(w * beta_lag,   na.rm = TRUE),
    n_high    = .N
  ), by = .(iso_year, iso_week)]
  
  mkt_weekly <- dt[, .(
    r_mkt = first(mkt_exc_wk),
    date  = max(date, na.rm = TRUE),
    n_universe = .N
  ), by = .(iso_year, iso_week)]
  
  bab_alt <- merge(
    low_port,
    high_port,
    by = c("iso_year", "iso_week")
  )
  
  bab_alt <- merge(
    bab_alt,
    mkt_weekly,
    by = c("iso_year", "iso_week")
  )
  
  setorder(bab_alt, iso_year, iso_week)
  
  bab_alt <- bab_alt[
    !is.na(beta_low) &
      !is.na(beta_high) &
      is.finite(beta_low) &
      is.finite(beta_high) &
      abs(beta_low) > 1e-6 &
      abs(beta_high) > 1e-6
  ]
  
  bab_alt[, r_bab := (1 / beta_low) * r_low -
            (1 / beta_high) * r_high]
  
  bab_alt[, gross_exposure := (1 / beta_low) + (1 / beta_high)]
  
  bab_alt <- bab_alt[
    !is.na(r_bab) &
      !is.na(r_mkt) &
      is.finite(r_bab) &
      is.finite(r_mkt)
  ]
  
  cat(sprintf(
    "BAB-serie konstrueret: %d uger\n",
    nrow(bab_alt)
  ))
  
  if (!is.null(save_file)) {
    saveRDS(bab_alt, save_file)
    cat("Gemt:", save_file, "\n")
  }
  
  perf <- calc_monthly_performance(
    bab_dt = bab_alt,
    universe_label = label,
    universe_n_label = n_label
  )
  
  cat(sprintf("Ann. afkast       : %+.2f%%\n", perf$ann_return_pct))
  cat(sprintf("CAPM-alpha, måned : %+.3f%%, t = %.2f\n",
              perf$capm_alpha_m_pct, perf$capm_alpha_t))
  cat(sprintf("Sharpe            : %.3f\n", perf$sharpe))
  cat(sprintf("Max drawdown      : %+.2f%%\n", perf$max_drawdown_pct))
  cat(sprintf("Gns. gross exposure: %.3f\n", perf$avg_gross_exposure))
  
  list(
    bab = bab_alt,
    performance = perf,
    data = dt
  )
}

# ------------------------------------------------------------------------------
# 3) Kør universer: Full, Top 1000, Top 500, Top 250
# ------------------------------------------------------------------------------

res_full <- run_bab_liquidity_universe(
  data_weekly_input = data_weekly,
  top_n = Inf,
  label = "Full universe",
  n_label = "Alle",
  save_file = "bab_liquidity_full.rds"
)

res_top1000 <- run_bab_liquidity_universe(
  data_weekly_input = data_weekly,
  top_n = 1000,
  label = "Top 1000",
  n_label = "1000",
  save_file = "bab_liquidity_top1000.rds"
)

res_top500 <- run_bab_liquidity_universe(
  data_weekly_input = data_weekly,
  top_n = 500,
  label = "Top 500",
  n_label = "500",
  save_file = "bab_liquidity_top500.rds"
)

res_top250 <- run_bab_liquidity_universe(
  data_weekly_input = data_weekly,
  top_n = 250,
  label = "Top 250",
  n_label = "250",
  save_file = "bab_liquidity_top250.rds"
)

# ------------------------------------------------------------------------------
# 4) Saml performance-tabel
# ------------------------------------------------------------------------------

liquidity_summary <- rbindlist(list(
  res_full$performance,
  res_top1000$performance,
  res_top500$performance,
  res_top250$performance
), fill = TRUE)

# Sæt rækkefølge
liquidity_summary[, universe := factor(
  universe,
  levels = c("Full universe", "Top 1000", "Top 500", "Top 250")
)]

setorder(liquidity_summary, universe)

liquidity_summary[, universe := as.character(universe)]

cat("\n=== Likviditets- og implementerbarhedstabel ===\n")
print(liquidity_summary)

fwrite(
  liquidity_summary,
  "likviditet_universe_summary.csv"
)

saveRDS(
  liquidity_summary,
  "likviditet_universe_summary.rds"
)

cat("\nGemt:\n")
cat("  -> likviditet_universe_summary.csv\n")
cat("  -> likviditet_universe_summary.rds\n")

# ------------------------------------------------------------------------------
# 5) Lav formateret tabel til opgaven
# ------------------------------------------------------------------------------

fmt_num <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.na(x),
    "\u2014",
    formatC(x, format = "f", digits = digits, big.mark = "")
  )
}

liquidity_table <- liquidity_summary[, .(
  Univers = universe,
  `Antal aktier` = antal_aktier,
  `Ann. afkast (%)` = fmt_num(ann_return_pct, 2),
  `CAPM-alpha (%)` = fmt_num(capm_alpha_m_pct, 2),
  `t(alpha)` = fmt_num(capm_alpha_t, 2),
  `Sharpe` = fmt_num(sharpe, 2)
)]

cat("\nFormateret tabel:\n")
print(liquidity_table)

fwrite(
  liquidity_table,
  "likviditet_universe_summary_formatted.csv"
)

cat("Gemt: likviditet_universe_summary_formatted.csv\n")

# ------------------------------------------------------------------------------
# 6) Tegn tabel som PNG
# ------------------------------------------------------------------------------

col_order <- names(liquidity_table)
plot_table <- copy(liquidity_table)

table_cells <- data.table()

for (i in seq_len(nrow(plot_table))) {
  for (j in seq_along(col_order)) {
    table_cells <- rbind(
      table_cells,
      data.table(
        row = i,
        col = j,
        col_name = col_order[j],
        value = as.character(plot_table[[col_order[j]]][i])
      )
    )
  }
}

header_cells <- data.table(
  row = 0,
  col = seq_along(col_order),
  col_name = col_order,
  value = col_order
)

all_cells <- rbind(header_cells, table_cells, fill = TRUE)

all_cells[, hjust_val := fifelse(col == 1, 0, 0.5)]
all_cells[, x_pos := col]
all_cells[col == 1, x_pos := col - 0.45]

p_table <- ggplot() +
  
  geom_rect(
    aes(
      xmin = 0.5,
      xmax = length(col_order) + 0.5,
      ymin = -0.5,
      ymax = 0.5
    ),
    fill = "#F2F2F2",
    color = NA
  ) +
  
  geom_segment(
    aes(
      x = 0.5,
      xend = length(col_order) + 0.5,
      y = seq(-0.5, nrow(plot_table) + 0.5, by = 1),
      yend = seq(-0.5, nrow(plot_table) + 0.5, by = 1)
    ),
    color = "#BDBDBD",
    linewidth = 0.3
  ) +
  
  geom_text(
    data = all_cells,
    aes(
      x = x_pos,
      y = row,
      label = value,
      hjust = hjust_val
    ),
    size = 3.6,
    fontface = ifelse(all_cells$row == 0, "bold", "plain")
  ) +
  
  scale_y_reverse(
    limits = c(nrow(plot_table) + 0.6, -0.8),
    expand = c(0, 0)
  ) +
  
  scale_x_continuous(
    limits = c(0.4, length(col_order) + 0.6),
    expand = c(0, 0)
  ) +
  
  theme_void(base_size = 12) +
  
  theme(
    plot.margin = margin(18, 18, 18, 18)
  )

ggsave(
  filename = "likviditet_universe_table.png",
  plot = p_table,
  width = 10,
  height = 3.2,
  dpi = 300
)

cat("Gemt: likviditet_universe_table.png\n")

# ------------------------------------------------------------------------------
# 7) Kumulativt afkastplot for de fire universer
# ------------------------------------------------------------------------------

cum_data <- rbindlist(list(
  res_full$bab[, .(date, r_bab, universe = "Full universe")],
  res_top1000$bab[, .(date, r_bab, universe = "Top 1000")],
  res_top500$bab[, .(date, r_bab, universe = "Top 500")],
  res_top250$bab[, .(date, r_bab, universe = "Top 250")]
), fill = TRUE)

cum_data[, universe := factor(
  universe,
  levels = c("Full universe", "Top 1000", "Top 500", "Top 250")
)]

setorder(cum_data, universe, date)

cum_data[, cum_return := cumprod(1 + r_bab) - 1, by = universe]

p_cum <- ggplot(
  cum_data,
  aes(
    x = date,
    y = cum_return * 100,
    color = universe,
    linetype = universe
  )
) +
  geom_line(linewidth = 0.75) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    color = "grey40"
  ) +
  labs(
    title = "BAB-afkast i forskellige large-cap-universer",
    subtitle = "Top-universer vælges ugentligt efter lagget market equity",
    x = NULL,
    y = "Kumuleret afkast (%)",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(
  filename = "likviditet_universe_cum_return.png",
  plot = p_cum,
  width = 12,
  height = 6,
  dpi = 300
)

cat("Gemt: likviditet_universe_cum_return.png\n")

# ------------------------------------------------------------------------------
# 8) Automatisk opsummering
# ------------------------------------------------------------------------------

cat("\n=== Automatisk opsummering ===\n")

for (i in seq_len(nrow(liquidity_summary))) {
  cat(sprintf(
    "%s: Ann. afkast = %+.2f%%, CAPM-alpha = %+.2f%%, Sharpe = %.2f, Max DD = %+.2f%%, Gross exposure = %.2f\n",
    liquidity_summary$universe[i],
    liquidity_summary$ann_return_pct[i],
    liquidity_summary$capm_alpha_m_pct[i],
    liquidity_summary$sharpe[i],
    liquidity_summary$max_drawdown_pct[i],
    liquidity_summary$avg_gross_exposure[i]
  ))
}

cat("\n=== STEP 7.X FÆRDIG — likviditets- og implementerbarhedstest gennemført ===\n")


# ==============================================================================
# BAB-FAKTOR FOR DELPERIODER — JAPAN
# ==============================================================================

library(data.table)
library(ggplot2)
library(sandwich)
library(lmtest)

cat("\n")
cat("===============================================================================\n")
cat("BAB-FAKTOR FOR DELPERIODER — JAPAN\n")
cat("===============================================================================\n\n")

# ==============================================================================
# 0) OUTPUT-FILNAVNE MED TIMESTAMP
# ==============================================================================

timestamp_bab_periods <- format(Sys.time(), "%Y%m%d_%H%M%S")

file_bab_periods_csv <- paste0(
  "bab_faktor_delperioder_",
  timestamp_bab_periods,
  ".csv"
)

file_bab_periods_png <- paste0(
  "bab_faktor_delperioder_",
  timestamp_bab_periods,
  ".png"
)

# ==============================================================================
# 1) SIKKERHEDSTJEK
# ==============================================================================

if (!exists("bab")) {
  stop("FEJL: Objektet 'bab' findes ikke. Kør først din primære BAB-kode.")
}

required_bab_cols <- c("date", "r_bab", "r_mkt")
missing_bab_cols <- setdiff(required_bab_cols, names(bab))

if (length(missing_bab_cols) > 0) {
  stop("FEJL: bab mangler kolonner: ", paste(missing_bab_cols, collapse = ", "))
}

# Tjek om antal aktier kan beregnes
has_stock_counts <- all(c("n_low", "n_high") %in% names(bab))

if (!has_stock_counts) {
  warning("bab mangler n_low og/eller n_high. Gennemsnitligt antal aktier sættes til NA.")
}

# Tjek om FF3/FF4 kan beregnes
has_ff <- exists("ff_data") && exists("ff_ok") && isTRUE(ff_ok)

if (!has_ff) {
  warning("ff_data eller ff_ok findes ikke. FF3- og FF4-alpha sættes til NA.")
}

# Arbejd kun på kopi, så det oprindelige bab-objekt ikke ændres
bab_periods <- copy(bab)
bab_periods[, date := as.IDate(date)]

cat("Sikkerhedstjek gennemført:\n")
cat("  Objektet 'bab' findes\n")
cat("  Nødvendige kolonner findes\n")
cat("  Der arbejdes kun på en kopi af bab\n")
cat(sprintf("  Antal aktier kan beregnes: %s\n", ifelse(has_stock_counts, "JA", "NEJ")))
cat(sprintf("  FF3/FF4 kan beregnes: %s\n\n", ifelse(has_ff, "JA", "NEJ")))

# ==============================================================================
# 2) AGGREGÉR UGENTLIGE BAB-AFKAST TIL MÅNEDLIGE AFKAST
# ==============================================================================
# BAB-afkastet i hovedkoden er på ugentlig frekvens.
# Her omdannes de ugentlige afkast til månedlige afkast ved compounding.
# ==============================================================================

bab_periods[, `:=`(
  yr  = as.integer(format(date, "%Y")),
  mth = as.integer(format(date, "%m"))
)]

if (has_stock_counts) {
  
  bab_mth_periods <- bab_periods[!is.na(r_bab) & !is.na(r_mkt), .(
    r_bab_m      = prod(1 + r_bab, na.rm = TRUE) - 1,
    r_mkt_m      = prod(1 + r_mkt, na.rm = TRUE) - 1,
    avg_n_stocks = mean(n_low + n_high, na.rm = TRUE),
    date         = max(date),
    n_wk         = .N
  ), by = .(yr, mth)]
  
} else {
  
  bab_mth_periods <- bab_periods[!is.na(r_bab) & !is.na(r_mkt), .(
    r_bab_m      = prod(1 + r_bab, na.rm = TRUE) - 1,
    r_mkt_m      = prod(1 + r_mkt, na.rm = TRUE) - 1,
    avg_n_stocks = NA_real_,
    date         = max(date),
    n_wk         = .N
  ), by = .(yr, mth)]
}

setorder(bab_mth_periods, yr, mth)

cat(sprintf("Månedlige BAB-observationer: %d\n", nrow(bab_mth_periods)))
cat(sprintf("Periode i månedlige data: %s til %s\n\n",
            min(bab_mth_periods$date),
            max(bab_mth_periods$date)))

# ==============================================================================
# 3) DEFINÉR DELPERIODER
# ==============================================================================

period_def <- data.table(
  periode = c("2003-2007", "2008-2010", "2011-2019", "2020-2022"),
  start   = as.IDate(c("2003-01-01", "2008-01-01", "2011-01-01", "2020-01-01")),
  end     = as.IDate(c("2007-12-31", "2010-12-31", "2019-12-31", "2022-12-31"))
)

# ==============================================================================
# 4) FUNKTION TIL AT BEREGNE BAB-PERFORMANCE FOR ÉN PERIODE
# ==============================================================================

calc_bab_period_stats <- function(dt, period_name) {
  
  dt <- copy(dt)
  dt <- dt[!is.na(r_bab_m) & !is.na(r_mkt_m)]
  
  if (nrow(dt) < 12) {
    warning("For få observationer i perioden: ", period_name)
    
    return(data.table(
      periode = period_name,
      n_måneder = nrow(dt),
      avg_n_stocks = NA_real_,
      merafkast_måned_pct = NA_real_,
      t_merafkast = NA_real_,
      capm_alpha_måned_pct = NA_real_,
      t_capm_alpha = NA_real_,
      ff3_alpha_måned_pct = NA_real_,
      t_ff3_alpha = NA_real_,
      ff4_alpha_måned_pct = NA_real_,
      t_ff4_alpha = NA_real_,
      beta_realiseret = NA_real_,
      ann_volatilitet_pct = NA_real_,
      sharpe_ratio = NA_real_
    ))
  }
  
  # --------------------------------------------------------------------------
  # Gennemsnitligt månedligt BAB-afkast med Newey-West t-statistik
  # --------------------------------------------------------------------------
  
  mod_mean <- lm(r_bab_m ~ 1, data = dt)
  nw_mean  <- coeftest(
    mod_mean,
    vcov = NeweyWest(mod_mean, lag = 3, prewhite = FALSE)
  )
  
  merafkast_måned_pct <- nw_mean[1, 1] * 100
  t_merafkast         <- nw_mean[1, 3]
  
  # --------------------------------------------------------------------------
  # CAPM-alpha: r_BAB = alpha + beta * r_mkt + fejlled
  # --------------------------------------------------------------------------
  
  mod_capm <- lm(r_bab_m ~ r_mkt_m, data = dt)
  nw_capm  <- coeftest(
    mod_capm,
    vcov = NeweyWest(mod_capm, lag = 3, prewhite = FALSE)
  )
  
  capm_alpha_måned_pct <- nw_capm[1, 1] * 100
  t_capm_alpha         <- nw_capm[1, 3]
  beta_realiseret      <- coef(mod_capm)[["r_mkt_m"]]
  
  # --------------------------------------------------------------------------
  # FF3- og FF4-alpha, hvis ff_data findes fra hovedkoden
  # --------------------------------------------------------------------------
  
  ff3_alpha_måned_pct <- NA_real_
  t_ff3_alpha         <- NA_real_
  ff4_alpha_måned_pct <- NA_real_
  t_ff4_alpha         <- NA_real_
  
  if (exists("ff_data") && exists("ff_ok") && isTRUE(ff_ok)) {
    
    dff <- merge(
      dt,
      ff_data[, .(yr, mth, mkt_rf, smb, hml, mom)],
      by = c("yr", "mth"),
      all = FALSE
    )
    
    dff <- dff[complete.cases(dff)]
    
    if (nrow(dff) >= 12) {
      
      mod_ff3 <- lm(r_bab_m ~ mkt_rf + smb + hml, data = dff)
      nw_ff3  <- coeftest(
        mod_ff3,
        vcov = NeweyWest(mod_ff3, lag = 3, prewhite = FALSE)
      )
      
      ff3_alpha_måned_pct <- nw_ff3[1, 1] * 100
      t_ff3_alpha         <- nw_ff3[1, 3]
      
      mod_ff4 <- lm(r_bab_m ~ mkt_rf + smb + hml + mom, data = dff)
      nw_ff4  <- coeftest(
        mod_ff4,
        vcov = NeweyWest(mod_ff4, lag = 3, prewhite = FALSE)
      )
      
      ff4_alpha_måned_pct <- nw_ff4[1, 1] * 100
      t_ff4_alpha         <- nw_ff4[1, 3]
    }
  }
  
  # --------------------------------------------------------------------------
  # Volatilitet og Sharpe-ratio
  # --------------------------------------------------------------------------
  
  r <- dt$r_bab_m
  
  ann_volatilitet_pct <- sd(r, na.rm = TRUE) * sqrt(12) * 100
  sharpe_ratio        <- mean(r, na.rm = TRUE) / sd(r, na.rm = TRUE) * sqrt(12)
  
  # Gennemsnitligt antal aktier i BAB-porteføljen
  avg_n_stocks_period <- mean(dt$avg_n_stocks, na.rm = TRUE)
  
  data.table(
    periode = period_name,
    n_måneder = nrow(dt),
    avg_n_stocks = avg_n_stocks_period,
    merafkast_måned_pct = merafkast_måned_pct,
    t_merafkast = t_merafkast,
    capm_alpha_måned_pct = capm_alpha_måned_pct,
    t_capm_alpha = t_capm_alpha,
    ff3_alpha_måned_pct = ff3_alpha_måned_pct,
    t_ff3_alpha = t_ff3_alpha,
    ff4_alpha_måned_pct = ff4_alpha_måned_pct,
    t_ff4_alpha = t_ff4_alpha,
    beta_realiseret = beta_realiseret,
    ann_volatilitet_pct = ann_volatilitet_pct,
    sharpe_ratio = sharpe_ratio
  )
}

# ==============================================================================
# 5) BEREGN BAB-FAKTOREN FOR ALLE DELPERIODER
# ==============================================================================

bab_period_results <- rbindlist(lapply(1:nrow(period_def), function(i) {
  
  p <- period_def[i]
  
  sub <- bab_mth_periods[
    date >= p$start & date <= p$end
  ]
  
  cat(sprintf("Beregner BAB for %s: %d måneder\n",
              p$periode,
              nrow(sub)))
  
  calc_bab_period_stats(sub, p$periode)
}))

# ==============================================================================
# 6) FORMATER RESULTATER
# ==============================================================================

bab_period_table <- copy(bab_period_results)

num_cols <- setdiff(names(bab_period_table), c("periode", "n_måneder"))

bab_period_table[, (num_cols) := lapply(.SD, function(x) round(x, 2)),
                 .SDcols = num_cols]

cat("\n")
cat("===============================================================================\n")
cat("RESULTATER — BAB-FAKTOR FOR DELPERIODER\n")
cat("===============================================================================\n\n")

print(bab_period_table)

# ==============================================================================
# 7) GEM NY CSV-FIL MED TIMESTAMP
# ==============================================================================

fwrite(bab_period_table, file_bab_periods_csv)

cat("\nGemt ny CSV-fil:\n")
cat(sprintf("  - %s\n", file_bab_periods_csv))

# ==============================================================================
# 8) LAV GRAFISK TABEL SOM PNG
# ==============================================================================

cat("\nLaver grafisk tabel som PNG...\n")

period_cols <- bab_period_table$periode
n_periods <- length(period_cols)

plot_stats <- copy(bab_period_results)

stats_list <- list(
  avg_n_stocks = plot_stats$avg_n_stocks,
  merafkast_måned_pct = plot_stats$merafkast_måned_pct,
  t_merafkast = plot_stats$t_merafkast,
  capm_alpha_måned_pct = plot_stats$capm_alpha_måned_pct,
  t_capm_alpha = plot_stats$t_capm_alpha,
  ff3_alpha_måned_pct = plot_stats$ff3_alpha_måned_pct,
  t_ff3_alpha = plot_stats$t_ff3_alpha,
  ff4_alpha_måned_pct = plot_stats$ff4_alpha_måned_pct,
  t_ff4_alpha = plot_stats$t_ff4_alpha,
  beta_realiseret = plot_stats$beta_realiseret,
  ann_volatilitet_pct = plot_stats$ann_volatilitet_pct,
  sharpe_ratio = plot_stats$sharpe_ratio
)

rdef_periods <- list(
  list(lbl = "Merafkast",          val = "merafkast_måned_pct",    tv = "t_merafkast",    t = TRUE),
  list(lbl = "CAPM-alfa",          val = "capm_alpha_måned_pct",   tv = "t_capm_alpha",   t = TRUE),
  list(lbl = "FF3-alfa",           val = "ff3_alpha_måned_pct",    tv = "t_ff3_alpha",    t = TRUE),
  list(lbl = "FF4-alfa",           val = "ff4_alpha_måned_pct",    tv = "t_ff4_alpha",    t = TRUE),
  list(lbl = "Beta (realiseret)",  val = "beta_realiseret",        tv = NA,               t = FALSE),
  list(lbl = "Volatilitet",        val = "ann_volatilitet_pct",    tv = NA,               t = FALSE),
  list(lbl = "Sharpe-ratio",       val = "sharpe_ratio",           tv = NA,               t = FALSE),
  list(lbl = "Antal aktier",       val = "avg_n_stocks",           tv = NA,               t = FALSE)
)

fmt_num_periods <- function(x, d = 2) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0) return("\u2014")
  x <- x[1]
  if (is.na(x)) "\u2014" else formatC(x, format = "f", digits = d, big.mark = "")
}

n_cells_periods <- sum(
  vapply(rdef_periods, function(r) if (isTRUE(r$t)) 2L else 1L, integer(1))
) * n_periods

gg_group_id <- integer(n_cells_periods)
gg_cn       <- integer(n_cells_periods)
gg_lbl      <- character(n_cells_periods)
gg_txt      <- character(n_cells_periods)
gg_bold     <- logical(n_cells_periods)
gg_is_t     <- logical(n_cells_periods)
gg_yp       <- numeric(n_cells_periods)

idx <- 0L
cur_y <- 0

for (ri in seq_along(rdef_periods)) {
  
  rd <- rdef_periods[[ri]]
  
  if (ri > 1) cur_y <- cur_y - 1.35
  
  y_est <- cur_y
  
  vals <- as.numeric(stats_list[[rd$val]])
  tvs  <- if (isTRUE(rd$t)) as.numeric(stats_list[[rd$tv]]) else rep(NA_real_, n_periods)
  
  for (ci in 1:n_periods) {
    
    idx <- idx + 1L
    
    val_sc <- vals[ci]
    tv_sc  <- tvs[ci]
    
    gg_group_id[idx] <- ri
    gg_cn[idx]       <- ci
    gg_lbl[idx]      <- rd$lbl
    
    if (rd$val == "avg_n_stocks") {
      gg_txt[idx] <- fmt_num_periods(val_sc, 0)
    } else {
      gg_txt[idx] <- fmt_num_periods(val_sc, 2)
    }
    
    gg_bold[idx]     <- isTRUE(rd$t) && !is.na(tv_sc) && abs(tv_sc) > 1.645
    gg_is_t[idx]     <- FALSE
    gg_yp[idx]       <- y_est
  }
  
  if (isTRUE(rd$t)) {
    
    cur_y <- cur_y - 0.55
    y_t <- cur_y
    
    for (ci in 1:n_periods) {
      
      idx <- idx + 1L
      
      tv_sc <- tvs[ci]
      t_str <- if (is.na(tv_sc)) "" else paste0("(", fmt_num_periods(tv_sc, 2), ")")
      
      gg_group_id[idx] <- ri
      gg_cn[idx]       <- ci
      gg_lbl[idx]      <- rd$lbl
      gg_txt[idx]      <- t_str
      gg_bold[idx]     <- !is.na(tv_sc) && abs(tv_sc) > 1.645
      gg_is_t[idx]     <- TRUE
      gg_yp[idx]       <- y_t
    }
  }
}

long_dt_periods <- data.table(
  group_id = gg_group_id,
  cn       = gg_cn,
  lbl      = gg_lbl,
  txt      = gg_txt,
  bold     = gg_bold,
  is_t     = gg_is_t,
  yp       = gg_yp
)

lbl_y_periods <- long_dt_periods[is_t == FALSE, .(
  ym  = mean(yp),
  lbl = unique(lbl)
), by = group_id]

setorder(lbl_y_periods, group_id)

# ==============================================================================
# 8.2) TEGN GRAFISK TABEL
# ==============================================================================

plot_long_dt_periods <- copy(long_dt_periods)
plot_lbl_y_periods   <- copy(lbl_y_periods)

Y_STRETCH <- 1.2
X_STRETCH <- 1.65   # Mere vandret afstand mellem periodekolonnerne

plot_long_dt_periods[, yp_plot := yp * Y_STRETCH]
plot_lbl_y_periods[, ym_plot := ym * Y_STRETCH]

# Giv periodekolonnerne mere afstand
plot_long_dt_periods[, x_plot := cn * X_STRETCH]

y_min <- min(plot_long_dt_periods$yp_plot)
y_max <- max(plot_long_dt_periods$yp_plot)

y_h       <- y_max + 2.00
y_h_rule  <- y_h - 0.90
y_top     <- y_h + 1.35
y_bot     <- y_min - 1.20

# Mindre venstre margen og mere plads mellem resultaterne
x_l       <- -2.05
x_r       <- n_periods * X_STRETCH + 0.85
x_label_l <- -1.90

col_text     <- "#0A0A0A"
col_tstat    <- "#5A5A5A"
col_rule     <- "#0A0A0A"
col_sub_rule <- "#BBBBBB"
col_last_bg  <- "#F7F7F7"

hdr_dt_periods <- data.table(
  cn  = 1:n_periods,
  lbl = period_cols
)

hdr_dt_periods[, x_plot := cn * X_STRETCH]

grp_sep_periods <- plot_long_dt_periods[, .(
  gmin = min(yp_plot)
), by = group_id]

setorder(grp_sep_periods, group_id)

grp_sep_periods[, ysep := gmin - 1.15]
grp_sep_periods <- grp_sep_periods[group_id < max(grp_sep_periods$group_id)]

p_tbl_periods <- ggplot() +
  
  geom_rect(
    aes(
      xmin = n_periods * X_STRETCH - 0.70,
      xmax = x_r,
      ymin = y_bot,
      ymax = y_top
    ),
    fill = col_last_bg,
    color = NA
  ) +
  
  geom_segment(
    aes(x = x_l, xend = x_r, y = y_top, yend = y_top),
    linewidth = 1.1,
    color = col_rule
  ) +
  geom_segment(
    aes(x = x_l, xend = x_r, y = y_bot, yend = y_bot),
    linewidth = 1.1,
    color = col_rule
  ) +
  
  geom_segment(
    aes(x = x_l, xend = x_r, y = y_h_rule, yend = y_h_rule),
    linewidth = 0.6,
    color = col_rule
  ) +
  
  geom_segment(
    data = grp_sep_periods,
    aes(x = x_l, xend = x_r, y = ysep, yend = ysep),
    linewidth = 0.3,
    color = col_sub_rule
  ) +
  
  geom_text(
    data = hdr_dt_periods,
    aes(x = x_plot, y = y_h, label = lbl),
    fontface = "bold",
    size = 8.3,
    color = col_text
  ) +
  
  geom_text(
    data = plot_lbl_y_periods,
    aes(x = x_label_l, y = ym_plot, label = lbl),
    hjust = 0,
    size = 7.3,
    color = col_text
  ) +
  
  geom_text(
    data = plot_long_dt_periods[is_t == FALSE & bold == FALSE],
    aes(x = x_plot, y = yp_plot, label = txt),
    size = 7.2,
    color = col_text
  ) +
  geom_text(
    data = plot_long_dt_periods[is_t == FALSE & bold == TRUE],
    aes(x = x_plot, y = yp_plot, label = txt),
    size = 7.2,
    color = col_text,
    fontface = "bold"
  ) +
  
  geom_text(
    data = plot_long_dt_periods[is_t == TRUE & bold == FALSE],
    aes(x = x_plot, y = yp_plot, label = txt),
    size = 6.3,
    color = col_tstat
  ) +
  geom_text(
    data = plot_long_dt_periods[is_t == TRUE & bold == TRUE],
    aes(x = x_plot, y = yp_plot, label = txt),
    size = 6.3,
    color = col_text,
    fontface = "bold"
  ) +
  
  scale_x_continuous(
    limits = c(x_l - 0.15, x_r + 0.15),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(y_bot - 0.60, y_top + 1.60),
    expand = c(0, 0)
  ) +
  theme_void(base_family = "") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(34, 28, 24, 28)
  ) +
  labs(
    title = NULL,
    subtitle = NULL
  )

fig_h_periods <- max(13.5, 6.5 + length(rdef_periods) * 1.10)

ggsave(
  filename = file_bab_periods_png,
  plot     = p_tbl_periods,
  width    = 17,
  height   = fig_h_periods,
  dpi      = 500
)

cat("\nGemt ny PNG-fil:\n")
cat(sprintf("  - %s\n", file_bab_periods_png))

cat("\n")
cat("===============================================================================\n")
cat("BAB-FAKTOR FOR DELPERIODER FÆRDIG\n")
cat("===============================================================================\n\n")

cat("Nye filer lavet:\n")
cat(sprintf("  - %s\n", file_bab_periods_csv))
cat(sprintf("  - %s\n", file_bab_periods_png))

cat("\nIngen eksisterende RDS-filer er ændret.\n")
cat("Ingen eksisterende PNG-filer er overskrevet, da outputfilen har timestamp.\n")
cat("Ingen eksisterende CSV-filer er overskrevet, da outputfilen har timestamp.\n")
cat("Originale objekter er ikke ændret, da analysen kun arbejder på kopier.\n")

