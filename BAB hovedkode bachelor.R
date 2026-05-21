# Datafilerne medfølger ikke i GitHub-repositoriet.

# Brugeren skal selv placere JPNall.csv og JPN_market.csv i arbejdsmappen

# eller ændre PATH_STOCKS og PATH_MARKET nedenfor.

# ==============================================================================
# BAB STRATEGI - JAPAN
# Frazzini & Pedersen (2014): "Betting Against Beta"
# ==============================================================================
# Data:
#   - JPNall.csv       : Daglige aktiedata (id, date, excntry, me, ret, ret_exc)
#   - JPN_market.csv   : Daglige markedsafkast Japan (date, ret_japan)
# ==============================================================================

library(data.table)

# ==============================================================================
# KONFIGURATION — TILPAS DISSE STIER
# ==============================================================================
# Sæt arbejdsmappen til den mappe, hvor dine data-filer ligger.
# JPNall.csv er stor (1,9 GB) — den behøver IKKE ligge i Bachelor-mappen,
# du kan blot angive den fulde sti herunder.

PATH_STOCKS <- "JPNall.csv"          # daglige aktiedata
PATH_MARKET <- "JPN_market.csv"      # daglige markedsafkast

market_summary <- market[, .(
  n_obs = .N,
  start_date = min(date),
  end_date = max(date),
  mean_daily_ret = mean(ret_japan, na.rm = TRUE),
  sd_daily_ret = sd(ret_japan, na.rm = TRUE),
  min_daily_ret = min(ret_japan, na.rm = TRUE),
  max_daily_ret = max(ret_japan, na.rm = TRUE),
  ann_return = mean(ret_japan, na.rm = TRUE) * 252,
  ann_volatility = sd(ret_japan, na.rm = TRUE) * sqrt(252)
)]

print(market_summary)
fwrite(market_summary, "jpn_market_summary.csv")

stocks_summary <- stocks[, .(
  n_obs = .N,
  n_stocks = uniqueN(id),
  start_date = min(date),
  end_date = max(date),
  mean_daily_ret = mean(ret, na.rm = TRUE),
  sd_daily_ret = sd(ret, na.rm = TRUE),
  min_daily_ret = min(ret, na.rm = TRUE),
  max_daily_ret = max(ret, na.rm = TRUE),
  mean_daily_ret_exc = mean(ret_exc, na.rm = TRUE),
  sd_daily_ret_exc = sd(ret_exc, na.rm = TRUE),
  min_daily_ret_exc = min(ret_exc, na.rm = TRUE),
  max_daily_ret_exc = max(ret_exc, na.rm = TRUE),
  mean_me = mean(me, na.rm = TRUE),
  median_me = median(me, na.rm = TRUE)
)]

print(stocks_summary)
fwrite(stocks_summary, "jpn_stocks_summary.csv")

# ==============================================================================
# STEP 0: CACHE — GENINDLÆS FORBEREGNEDE RDS-FILER
# ==============================================================================
# Hvis alle RDS-filer fra tidligere kørsler findes i arbejdsmappen, så
# springer vi de dyre trin (CSV-indlæsning, beta-estimering, portefølje-
# konstruktion) helt over og henter objekterne direkte fra disk.
#
# Sæt USE_CACHE <- FALSE hvis du vil tvinge en komplet gen-beregning.
# ==============================================================================

USE_CACHE <- TRUE

RDS_FILES <- c(
  data_daily        = "data_daily_clean.rds",
  data_weekly_clean = "data_weekly_clean.rds",
  rf_daily          = "rf_daily.rds",
  data_weekly       = "data_weekly_beta.rds",   # STEP 2 output (med beta)
  data_bab          = "data_bab.rds",
  bab               = "bab_portfolios.rds"
)

cat("\n=== STEP 0: Cache check ===\n")
rds_exists <- file.exists(RDS_FILES)
names(rds_exists) <- names(RDS_FILES)
for (i in seq_along(RDS_FILES)) {
  sz <- if (rds_exists[i]) {
    sprintf(" (%.1f MB)", file.info(RDS_FILES[i])$size / 1024^2)
  } else ""
  cat(sprintf("  %-20s %-28s %s%s\n",
              names(RDS_FILES)[i], RDS_FILES[i],
              if (rds_exists[i]) "FUNDET" else "mangler", sz))
}
USE_CACHE <- USE_CACHE && all(rds_exists)

if (USE_CACHE) {
  cat("\nAlle cache-filer findes -> indlæser fra disk (springer STEP 1-3 over)\n")
  t0 <- Sys.time()
  for (nm in names(RDS_FILES)) {
    assign(nm, readRDS(RDS_FILES[nm]))
    obj <- get(nm)
    cat(sprintf("  %-18s %10d rækker, %2d kolonner\n",
                nm, nrow(obj), ncol(obj)))
  }
  cat(sprintf("Cache indlæst på %.1f sekunder\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  
  # ---- Data-validering: tjek for fejl i de indlæste objekter ----
  cat("\n=== STEP 0.1: Data-validering ===\n")
  chk_ok <- TRUE
  chk <- function(cond, msg) {
    if (!isTRUE(cond)) {
      cat("  [FEJL] ", msg, "\n", sep = "")
      chk_ok <<- FALSE
    } else {
      cat("  [ OK ] ", msg, "\n", sep = "")
    }
  }
  
  chk(nrow(data_daily) > 1e6,             "data_daily har > 1M rækker")
  chk(all(c("id","date","ret","ret_exc","me") %in% names(data_daily)),
      "data_daily har forventede kolonner (id, date, ret, ret_exc, me)")
  chk(!anyNA(data_daily$date),            "data_daily$date er uden NA")
  chk(inherits(data_daily$date, "Date") || inherits(data_daily$date, "IDate"),
      "data_daily$date er Date/IDate")
  
  chk(nrow(data_weekly) > 1e5,            "data_weekly har > 100k rækker")
  chk("beta_lag" %in% names(data_weekly), "data_weekly indeholder beta_lag")
  chk(all(c("ret_exc_wk","iso_year","iso_week") %in% names(data_weekly)),
      "data_weekly har ret_exc_wk, iso_year, iso_week")
  if ("beta_lag" %in% names(data_weekly)) {
    bl <- data_weekly$beta_lag
    bl <- bl[!is.na(bl)]
    chk(length(bl) > 0 && median(bl) > 0.3 && median(bl) < 2,
        sprintf("beta_lag median = %.3f (forventet 0.3-2.0)", median(bl)))
  }
  
  chk(nrow(data_bab) > 1e5,               "data_bab har > 100k rækker")
  chk("decil" %in% names(data_bab) || TRUE,
      "data_bab eksisterer (decil tildeles i STEP 4)")
  if ("ret_exc_wk" %in% names(data_bab)) {
    rr <- range(data_bab$ret_exc_wk, na.rm = TRUE)
    chk(abs(rr[1]) < 2 & abs(rr[2]) < 2,
        sprintf("ret_exc_wk range = [%.3f, %.3f] (rimelig)", rr[1], rr[2]))
  }
  
  chk(nrow(bab) > 100,                    "bab har > 100 ugentlige obs")
  chk(all(c("r_bab","r_mkt") %in% names(bab)),
      "bab har r_bab og r_mkt")
  if ("r_bab" %in% names(bab)) {
    rr <- range(bab$r_bab, na.rm = TRUE)
    chk(abs(rr[1]) < 1 & abs(rr[2]) < 1,
        sprintf("r_bab range = [%.3f, %.3f] (rimelig)", rr[1], rr[2]))
  }
  
  chk(nrow(rf_daily) > 1000,              "rf_daily har > 1000 rækker")
  
  if (chk_ok) {
    cat("\nALLE CHECKS OK — cache er sund, fortsætter til STEP 4\n")
  } else {
    cat("\n!!! Fejl i cached data — anbefaler USE_CACHE <- FALSE og gen-beregning\n")
  }
} else {
  cat("\nCache ufuldstændig -> kører STEP 1-3 fra bunden\n")
}

# ==============================================================================
# STEP 1: RENS OG KLARGØR DATA
# ==============================================================================
# Formål:
#   1. Indlæs og parse rådata
#   2. Rens manglende værdier
#   3. Udled risikofri rente (rf = ret - ret_exc)
#   4. Flet aktiedata med markedsdata
#   5. Filtrer aktier med for få observationer
#   6. Konstruér ugentlige afkast til brug i volatilitetsestimering
# ==============================================================================

if (!USE_CACHE) {
  
  cat("=== STEP 1: Indlæser data ===\n")
  
  # --- 1.1 Indlæs aktiedata ---
  # Brug data.table for hastighed (28 mio. rækker)
  stocks_raw <- fread(
    PATH_STOCKS,
    colClasses = list(
      character = c("date", "excntry"),
      numeric   = c("me", "ret", "ret_exc")
    )
  )
  
  cat(sprintf("Aktiedata indlæst: %d rækker, %d aktier\n",
              nrow(stocks_raw),
              uniqueN(stocks_raw$id)))
  
  # --- 1.2 Indlæs markedsdata ---
  market_raw <- fread(PATH_MARKET)
  
  cat(sprintf("Markedsdata indlæst: %d rækker\n", nrow(market_raw)))
  
  # ==============================================================================
  # STEP 1.3: Parse datoer og grundlæggende rensning
  # ==============================================================================
  
  cat("\n=== STEP 1.3: Rensning ===\n")
  
  # Arbejd på kopi så rådata bevares uberørt
  stocks <- copy(stocks_raw)
  
  # Parse dato
  stocks[, date := as.IDate(date)]
  
  # Fjern rækker uden afkast eller markedsværdi, samt negative me
  stocks <- stocks[!is.na(ret) & !is.na(ret_exc) & !is.na(me) & me > 0]
  
  # Sorter: aktie, dernæst dato
  setorder(stocks, id, date)
  
  # Rens markedsdata
  market <- copy(market_raw)
  market[, date := as.IDate(date)]
  market <- market[!is.na(ret_japan)]
  setorder(market, date)
  
  cat(sprintf("Efter rensning: %d rækker, %d aktier\n",
              nrow(stocks),
              uniqueN(stocks$id)))
  
  # ==============================================================================
  # STEP 1.4: Udled risikofri rente
  # ==============================================================================
  # I data er:
  #   ret     = brutto dagligt afkast
  #   ret_exc = ret - rf  (excess return)
  # => rf = ret - ret_exc
  # ==============================================================================
  
  stocks[, rf := ret - ret_exc]
  
  # Én rf per dag, fælles for alle aktier (gennemsnit på tværs af aktier)
  rf_daily <- stocks[, .(rf = mean(rf, na.rm = TRUE)), by = date]
  
  cat(sprintf("Risikofri rente estimeret. Gennemsnit: %.4f%% per dag\n",
              mean(rf_daily$rf, na.rm = TRUE) * 100))
  
  # ==============================================================================
  # STEP 1.5: Flet aktiedata med markedsdata
  # ==============================================================================
  
  cat("\n=== STEP 1.5: Fletter aktie- og markedsdata ===\n")
  
  # Flet markedsafkast ind
  data_daily <- market[stocks, on = "date"]
  
  # Flet daglig rf ind
  data_daily <- rf_daily[data_daily, on = "date"]
  
  # Beregn markedets excess return og behold kun regnte kolonner
  data_daily[, mkt_exc := ret_japan - rf]
  data_daily <- data_daily[, .(id, date, me, ret, ret_exc, rf, ret_japan, mkt_exc)]
  
  # Tjek for manglende markedsdata
  n_missing_mkt <- sum(is.na(data_daily$ret_japan))
  cat(sprintf("Rækker uden markedsdata: %d (%.1f%%)\n",
              n_missing_mkt,
              100 * n_missing_mkt / nrow(data_daily)))
  
  # Fjern rækker uden markedsafkast
  data_daily <- data_daily[!is.na(ret_japan)]
  
  cat(sprintf("Endelig daglig panel: %d rækker, %d aktier\n",
              nrow(data_daily),
              uniqueN(data_daily$id)))
  
  # ==============================================================================
  # STEP 1.6: Filtrer aktier med for få observationer
  # ==============================================================================
  # F&P kræver minimum 120 daglige obs. i 1-årsvindue til korrelation.
  # Vi beholder kun aktier med mindst 252 daglige observationer totalt.
  # ==============================================================================
  
  cat("\n=== STEP 1.6: Filtrerer aktier med for få obs ===\n")
  
  min_obs      <- 252
  obs_per_stock <- data_daily[, .N, by = id]
  active_stocks <- obs_per_stock[N >= min_obs, id]
  data_daily    <- data_daily[id %in% active_stocks]
  
  cat(sprintf("Aktier med >= %d obs.: %d (fjernet %d)\n",
              min_obs,
              uniqueN(data_daily$id),
              uniqueN(stocks$id) - uniqueN(data_daily$id)))
  
  # ==============================================================================
  # STEP 1.6c: Winsorization af ekstreme daglige excess returns
  # ==============================================================================
  # Formål:
  #   Håndterer ekstreme afkast konsistent i hele sampleperioden.
  #   Dette erstatter ad hoc-filtrering i enkelte delperioder.
  #   Afkast winsoriseres månedligt ved 0.1% og 99.9%.
  # ==============================================================================
  
  cat("\n=== STEP 1.6c: Winsorization af ekstreme afkast ===\n")
  
  data_daily[, yr  := as.integer(format(date, "%Y"))]
  data_daily[, mth := as.integer(format(date, "%m"))]
  
  winsorize_vec <- function(x, p_low = 0.001, p_high = 0.999) {
    q <- quantile(x, probs = c(p_low, p_high), na.rm = TRUE, type = 7)
    pmin(pmax(x, q[1]), q[2])
  }
  
  data_daily[, ret_exc_raw := ret_exc]
  data_daily[, ret_raw     := ret]
  
  data_daily[, ret_exc := winsorize_vec(ret_exc), by = .(yr, mth)]
  data_daily[, ret     := winsorize_vec(ret),     by = .(yr, mth)]
  
  # Dokumentér effekten
  winsor_summary <- data_daily[, .(
    n_obs = .N,
    n_ret_exc_changed = sum(ret_exc != ret_exc_raw, na.rm = TRUE),
    pct_ret_exc_changed = 100 * mean(ret_exc != ret_exc_raw, na.rm = TRUE),
    min_ret_exc_raw = min(ret_exc_raw, na.rm = TRUE),
    max_ret_exc_raw = max(ret_exc_raw, na.rm = TRUE),
    min_ret_exc_wins = min(ret_exc, na.rm = TRUE),
    max_ret_exc_wins = max(ret_exc, na.rm = TRUE)
  )]
  
  print(winsor_summary)
  fwrite(winsor_summary, "winsorization_summary.csv")
  
  # Tjek rå aktiedata før rensning
  
  summary(stocks_raw$ret_exc)
  
  summary(stocks_raw$ret)
  
  # Find meget ekstreme observationer i rådata
  
  stocks_raw[ret_exc > 1 | ret > 1, .(id, date, ret, ret_exc, me)][1:20]
  
  # Tjek efter rensning, men før winsorization
  
  data_daily[ret_exc_raw > 1 | ret_raw > 1, .(id, date, ret_raw, ret_exc_raw, me)][1:20]
  
  # Tjek om ekstremen ligger i ret og ikke ret_exc
  
  data_daily[ret_raw > 1, .(id, date, ret_raw, ret_exc_raw, me)][1:20]
  
  # ==============================================================================
  # STEP 1.7: Konstruér ugentlige afkast
  # ==============================================================================
  # Vi komprimerer daglige afkast til ugentlige ved at sammensætte dem:
  #   r_wk = prod(1 + r_daily) - 1
  # Ugen identificeres via ISO-ugenummer og ISO-år,
  # så uger der krydser årsskiftet håndteres korrekt.
  # Markedsværdi sættes til ultimo-værdien for ugen (fredag).
  # ==============================================================================
  
  cat("\n=== STEP 1.7: Konstruerer ugentlige afkast ===\n")
  
  # Tilføj ISO-uge og ISO-år som nye kolonner
  # %G = ISO 8601 ugebaseret år (håndterer årsskiftet korrekt uden lubridate)
  data_daily[, iso_year := as.integer(format(date, "%G"))]
  data_daily[, iso_week := isoweek(date)]
  
  # Aggregér til ugentligt panel
  data_weekly <- data_daily[, .(
    ret_wk     = prod(1 + ret,       na.rm = TRUE) - 1,
    ret_exc_wk = prod(1 + ret_exc,   na.rm = TRUE) - 1,
    ret_mkt_wk = prod(1 + ret_japan, na.rm = TRUE) - 1,
    rf_wk      = prod(1 + rf,        na.rm = TRUE) - 1,
    mkt_exc_wk = prod(1 + mkt_exc,   na.rm = TRUE) - 1,
    me         = last(me),           # ultimo markedsværdi for ugen
    date       = max(date),          # fredag (eller sidste handelsdag)
    n_days     = .N
  ), by = .(id, iso_year, iso_week)]
  
  setorder(data_weekly, id, iso_year, iso_week)
  
  cat(sprintf("Ugentlig panel: %d rækker, %d aktier, %d uger\n",
              nrow(data_weekly),
              uniqueN(data_weekly$id),
              uniqueN(data_weekly[, paste(iso_year, iso_week)])))
  
  # ==============================================================================
  # STEP 1.8: Overblik og validering
  # ==============================================================================
  
  cat("\n=== STEP 1.8: Validering ===\n")
  
  cat(sprintf("Daglig data periode  : %s til %s\n",
              min(data_daily$date), max(data_daily$date)))
  cat(sprintf("Ugentlig data periode: %d uge %d til %d uge %d\n",
              data_weekly[iso_year == min(iso_year), min(iso_year)],
              data_weekly[iso_year == min(iso_year), min(iso_week)],
              data_weekly[iso_year == max(iso_year), max(iso_year)],
              data_weekly[iso_year == max(iso_year), max(iso_week)]))
  
  # Tjek for duplikater
  dupl_daily <- sum(duplicated(data_daily, by = c("id", "date")))
  cat(sprintf("Duplikerede (id, date) i daglig panel: %d\n", dupl_daily))
  
  # Afkast-fordeling
  cat("\nSammenfatning af daglige excess afkast:\n")
  print(summary(data_daily$ret_exc))
  
  cat("\nSammenfatning af daglige markedsafkast:\n")
  print(summary(data_daily$mkt_exc))
  
  # ==============================================================================
  # STEP 1.9: Gem rensede datasæt
  # ==============================================================================
  
  cat("\n=== STEP 1.9: Gemmer rensede datasæt ===\n")
  
  saveRDS(data_daily,  "data_daily_clean.rds")
  saveRDS(data_weekly, "data_weekly_clean.rds")
  saveRDS(rf_daily,    "rf_daily.rds")
  
  cat("Gemt:\n")
  cat("  -> data_daily_clean.rds\n")
  cat("  -> data_weekly_clean.rds\n")
  cat("  -> rf_daily.rds\n")
  cat("\n=== STEP 1 FÆRDIG ===\n")
  cat("Klar til Step 2: Beta-estimering (Frazzini & Pedersen metode)\n")
  
} # end if (!USE_CACHE) — STEP 1

# ==============================================================================
# STEP 2: BETA-ESTIMERING (Frazzini & Pedersen metode — ugentlig)
# ==============================================================================
# Alt beregnes på UGENTLIGE afkast for at undgå den lange køretid ved daglig data.
#
#   β_i = ρ_{i,m} × (σ_i / σ_m)
#
#   ρ : rolling korrelation fra ugentlige afkast
#       vindue = 260 uger (5 år), minimum 52 uger
#
#   σ : rolling volatilitet fra ugentlige afkast
#       vindue = 52 uger (1 år), minimum 26 uger 
#
# Shrinkage mod 1 (F&P):
#   β̂_i = 0.6 × β_raw + 0.4
#
# LOOK-AHEAD BIAS:
#   Beta beregnet i uge t bruger data op til og med uge t.
#   Til porteføljedannelse bruges beta_lag = shift(beta, 1) — altså
#   beta fra uge t-1 — så vi aldrig handler på viden fra den aktuelle uge.
# ==============================================================================

if (!USE_CACHE) {
  
  cat("\n=== STEP 2: Beta-estimering (ugentlig) ===\n")
  
  # Sikr sortering inden rolling beregninger
  setorder(data_weekly, id, iso_year, iso_week)
  
  # ==============================================================================
  # STEP 2.1: Rolling korrelation (ugentlig, 52 uger, min 26)
  # ==============================================================================
  
  cat("Beregner ugentlig rolling korrelation (52 uger)...\n")
  
  data_weekly[, rho := frollapply(
    seq_len(.N), n = 260,
    FUN = function(idx) {
      xi <- ret_exc_wk[idx]
      xm <- mkt_exc_wk[idx]
      if (sum(!is.na(xi) & !is.na(xm)) < 156) return(NA_real_)
      cor(xi, xm, use = "complete.obs")
    },
    fill = NA
  ), by = id]
  
  cat(sprintf("Korrelation beregnet. Andel med gyldig rho: %.1f%%\n",
              100 * mean(!is.na(data_weekly$rho))))
  
  # ==============================================================================
  # STEP 2.2: Rolling volatilitet (ugentlig, 260 uger, min 156)
  # ==============================================================================
  
  cat("Beregner ugentlig rolling volatilitet (260 uger)...\n")
  
  # Aktie-volatilitet (per aktie)
  data_weekly[, sigma_i := frollapply(
    ret_exc_wk, n = 52,
    FUN = function(x) if (sum(!is.na(x)) < 26) NA_real_ else sd(x, na.rm = TRUE),
    fill = NA
  ), by = id]
  
  # Markeds-volatilitet — samme for alle aktier, beregnes kun én gang
  # Brug first() for at sikre præcis én række per uge (undgår cartesian join-fejl)
  mkt_sigma <- data_weekly[, .(mkt_exc_wk = first(mkt_exc_wk)),
                           by = .(iso_year, iso_week)]
  setorder(mkt_sigma, iso_year, iso_week)
  mkt_sigma[, sigma_m := frollapply(
    mkt_exc_wk, n = 52,
    FUN = function(x) if (sum(!is.na(x)) < 26) NA_real_ else sd(x, na.rm = TRUE),
    fill = NA
  )]
  
  # Flet sigma_m ind — nu garanteret én række per uge i mkt_sigma
  data_weekly <- mkt_sigma[, .(iso_year, iso_week, sigma_m)][
    data_weekly, on = .(iso_year, iso_week)
  ]
  
  cat(sprintf("Volatilitet beregnet. Andel med gyldig sigma_i: %.1f%%\n",
              100 * mean(!is.na(data_weekly$sigma_i))))
  
  # ==============================================================================
  # STEP 2.3: Beregn raw beta og shrinkage
  # ==============================================================================
  
  cat("Beregner beta og shrinkage...\n")
  
  data_weekly[, beta_raw := rho * (sigma_i / sigma_m)]
  data_weekly[, beta     := 0.6 * beta_raw + 0.4]
  
  # ==============================================================================
  # STEP 2.4: Lag beta én uge for at undgå look-ahead bias
  # ==============================================================================
  # Beta estimeret i uge t må IKKE bruges til at handle i uge t.
  # Vi lagger beta én uge: beta_lag er hvad vi kender ved starten af uge t,
  # og det er den vi bruger til porteføljedannelse i Step 3.
  # ==============================================================================
  
  data_weekly[, beta_lag := shift(beta, n = 1, type = "lag"), by = id]
  
  cat(sprintf("Look-ahead bias korrigeret: bruger beta fra uge t-1 til handel i uge t\n"))
  
  # ==============================================================================
  # STEP 2.5: Validering
  # ==============================================================================
  
  cat("\n=== STEP 2.5: Validering af beta ===\n")
  
  cat(sprintf("Rækker med gyldig beta_lag: %d (%.1f%%)\n",
              sum(!is.na(data_weekly$beta_lag)),
              100 * mean(!is.na(data_weekly$beta_lag))))
  
  cat("\nFordeling af beta_lag (shrunk, lagget):\n")
  print(summary(data_weekly$beta_lag))
  
  cat("\nFordeling af beta_raw:\n")
  print(summary(data_weekly$beta_raw))
  
  # Vægtet gennemsnitsbeta bør ligge tæt på 1
  cat(sprintf("\nVægtet gns. beta_lag (value-weighted): %.3f\n",
              data_weekly[!is.na(beta_lag) & !is.na(me),
                          sum(beta_lag * me) / sum(me)]))
  
  # ==============================================================================
  # STEP 2.6: Gem
  # ==============================================================================
  
  saveRDS(data_weekly, "data_weekly_beta.rds")
  
  cat("\nGemt: -> data_weekly_beta.rds\n")
  cat("\n=== STEP 2 FÆRDIG ===\n")
  cat("Klar til Step 3: BAB porteføljekonstruktion\n")
  
} # end if (!USE_CACHE) — STEP 2

# ==============================================================================
# STEP 3: BAB PORTEFØLJEKONSTRUKTION (Frazzini & Pedersen metode)
# ==============================================================================
# Hver uge rangeres aktier efter beta_lag (beta fra forrige uge).
# Vi danner to porteføljer:
#   - Lav-beta (L): aktier med beta <= median denne uge
#   - Høj-beta (H): aktier med beta >  median denne uge
#
# Vægtning (rank-vægtet som i F&P):
#   z_i = rank_i - mean(rank)          (centreret rang)
#   w_L_i = -z_i / Σ(-z_j)  for j i L  (lav beta = lav rang = negativ z → flip)
#   w_H_i =  z_i / Σ( z_j)  for j i H  (høj beta = høj rang = positiv z)
#
# BAB-afkast (skaleret til beta = 1 på begge ben):
#   r_BAB = (1/β_L) × r_L − (1/β_H) × r_H
#
# Timing (ingen look-ahead bias):
#   beta_lag[t] dannet ved udgangen af uge t-1
#   ret_exc_wk[t] er afkastet tjent I uge t
#   → vi handler på beta_lag og modtager ret_exc_wk samme uge ✓
# ==============================================================================

if (!USE_CACHE) {
  
  cat("\n=== STEP 3: BAB porteføljekonstruktion ===\n")
  
  # --- 3.1 Klargør arbejdsdata ---
  # Behold kun uger hvor vi har gyldig lagget beta og ugentligt afkast
  data_bab <- data_weekly[!is.na(beta_lag) & !is.na(ret_exc_wk) & !is.na(me)]
  
  cat(sprintf("Arbejdsdata: %d rækker, %d aktier, %d uger\n",
              nrow(data_bab),
              uniqueN(data_bab$id),
              uniqueN(data_bab[, paste(iso_year, iso_week)])))
  
  # --- 3.2 Rangér og beregn rank-vægte per uge ---
  # rank() per uge — ties.method = "average" giver gennemsnitsrang ved bindinger
  data_bab[, rnk      := rank(beta_lag, ties.method = "average"), by = .(iso_year, iso_week)]
  data_bab[, z        := rnk - mean(rnk),                         by = .(iso_year, iso_week)]
  data_bab[, beta_med := median(beta_lag),                         by = .(iso_year, iso_week)]
  
  # Porteføljetildeling
  data_bab[, portfolio := fifelse(beta_lag <= beta_med, "low", "high")]
  
  # Rank-vægte inden for hver portefølje
  # Lav-beta:  w ∝ −z  (lav rang → negativ z → positiv vægt efter flip)
  # Høj-beta:  w ∝  z  (høj rang → positiv z)
  data_bab[portfolio == "low",
           w := -z / sum(-z),
           by = .(iso_year, iso_week)]
  
  data_bab[portfolio == "high",
           w :=  z / sum( z),
           by = .(iso_year, iso_week)]
  
  # --- 3.3 Beregn porteføljeafkast og porteføljiebeta per uge ---
  low_port <- data_bab[portfolio == "low", .(
    r_low    = sum(w * ret_exc_wk, na.rm = TRUE),
    beta_low = sum(w * beta_lag,   na.rm = TRUE),
    n_low    = .N
  ), by = .(iso_year, iso_week)]
  
  high_port <- data_bab[portfolio == "high", .(
    r_high    = sum(w * ret_exc_wk, na.rm = TRUE),
    beta_high = sum(w * beta_lag,   na.rm = TRUE),
    n_high    = .N
  ), by = .(iso_year, iso_week)]
  
  # Markedsafkast per uge (til sammenligning)
  mkt_weekly <- data_bab[, .(
    r_mkt = first(mkt_exc_wk)
  ), by = .(iso_year, iso_week)]
  
  # --- 3.4 Sammensæt BAB-portefølje ---
  bab <- merge(low_port, high_port, by = c("iso_year", "iso_week"))
  bab <- merge(bab,      mkt_weekly, by = c("iso_year", "iso_week"))
  setorder(bab, iso_year, iso_week)
  
  # BAB-afkast: skaler hvert ben til beta = 1
  bab[, r_bab := (1 / beta_low) * r_low - (1 / beta_high) * r_high]
  
  # Tilføj dato (fredag/sidste handelsdag i ugen) til tidsserieplot
  bab <- data_bab[, .(date = max(date)), by = .(iso_year, iso_week)][bab, on = .(iso_year, iso_week)]
  
  cat(sprintf("BAB-portefølje konstrueret: %d uger\n", nrow(bab)))
  
  # ==============================================================================
  # STEP 3.5: Validering
  # ==============================================================================
  
  cat("\n=== STEP 3.5: Validering ===\n")
  
  cat(sprintf("Gns. antal aktier i lav-beta ben : %.0f\n", mean(bab$n_low)))
  cat(sprintf("Gns. antal aktier i høj-beta ben : %.0f\n", mean(bab$n_high)))
  cat(sprintf("Gns. beta_low                    : %.3f\n", mean(bab$beta_low,  na.rm = TRUE)))
  cat(sprintf("Gns. beta_high                   : %.3f\n", mean(bab$beta_high, na.rm = TRUE)))
  
  # Tjek at lav-beta faktisk er lavere end høj-beta
  stopifnot(mean(bab$beta_low) < mean(bab$beta_high))
  cat("OK: beta_low < beta_high i gennemsnit ✓\n")
  
  cat(sprintf("\nBAB ugentligt gns. afkast : %.4f%%\n", mean(bab$r_bab, na.rm = TRUE) * 100))
  cat(sprintf("BAB annualiseret afkast   : %.2f%%\n",  mean(bab$r_bab, na.rm = TRUE) * 52 * 100))
  cat(sprintf("BAB annualiseret Sharpe   : %.3f\n",
              (mean(bab$r_bab, na.rm = TRUE) / sd(bab$r_bab, na.rm = TRUE)) * sqrt(52)))
  
  # ==============================================================================
  # STEP 3.6: Gem
  # ==============================================================================
  
  saveRDS(bab,      "bab_portfolios.rds")
  saveRDS(data_bab, "data_bab.rds")
  
  cat("\nGemt:\n")
  cat("  -> bab_portfolios.rds  (ugentlige BAB-afkast)\n")
  cat("  -> data_bab.rds        (aktieniveau med vægte)\n")
  cat("\n=== STEP 3 FÆRDIG ===\n")
  cat("Klar til Step 4: Performance analyse og resultater\n")
  
} # end if (!USE_CACHE) — STEP 3

# ==============================================================================
# STEP 4: PERFORMANCE ANALYSE OG RESULTATER
# ==============================================================================
# Indhold:
#   4.1  Newey-West t-test: er BAB-afkastet signifikant forskelligt fra nul?
#   4.2  CAPM-alfa: er afkastet signifikant efter korrektion for markedsrisiko?
#   4.3  Nøgletal: afkast, Sharpe, Sortino, max drawdown
#   4.4  Decilporteføljer: afkast og beta for hvert betadecil (som F&P Table 3)
#   4.5  Visualisering: kumuleret afkast + decilplot
# ==============================================================================

library(sandwich)   # Newey-West kovariansmatrix
library(lmtest)     # coeftest()
library(ggplot2)

cat("\n=== STEP 4: Performance analyse ===\n")

# ==============================================================================
# STEP 4.1: Newey-West t-test på BAB-afkast
# ==============================================================================
# Vi tester H0: E[r_BAB] = 0 med Newey-West standardfejl (lag = 4 uger)
# for at korrigere for eventuel autokorrelation og heteroskedasticitet.
# ==============================================================================

cat("\n--- 4.1: Newey-West t-test ---\n")

bab_clean <- bab[!is.na(r_bab)]

# Interceptmodel: r_BAB ~ 1
mod_int <- lm(r_bab ~ 1, data = bab_clean)
nw_int  <- coeftest(mod_int, vcov = NeweyWest(mod_int, lag = 4))

cat(sprintf("Gns. ugentligt BAB-afkast : %+.4f%%\n", coef(mod_int) * 100))
cat(sprintf("Newey-West t-statistik    : %.3f\n",     nw_int[1, 3]))
cat(sprintf("p-værdi (2-sidet)         : %.4f\n",     nw_int[1, 4]))
cat(sprintf("Signifikant (5%%)          : %s\n",
            ifelse(nw_int[1, 4] < 0.05, "JA", "NEJ")))

# ==============================================================================
# STEP 4.2: CAPM-alfa (Newey-West)
# ==============================================================================
# r_BAB = α + β_mkt × r_mkt + ε
# En positiv og signifikant α viser at BAB ikke blot er markedseksponering.
# ==============================================================================

cat("\n--- 4.2: CAPM-alfa ---\n")

mod_capm <- lm(r_bab ~ r_mkt, data = bab_clean)
nw_capm  <- coeftest(mod_capm, vcov = NeweyWest(mod_capm, lag = 4))

cat(sprintf("CAPM alfa (ugentlig)      : %+.4f%%\n", nw_capm[1, 1] * 100))
cat(sprintf("Alfa t-statistik          : %.3f\n",     nw_capm[1, 3]))
cat(sprintf("Alfa p-værdi              : %.4f\n",     nw_capm[1, 4]))
cat(sprintf("Markeds-beta (β_mkt)      : %.3f\n",     nw_capm[2, 1]))
cat(sprintf("R²                        : %.3f\n",     summary(mod_capm)$r.squared))

# ==============================================================================
# STEP 4.3: Nøgletal
# ==============================================================================

cat("\n--- 4.3: Nøgletal (annualiseret, 52 uger) ---\n")

r     <- bab_clean$r_bab
mu    <- mean(r)
sigma <- sd(r)

# Annualisering
ann_ret    <- mu    * 52
ann_vol    <- sigma * sqrt(52)
ann_sharpe <- ann_ret / ann_vol

# Sortino (nedadgående volatilitet)
downside   <- r[r < 0]
ann_sortino <- ann_ret / (sd(downside) * sqrt(52))

# Maximum drawdown
cum_ret  <- cumprod(1 + r)
peak     <- cummax(cum_ret)
drawdown <- (cum_ret - peak) / peak
max_dd   <- min(drawdown)

# Hit rate: andel af uger med positivt afkast
hit_rate <- mean(r > 0)

cat(sprintf("Annualiseret afkast  : %+.2f%%\n", ann_ret    * 100))
cat(sprintf("Annualiseret vol.    :  %.2f%%\n", ann_vol    * 100))
cat(sprintf("Sharpe ratio         :  %.3f\n",   ann_sharpe))
cat(sprintf("Sortino ratio        :  %.3f\n",   ann_sortino))
cat(sprintf("Max drawdown         : %+.2f%%\n", max_dd     * 100))
cat(sprintf("Hit rate             :  %.1f%%\n", hit_rate   * 100))

# ==============================================================================
# STEP 4.4: Decilporteføljer (som F&P Table 3)
# ==============================================================================
# Aktier rangeres i 10 deciler efter beta_lag hver uge.
# Vi beregner equal-weighted afkast og gennemsnitsbeta per decil.
# ==============================================================================

cat("\n--- 4.4: Decilporteføljer ---\n")

# Brug ntile() via rank — giver direkte integer-deciler uden factor-problemer
data_bab[, decil := as.integer(cut(
  rank(beta_lag, ties.method = "first"),
  breaks = quantile(rank(beta_lag, ties.method = "first"),
                    probs = 0:10/10, na.rm = TRUE),
  labels = 1:10,
  include.lowest = TRUE
)), by = .(iso_year, iso_week)]

decil_ts <- data_bab[!is.na(decil), .(
  r_decil   = mean(ret_exc_wk, na.rm = TRUE),
  beta_ante = mean(beta_lag,   na.rm = TRUE)
), by = .(decil, iso_year, iso_week)]

decil_summary <- decil_ts[, .(
  ann_ret  = mean(r_decil, na.rm = TRUE) * 52 * 100,
  beta_avg = mean(beta_ante, na.rm = TRUE)
), by = decil][order(as.integer(decil))]

cat("\nDecil  Beta    Ann. afkast\n")
cat("-------------------------------\n")
decil_summary[, cat(sprintf("  D%02d   %.3f   %+.2f%%\n",
                            as.integer(decil), beta_avg, ann_ret),
                    sep = ""), by = decil]

# ==============================================================================
# STEP 4.5: Visualisering
# ==============================================================================

cat("\n--- 4.5: Gemmer grafer ---\n")

# --- Graf 1: Kumuleret afkast for BAB, lav-beta og høj-beta ---
bab_clean[, `:=`(
  cum_bab  = cumprod(1 + r_bab)  - 1,
  cum_low  = cumprod(1 + r_low)  - 1,
  cum_high = cumprod(1 + r_high) - 1
)]

plot_data <- melt(
  bab_clean[, .(date, cum_bab, cum_low, cum_high)],
  id.vars       = "date",
  variable.name = "Portefølje",
  value.name    = "Kumuleret afkast"
)

plot_data[, Portefølje := fcase(
  Portefølje == "cum_bab",  "BAB",
  Portefølje == "cum_low",  "Lav beta",
  Portefølje == "cum_high", "Høj beta"
)]

p1 <- ggplot(
  plot_data,
  aes(x = date, y = `Kumuleret afkast` * 100,
      color = Portefølje, linetype = Portefølje)
) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c(
    "BAB"      = "firebrick",
    "Lav beta" = "steelblue",
    "Høj beta" = "steelblue4"
  )) +
  scale_linetype_manual(values = c(
    "BAB"      = "solid",
    "Lav beta" = "solid",
    "Høj beta" = "solid"
  )) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title    = "Kumuleret afkast — BAB strategi Japan",
    subtitle = "Ugentlige excess returns for BAB samt lav- og høj-beta-ben",
    x        = NULL,
    y        = "Kumuleret afkast (%)",
    color    = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("BAB_kumuleret_afkast.png", plot = p1, width = 12, height = 6, dpi = 300)
cat("Gemt: BAB_kumuleret_afkast.png\n")

# --- Graf 2: Annualiseret afkast per betadecil ---
p2 <- ggplot(decil_summary, aes(x = as.integer(decil), y = ann_ret)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = 1:10, labels = paste0("D", 1:10)) +
  labs(
    title    = "Annualiseret afkast per betadecil — Japan",
    subtitle = "D1 = laveste beta, D10 = højeste beta (equal-weighted)",
    x        = "Beta-decil",
    y        = "Annualiseret afkast (%)"
  ) +
  theme_minimal(base_size = 12)

ggsave("BAB_decil_afkast.png", plot = p2, width = 10, height = 5, dpi = 300)
cat("Gemt: BAB_decil_afkast.png\n")

# ==============================================================================
# STEP 4.6: Samlede resultater
# ==============================================================================

cat("\n=== SAMLEDE RESULTATER ===\n")
cat(sprintf("Periode              : %s til %s\n", min(bab_clean$date), max(bab_clean$date)))
cat(sprintf("Antal uger           : %d\n",         nrow(bab_clean)))
cat(sprintf("Ann. BAB-afkast      : %+.2f%%\n",    ann_ret    * 100))
cat(sprintf("Sharpe ratio         :  %.3f\n",       ann_sharpe))
cat(sprintf("CAPM alfa (ugentlig) : %+.4f%% (p = %.4f)\n",
            nw_capm[1, 1] * 100, nw_capm[1, 4]))
cat(sprintf("Max drawdown         : %+.2f%%\n",    max_dd     * 100))
cat("\n=== STEP 4 FÆRDIG — ANALYSE KOMPLET ===\n")

# ==============================================================================
# STEP 4.7: Tabel 3 — Betting Against Beta: Japan (akademisk stil)
# ==============================================================================
# Rækker:
#   - Excess Return    : månedligt gennemsnitligt excess return (%)
#   - CAPM Alpha       : månedlig CAPM-alpha (%)
#   - FF3 Alpha        : månedlig Fama-French 3-faktor alpha (%)
#   - FF4 Alpha        : månedlig Carhart 4-faktor alpha (%)
#   - Beta (ex ante)   : gennemsnitlig pre-formation beta_lag
#   - Beta (realized)  : realiseret CAPM-beta
#   - Volatility       : annualiseret volatilitet, baseret på månedlige afkast (%)
#   - Sharpe Ratio     : annualiseret Sharpe-ratio
#   - Antal aktier     : gennemsnitligt antal aktier pr. portefølje
#
# Output:
#   - BAB_tabel3.png
#   - table3_japan.csv
# ==============================================================================

cat("\n=== STEP 4.7: Tabel 3 — Betting Against Beta: Japan ===\n")

library(data.table)
library(ggplot2)
library(sandwich)
library(lmtest)

# ------------------------------------------------------------------------------
# 0) Sikkerhed: sørg for at centrale objekter er data.tables
# ------------------------------------------------------------------------------

setDT(bab_clean)
setDT(data_bab)

# Arbejd på en kopi, så vi ikke utilsigtet ændrer data_bab brugt senere
data_bab_tbl <- copy(data_bab)

# ------------------------------------------------------------------------------
# 1) Konstruér månedlige BAB-afkast
# ------------------------------------------------------------------------------

bab_clean[, `:=`(
  yr  = as.integer(format(date, "%Y")),
  mth = as.integer(format(date, "%m"))
)]

bab_mth <- bab_clean[!is.na(r_bab), .(
  r_bab_m = prod(1 + r_bab, na.rm = TRUE) - 1,
  r_mkt_m = prod(1 + r_mkt, na.rm = TRUE) - 1
), by = .(yr, mth)]

setorder(bab_mth, yr, mth)

# ------------------------------------------------------------------------------
# 2) Tildel deciler efter beta_lag hver uge
# ------------------------------------------------------------------------------

data_bab_tbl <- data_bab_tbl[
  !is.na(beta_lag) &
    !is.na(ret_exc_wk) &
    !is.na(mkt_exc_wk) &
    is.finite(beta_lag) &
    is.finite(ret_exc_wk) &
    is.finite(mkt_exc_wk)
]

setorder(data_bab_tbl, iso_year, iso_week, beta_lag)

# Robust ntile/decil-opdeling:
# Decil 1 = laveste beta
# Decil 10 = højeste beta
data_bab_tbl[, n_week := .N, by = .(iso_year, iso_week)]

data_bab_tbl[, beta_rank := frank(
  beta_lag,
  ties.method = "first"
), by = .(iso_year, iso_week)]

data_bab_tbl[, decil := pmin(
  10L,
  floor((beta_rank - 1) * 10 / n_week) + 1L
)]

# ------------------------------------------------------------------------------
# 3) Ugentlige decilafkast og gennemsnitligt antal aktier
# ------------------------------------------------------------------------------

decil_wk <- data_bab_tbl[!is.na(decil), .(
  ret_wk = mean(ret_exc_wk, na.rm = TRUE),
  n_stk  = .N
), by = .(decil, iso_year, iso_week)]

# Gennemsnitligt antal aktier pr. decil pr. uge
n_decil_summary <- decil_wk[, .(
  n_avg = mean(n_stk, na.rm = TRUE)
), by = decil][order(decil)]

# Gennemsnitligt antal aktier i BAB-benene
n_bab_summary <- bab_clean[, .(
  n_avg      = mean(n_low + n_high, na.rm = TRUE),
  n_low_avg  = mean(n_low, na.rm = TRUE),
  n_high_avg = mean(n_high, na.rm = TRUE)
)]

cat("\nGennemsnitligt antal aktier pr. decil:\n")
print(n_decil_summary)

cat("\nGennemsnitligt antal aktier i BAB-ben:\n")
print(n_bab_summary)

# ------------------------------------------------------------------------------
# 4) Aggreger decilafkast fra ugentligt til månedligt
# ------------------------------------------------------------------------------

wk_dates <- data_bab_tbl[, .(
  date = max(date, na.rm = TRUE)
), by = .(iso_year, iso_week)]

wk_dates[, `:=`(
  yr  = as.integer(format(date, "%Y")),
  mth = as.integer(format(date, "%m"))
)]

decil_wk <- merge(
  decil_wk,
  wk_dates[, .(iso_year, iso_week, yr, mth)],
  by = c("iso_year", "iso_week"),
  all.x = TRUE
)

decil_mth <- decil_wk[, .(
  ret_mth = prod(1 + ret_wk, na.rm = TRUE) - 1,
  n_wk    = .N
), by = .(decil, yr, mth)]

setorder(decil_mth, decil, yr, mth)

# Tilføj markedets månedlige excess return
mkt_mth <- bab_mth[, .(yr, mth, r_mkt_m)]

decil_mth <- merge(
  decil_mth,
  mkt_mth,
  by = c("yr", "mth"),
  all.x = TRUE
)

# Sanity check
rng <- range(decil_mth$ret_mth, na.rm = TRUE)

cat(sprintf(
  "\nDecil-månedsafkast range: [%.4f, %.4f]\n",
  rng[1], rng[2]
))

if (any(abs(decil_mth$ret_mth) > 2, na.rm = TRUE)) {
  warning("Decil-månedsafkast over 200%. Kontrollér data og outliers.")
}

# ------------------------------------------------------------------------------
# 5) Hent Kenneth French Japan-faktorer
# ------------------------------------------------------------------------------

download_ff_zip <- function(url, local_fallback = NULL) {
  
  zip_path <- tempfile(fileext = ".zip")
  
  # Forsøg 1: curl
  if (requireNamespace("curl", quietly = TRUE)) {
    ok <- tryCatch({
      h <- curl::new_handle()
      curl::handle_setopt(
        h,
        followlocation = TRUE,
        ssl_verifypeer = 0L,
        useragent = "Mozilla/5.0 (R)"
      )
      curl::curl_download(url, zip_path, handle = h, quiet = TRUE)
      file.exists(zip_path) && file.info(zip_path)$size > 1000
    }, error = function(e) FALSE)
    
    if (isTRUE(ok)) return(zip_path)
  }
  
  # Forsøg 2: base R download.file
  old_to <- getOption("timeout")
  options(timeout = 300)
  
  ok <- tryCatch({
    suppressWarnings(
      utils::download.file(
        url,
        zip_path,
        quiet = TRUE,
        mode = "wb",
        method = "libcurl"
      )
    )
    file.exists(zip_path) && file.info(zip_path)$size > 1000
  }, error = function(e) FALSE)
  
  options(timeout = old_to)
  
  if (isTRUE(ok)) return(zip_path)
  
  # Forsøg 3: lokal fallback
  if (!is.null(local_fallback) && file.exists(local_fallback)) {
    return(local_fallback)
  }
  
  NULL
}

parse_ff_zip <- function(url, col_names, label, local_fallback = NULL) {
  
  out <- tryCatch({
    
    # Hvis lokal CSV findes, bruges den direkte
    if (!is.null(local_fallback) &&
        file.exists(local_fallback) &&
        grepl("\\.csv$|\\.CSV$", local_fallback)) {
      
      csv_path <- local_fallback
      
    } else {
      
      zip_path <- download_ff_zip(url, local_fallback = local_fallback)
      
      if (is.null(zip_path)) {
        stop(
          "kunne ikke hente filen. Tjek internet eller placér ",
          basename(url),
          " som .zip eller .csv i arbejdsmappen."
        )
      }
      
      if (grepl("\\.csv$|\\.CSV$", zip_path)) {
        csv_path <- zip_path
      } else {
        tmp_dir <- tempfile("ff_")
        dir.create(tmp_dir)
        
        files <- tryCatch(
          unzip(zip_path, exdir = tmp_dir),
          warning = function(w) character(0),
          error   = function(e) character(0)
        )
        
        csv_path <- files[grepl("\\.csv$|\\.CSV$", files)][1]
        
        if (length(csv_path) == 0 || is.na(csv_path)) {
          stop("ingen gyldig CSV i zip-filen.")
        }
      }
    }
    
    raw <- readLines(csv_path, warn = FALSE)
    
    # Kenneth French-filer har datalinjer, der starter med YYYYMM
    start_i <- which(grepl("^\\s*[0-9]{6}\\s*,", raw))[1]
    
    if (is.na(start_i)) {
      stop("fandt ikke datastart i CSV.")
    }
    
    end_i <- which(grepl("^\\s*$|Annual Factors|^\\s*--", raw))
    end_i <- end_i[end_i > start_i][1] - 1L
    
    if (is.na(end_i)) {
      end_i <- length(raw)
    }
    
    dt <- fread(
      text = raw[start_i:end_i],
      header = FALSE,
      col.names = col_names,
      strip.white = TRUE
    )
    
    dt <- dt[
      grepl("^[0-9]{6}$", gsub("\\s", "", as.character(yyyymm)))
    ]
    
    for (nm in setdiff(col_names, "yyyymm")) {
      dt[, (nm) := as.numeric(get(nm)) / 100]
    }
    
    dt[, `:=`(
      yr  = as.integer(yyyymm) %/% 100L,
      mth = as.integer(yyyymm) %% 100L
    )]
    
    dt[, yyyymm := NULL]
    
    cat(sprintf("  %s indlæst: %d måneder\n", label, nrow(dt)))
    
    dt
    
  }, error = function(e) {
    cat(sprintf("  %s fejlede: %s\n", label, conditionMessage(e)))
    NULL
  })
  
  out
}

find_local_ff <- function(candidates) {
  for (nm in candidates) {
    if (file.exists(nm)) return(nm)
  }
  NULL
}

cat("\nHenter Kenneth French Japan-faktorer...\n")

ff3_url <- "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/Japan_3_Factors_CSV.zip"
mom_url <- "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/Japan_Mom_Factor_CSV.zip"

ff3_local <- find_local_ff(c(
  "Japan_3_Factors_CSV.zip",
  "Japan_3_Factors.CSV",
  "Japan_3_Factors.csv",
  "Japan_3_Factors_CSV.csv"
))

mom_local <- find_local_ff(c(
  "Japan_Mom_Factor_CSV.zip",
  "Japan_Mom_Factor.CSV",
  "Japan_Mom_Factor.csv",
  "Japan_MOM_Factor.csv",
  "Japan_Mom_Factor_CSV.csv"
))

if (!is.null(ff3_local)) {
  cat("  lokal FF3 fundet: ", ff3_local, "\n", sep = "")
}

if (!is.null(mom_local)) {
  cat("  lokal MOM fundet: ", mom_local, "\n", sep = "")
}

ff3_dt <- parse_ff_zip(
  ff3_url,
  c("yyyymm", "mkt_rf", "smb", "hml", "rf"),
  "FF3 Japan",
  local_fallback = ff3_local
)

mom_dt <- parse_ff_zip(
  mom_url,
  c("yyyymm", "mom"),
  "MOM Japan",
  local_fallback = mom_local
)

ff_data <- NULL
ff_ok   <- FALSE

if (!is.null(ff3_dt) && !is.null(mom_dt)) {
  
  ff_data <- merge(
    ff3_dt[, .(yr, mth, mkt_rf, smb, hml, rf)],
    mom_dt[, .(yr, mth, mom)],
    by = c("yr", "mth"),
    all = FALSE
  )
  
  ff_ok <- nrow(ff_data) > 100L
  
  if (ff_ok) {
    ff_min <- ff_data[which.min(yr * 100 + mth)]
    ff_max <- ff_data[which.max(yr * 100 + mth)]
    
    cat(sprintf(
      "FF-paneldata: %d måneder (%d-%02d til %d-%02d)\n",
      nrow(ff_data),
      ff_min$yr, ff_min$mth,
      ff_max$yr, ff_max$mth
    ))
  }
}

# ------------------------------------------------------------------------------
# 6) Hjælpefunktion: månedlige porteføljestatistikker
# ------------------------------------------------------------------------------

calc_stats_mth <- function(ret_m, mkt_m, yr_v, mth_v) {
  
  df <- data.table(
    ret = ret_m,
    mkt = mkt_m,
    yr  = yr_v,
    mth = mth_v
  )
  
  df <- df[
    !is.na(ret) &
      !is.na(mkt) &
      is.finite(ret) &
      is.finite(mkt)
  ]
  
  if (nrow(df) < 12) {
    return(list(
      exc = NA_real_, exc_t = NA_real_,
      alp = NA_real_, alp_t = NA_real_,
      ff3_a = NA_real_, ff3_t = NA_real_,
      ff4_a = NA_real_, ff4_t = NA_real_,
      b_r = NA_real_,
      vol = NA_real_,
      sharpe = NA_real_
    ))
  }
  
  # Excess return
  m1 <- lm(ret ~ 1, data = df)
  
  nw1 <- coeftest(
    m1,
    vcov = NeweyWest(m1, lag = 3, prewhite = FALSE)
  )
  
  exc   <- nw1[1, 1] * 100
  exc_t <- nw1[1, 3]
  
  # CAPM alpha
  m2 <- lm(ret ~ mkt, data = df)
  
  nw2 <- coeftest(
    m2,
    vcov = NeweyWest(m2, lag = 3, prewhite = FALSE)
  )
  
  alp   <- nw2[1, 1] * 100
  alp_t <- nw2[1, 3]
  b_r   <- coef(m2)[["mkt"]]
  
  # FF3 og FF4 alpha
  ff3_a <- NA_real_
  ff3_t <- NA_real_
  ff4_a <- NA_real_
  ff4_t <- NA_real_
  
  if (isTRUE(ff_ok) && !is.null(ff_data)) {
    
    dff <- merge(
      df,
      ff_data[, .(yr, mth, mkt_rf, smb, hml, mom)],
      by = c("yr", "mth"),
      all = FALSE
    )
    
    dff <- dff[complete.cases(dff)]
    
    if (nrow(dff) >= 24) {
      
      m3 <- lm(ret ~ mkt_rf + smb + hml, data = dff)
      
      nw3 <- coeftest(
        m3,
        vcov = NeweyWest(m3, lag = 3, prewhite = FALSE)
      )
      
      ff3_a <- nw3[1, 1] * 100
      ff3_t <- nw3[1, 3]
      
      m4 <- lm(ret ~ mkt_rf + smb + hml + mom, data = dff)
      
      nw4 <- coeftest(
        m4,
        vcov = NeweyWest(m4, lag = 3, prewhite = FALSE)
      )
      
      ff4_a <- nw4[1, 1] * 100
      ff4_t <- nw4[1, 3]
    }
  }
  
  # Annualiseret volatilitet og Sharpe-ratio
  vol <- sd(df$ret, na.rm = TRUE) * sqrt(12) * 100
  
  sharpe <- mean(df$ret, na.rm = TRUE) /
    sd(df$ret, na.rm = TRUE) * sqrt(12)
  
  list(
    exc = exc,
    exc_t = exc_t,
    alp = alp,
    alp_t = alp_t,
    ff3_a = ff3_a,
    ff3_t = ff3_t,
    ff4_a = ff4_a,
    ff4_t = ff4_t,
    b_r = b_r,
    vol = vol,
    sharpe = sharpe
  )
}

# ------------------------------------------------------------------------------
# 7) Beregn statistik for P1-P10 og BAB
# ------------------------------------------------------------------------------

ante_avg <- data_bab_tbl[!is.na(decil) & !is.na(beta_lag), .(
  b_ante = mean(beta_lag, na.rm = TRUE)
), by = decil][order(decil)]

n_port <- 11

stat_keys <- c(
  "exc", "exc_t",
  "alp", "alp_t",
  "ff3_a", "ff3_t",
  "ff4_a", "ff4_t",
  "b_r", "b_ante",
  "vol", "sharpe",
  "n_avg"
)

stats <- setNames(
  lapply(stat_keys, function(k) rep(NA_real_, n_port)),
  stat_keys
)

# P1-P10
for (d in 1:10) {
  
  sub <- decil_mth[decil == d]
  
  st <- calc_stats_mth(
    sub$ret_mth,
    sub$r_mkt_m,
    sub$yr,
    sub$mth
  )
  
  st$b_ante <- ante_avg[decil == d, b_ante]
  st$n_avg  <- n_decil_summary[decil == d, n_avg]
  
  for (k in stat_keys) {
    v <- st[[k]]
    stats[[k]][d] <- suppressWarnings(as.numeric(unname(v))[1])
  }
}

# BAB
st_bab <- calc_stats_mth(
  bab_mth$r_bab_m,
  bab_mth$r_mkt_m,
  bab_mth$yr,
  bab_mth$mth
)

# BAB har ex ante-beta tæt på 0 ved konstruktion
st_bab$b_ante <- 0

# For BAB angives samlet gennemsnitligt antal aktier i low- og high-beta-benene
st_bab$n_avg <- n_bab_summary$n_avg

for (k in stat_keys) {
  v <- st_bab[[k]]
  stats[[k]][11] <- suppressWarnings(as.numeric(unname(v))[1])
}

port_lbl <- c(paste0("P", 1:10), "BAB")

has_ff <- isTRUE(ff_ok) &&
  !all(is.na(stats$ff3_a)) &&
  !all(is.na(stats$ff4_a))

cat("\nTabel 3 — månedlige tal:\n")
for (k in stat_keys) {
  cat(sprintf(
    "  %-8s : %s\n",
    k,
    paste(
      formatC(stats[[k]], format = "f", digits = 3, width = 7),
      collapse = " "
    )
  ))
}

# Eksportdata
tbl <- data.table(portfolio = port_lbl)

for (k in stat_keys) {
  tbl[, (k) := stats[[k]]]
}

# ------------------------------------------------------------------------------
# 8) Byg long-format til grafisk tabel
# ------------------------------------------------------------------------------

rdef <- list(
  list(lbl = "Excess return", val = "exc", tv = "exc_t", t = TRUE),
  list(lbl = "CAPM-alpha",    val = "alp", tv = "alp_t", t = TRUE)
)

if (has_ff) {
  rdef[[length(rdef) + 1]] <- list(
    lbl = "FF3-alpha",
    val = "ff3_a",
    tv = "ff3_t",
    t = TRUE
  )
  
  rdef[[length(rdef) + 1]] <- list(
    lbl = "FF4-alpha",
    val = "ff4_a",
    tv = "ff4_t",
    t = TRUE
  )
}

rdef[[length(rdef) + 1]] <- list(
  lbl = "Beta (ex ante)",
  val = "b_ante",
  tv = NA,
  t = FALSE
)

rdef[[length(rdef) + 1]] <- list(
  lbl = "Beta (realized)",
  val = "b_r",
  tv = NA,
  t = FALSE
)

rdef[[length(rdef) + 1]] <- list(
  lbl = "Volatility",
  val = "vol",
  tv = NA,
  t = FALSE
)

rdef[[length(rdef) + 1]] <- list(
  lbl = "Sharpe ratio",
  val = "sharpe",
  tv = NA,
  t = FALSE
)

rdef[[length(rdef) + 1]] <- list(
  lbl = "Antal aktier, gns.",
  val = "n_avg",
  tv = NA,
  t = FALSE
)

# ------------------------------------------------------------------------------
# 9) Formateringsfunktioner
# ------------------------------------------------------------------------------

fmt_num <- function(x, d = 2) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0) return("\u2014")
  x <- x[1]
  
  if (is.na(x)) {
    "\u2014"
  } else {
    formatC(x, format = "f", digits = d, big.mark = "")
  }
}

sig_star <- function(t) {
  t <- suppressWarnings(as.numeric(t))
  if (length(t) == 0 || is.na(t)) return("")
  
  if (abs(t) >= 1.96) {
    "**"
  } else if (abs(t) >= 1.645) {
    "*"
  } else {
    ""
  }
}

# ------------------------------------------------------------------------------
# 10) Pre-allokér celler til ggplot-tabellen
# ------------------------------------------------------------------------------

n_cells <- sum(vapply(
  rdef,
  function(r) if (isTRUE(r$t)) 2L else 1L,
  integer(1)
)) * n_port

gg_group_id <- integer(n_cells)
gg_cn       <- integer(n_cells)
gg_lbl      <- character(n_cells)
gg_txt      <- character(n_cells)
gg_bold     <- logical(n_cells)
gg_is_t     <- logical(n_cells)
gg_yp       <- numeric(n_cells)

idx <- 0L
cur_y <- 0

for (ri in seq_along(rdef)) {
  
  rd <- rdef[[ri]]
  
  if (ri > 1) {
    cur_y <- cur_y - 1.35
  }
  
  y_est <- cur_y
  
  vals <- as.numeric(stats[[rd$val]])
  
  tvs <- if (isTRUE(rd$t)) {
    as.numeric(stats[[rd$tv]])
  } else {
    rep(NA_real_, n_port)
  }
  
  if (length(vals) != n_port) vals <- rep(NA_real_, n_port)
  if (length(tvs)  != n_port) tvs  <- rep(NA_real_, n_port)
  
  # Estimat-række
  for (ci in 1:n_port) {
    
    idx <- idx + 1L
    
    val_sc <- vals[ci]
    tv_sc  <- tvs[ci]
    
    digits_i <- if (identical(rd$val, "n_avg")) 0 else 2
    
    txt_i <- fmt_num(val_sc, digits_i)
    
    # Asterisker tilføjes kun til rækker med t-statistik
    if (isTRUE(rd$t)) {
      txt_i <- paste0(txt_i, sig_star(tv_sc))
    }
    
    gg_group_id[idx] <- ri
    gg_cn[idx]       <- as.integer(ci)
    gg_lbl[idx]      <- as.character(rd$lbl)[1]
    gg_txt[idx]      <- txt_i
    gg_bold[idx]     <- isTRUE(rd$t) && !is.na(tv_sc) && abs(tv_sc) >= 1.96
    gg_is_t[idx]     <- FALSE
    gg_yp[idx]       <- y_est
  }
  
  # t-stat-række
  if (isTRUE(rd$t)) {
    
    cur_y <- cur_y - 0.55
    y_t <- cur_y
    
    for (ci in 1:n_port) {
      
      idx <- idx + 1L
      
      tv_sc <- tvs[ci]
      
      t_str <- if (is.na(tv_sc)) {
        ""
      } else {
        paste0("(", fmt_num(tv_sc, 2), ")")
      }
      
      gg_group_id[idx] <- ri
      gg_cn[idx]       <- as.integer(ci)
      gg_lbl[idx]      <- as.character(rd$lbl)[1]
      gg_txt[idx]      <- t_str
      gg_bold[idx]     <- !is.na(tv_sc) && abs(tv_sc) >= 1.96
      gg_is_t[idx]     <- TRUE
      gg_yp[idx]       <- y_t
    }
  }
}

long_df <- data.frame(
  group_id = gg_group_id,
  cn       = gg_cn,
  lbl      = gg_lbl,
  txt      = gg_txt,
  bold     = gg_bold,
  is_t     = gg_is_t,
  yp       = gg_yp,
  stringsAsFactors = FALSE
)

long_dt <- as.data.table(long_df)

lbl_y <- long_dt[is_t == FALSE, .(
  ym = mean(yp),
  lbl = unique(lbl)
), by = group_id]

setorder(lbl_y, group_id)

# ------------------------------------------------------------------------------
# 11) Tegn akademisk tabel
# ------------------------------------------------------------------------------

y_min <- min(long_dt$yp)
y_max <- max(long_dt$yp)

y_h      <- y_max + 1.35
y_h_sub  <- y_h - 0.55
y_h_rule <- y_h_sub - 0.40
y_top    <- y_h + 0.85
y_bot    <- y_min - 0.70

x_l       <- -2.55
x_r       <- 11.55
x_bab_l   <- 10.55
x_label_l <- -2.40

col_text     <- "#0A0A0A"
col_tstat    <- "#5A5A5A"
col_sub      <- "#666666"
col_rule     <- "#0A0A0A"
col_sub_rule <- "#BBBBBB"
col_bab_bg   <- "#F7F7F7"

hdr_dt <- data.table(
  cn  = 1:n_port,
  lbl = c(paste0("P", 1:10), "BAB"),
  sub = c("Lav \u03b2", rep("", 8), "H\u00f8j \u03b2", "")
)

grp_sep <- long_dt[, .(
  gmin = min(yp)
), by = group_id]

setorder(grp_sep, group_id)

grp_sep[, ysep := gmin - 0.58]

grp_sep <- grp_sep[group_id < max(grp_sep$group_id)]

mnd_dk <- c(
  "jan", "feb", "mar", "apr", "maj", "jun",
  "jul", "aug", "sep", "okt", "nov", "dec"
)

sample_start <- bab_mth[which.min(yr * 100 + mth)]
sample_end   <- bab_mth[which.max(yr * 100 + mth)]

sample_txt <- sprintf(
  "%s %d til %s %d",
  mnd_dk[sample_start$mth],
  sample_start$yr,
  mnd_dk[sample_end$mth],
  sample_end$yr
)

footnote <- paste0(
  "Newey\u2013West t-statistik i parentes (lag = 3 m\u00e5neder). ",
  "** angiver signifikans p\u00e5 5%-niveau; * angiver signifikans p\u00e5 10%-niveau. ",
  "Antal aktier angiver gennemsnitligt antal aktier pr. uge i hver decil; ",
  "for BAB angives gennemsnitligt samlet antal aktier i low- og high-beta-benene. ",
  "Volatilitet og Sharpe-ratio er annualiserede. ",
  "Faktordata: Kenneth French."
)

p_tbl <- ggplot() +
  
  # BAB-kolonne
  geom_rect(
    aes(xmin = x_bab_l, xmax = x_r, ymin = y_bot, ymax = y_top),
    fill = col_bab_bg,
    color = NA
  ) +
  
  # Top- og bundregler
  geom_segment(
    aes(x = x_l, xend = x_r, y = y_top, yend = y_top),
    linewidth = 1.0,
    color = col_rule
  ) +
  geom_segment(
    aes(x = x_l, xend = x_r, y = y_bot, yend = y_bot),
    linewidth = 1.0,
    color = col_rule
  ) +
  
  # Linje under header
  geom_segment(
    aes(x = x_l, xend = x_r, y = y_h_rule, yend = y_h_rule),
    linewidth = 0.5,
    color = col_rule
  ) +
  
  # Separatorer mellem rækkegrupper
  geom_segment(
    data = grp_sep,
    aes(x = x_l, xend = x_r, y = ysep, yend = ysep),
    linewidth = 0.25,
    color = col_sub_rule
  ) +
  
  # Header
  geom_text(
    data = hdr_dt,
    aes(x = cn, y = y_h, label = lbl),
    fontface = "bold",
    size = 3.4,
    color = col_text
  ) +
  
  geom_text(
    data = hdr_dt[sub != ""],
    aes(x = cn, y = y_h_sub, label = sub),
    fontface = "italic",
    size = 2.8,
    color = col_sub
  ) +
  
  # Rækkenavne
  geom_text(
    data = lbl_y,
    aes(x = x_label_l, y = ym, label = lbl),
    hjust = 0,
    size = 3.15,
    color = col_text
  ) +
  
  # Estimater, normal
  geom_text(
    data = long_dt[is_t == FALSE & bold == FALSE],
    aes(x = cn, y = yp, label = txt),
    size = 3.05,
    color = col_text
  ) +
  
  # Estimater, fed
  geom_text(
    data = long_dt[is_t == FALSE & bold == TRUE],
    aes(x = cn, y = yp, label = txt),
    size = 3.05,
    color = col_text,
    fontface = "bold"
  ) +
  
  # t-statistikker, normal
  geom_text(
    data = long_dt[is_t == TRUE & bold == FALSE],
    aes(x = cn, y = yp, label = txt),
    size = 2.70,
    color = col_tstat
  ) +
  
  # t-statistikker, fed
  geom_text(
    data = long_dt[is_t == TRUE & bold == TRUE],
    aes(x = cn, y = yp, label = txt),
    size = 2.70,
    color = col_text,
    fontface = "bold"
  ) +
  
  scale_x_continuous(
    limits = c(x_l - 0.15, x_r + 0.15),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(y_bot - 0.35, y_top + 1.10),
    expand = c(0, 0)
  ) +
  
  theme_void(base_family = "") +
  
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(22, 18, 16, 18),
    plot.title = element_text(
      face = "bold",
      size = 12.5,
      hjust = 0,
      color = col_text,
      margin = margin(b = 3)
    ),
    plot.subtitle = element_text(
      size = 9.0,
      hjust = 0,
      color = "#444444",
      margin = margin(b = 12)
    ),
    plot.caption = element_text(
      size = 7.6,
      hjust = 0,
      color = "#666666",
      margin = margin(t = 12)
    )
  ) +
  
  labs(
    title = "Betting Against Beta \u2014 Japan",
    subtitle = paste0(
      "M\u00e5nedlige excess returns og alphaer i procent. ",
      "Sampleperiode: ", sample_txt, ". ",
      "Porteføljer P1-P10 er sorteret efter ex ante-beta; ",
      "P1 har lavest beta og P10 har h\u00f8jest beta."
    ),
    caption = footnote
  )

# Vis tabel i RStudio
print(p_tbl)

# ------------------------------------------------------------------------------
# 12) Gem output
# ------------------------------------------------------------------------------

n_row_grps <- length(rdef)
fig_h <- 4.0 + n_row_grps * 0.45

ggsave(
  filename = "BAB_tabel3.png",
  plot     = p_tbl,
  width    = 14.5,
  height   = max(fig_h, 6.8),
  dpi      = 300
)

cat("Gemt: BAB_tabel3.png\n")

fwrite(tbl, "table3_japan.csv")

cat("Gemt: table3_japan.csv\n")

# Gem også de centrale hjælpeobjekter til evt. bilag
fwrite(n_decil_summary, "table3_avg_number_stocks_decils.csv")
fwrite(n_bab_summary,   "table3_avg_number_stocks_bab.csv")

cat("Gemt: table3_avg_number_stocks_decils.csv\n")
cat("Gemt: table3_avg_number_stocks_bab.csv\n")

cat("\n=== STEP 4.7 FÆRDIG ===\n")


# ==============================================================================
# Figur: Kumulativt afkast for BAB, lav-beta og høj-beta
# ==============================================================================

bab_clean[, `:=`(
  cum_bab  = cumprod(1 + r_bab)  - 1,
  cum_low  = cumprod(1 + r_low)  - 1,
  cum_high = cumprod(1 + r_high) - 1
)]

plot_data <- melt(
  bab_clean[, .(date, cum_bab, cum_low, cum_high)],
  id.vars = "date",
  variable.name = "Portfolio",
  value.name = "Cumulative_return"
)

plot_data[, Portfolio := fcase(
  Portfolio == "cum_bab",  "BAB",
  Portfolio == "cum_low",  "Low beta",
  Portfolio == "cum_high", "High beta"
)]

p_bab_cum <- ggplot(
  plot_data,
  aes(x = date, y = Cumulative_return * 100,
      color = Portfolio, linetype = Portfolio)
) +
  geom_line(linewidth = 0.75) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title    = "Kumulativt afkast for BAB-strategien i Japan",
    subtitle = "Ugentlige excess returns, 1991-2022",
    x        = NULL,
    y        = "Kumulativt afkast (%)",
    color    = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "BAB_kumuleret_afkast_hovedanalyse.png",
  plot = p_bab_cum,
  width = 12,
  height = 6,
  dpi = 300
)
