# Response-function regressions (proposal slides 16-17).
#   CROP  : ln(yield) ~ GDD + Heat + Frost + Precip + Precip^2 | NUTS3^crop + year   (per crop)
#   ENERGY: ln(energy per capita) ~ HDD (+ CDD) + ln(price)    | country + year       (per fuel)
# Reads 9.crop_panel_nuts3_estimation.csv, 16.energy_panel_estimation.csv.
# Writes coefficient tables to Amoc/results/. Author: Marco Bova
suppressMessages({library(fixest); library(data.table)})
d   <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
res <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/results")
dir.create(res, showWarnings = FALSE, recursive = TRUE)
save_tab <- function(x, f) writeLines(capture.output(print(x)), file.path(res, f))

# Hold the estimation sample constant across the columns of a build-up table. A fixed-effects
# spec drops not only rows with missing values (as OLS does) but also FIXED-EFFECT SINGLETONS -
# a region or year observed only once, which a dummy absorbs perfectly and which would inflate
# the reported fit. Its sample is therefore a strict subset of the OLS sample, so a raw OLS-vs-FE
# table compares columns estimated on slightly different rows. common_obs() returns the rows
# retained by ALL the fixed-effect specs in a table (obs() = the observations fixest actually
# used), iterated to a fixed point in case dropping one group orphans another. Refitting every
# column - including OLS - on this common set gives one shared N, so differences across columns
# reflect ONLY the fixed-effect structure, not a moving sample.
common_obs <- function(dt, fe_forms, cl) {
  repeat {
    kept <- Reduce(intersect, lapply(fe_forms, function(f) obs(feols(f, dt, cluster = cl))))
    if (length(kept) == nrow(dt)) return(dt)
    dt <- dt[kept]
  }
}

# Restricted wild-cluster bootstrap (Rademacher weights, null imposed; Cameron-Gelbach-Miller
# 2008, Roodman et al. 2019). Few-cluster inference for the benchmark: with only 15 (crop) / 29
# (energy) country clusters, cluster-robust asymptotics are unreliable, and CR2/Satterthwaite is
# infeasible on the ~800-region trend design. This partials the fixed effects out first (FWL) so
# each bootstrap draw is matrix algebra, not a model refit. Returns benchmark coefs, CRV1
# cluster-robust SE/t, and the bootstrap p-value per term.
wcb <- function(dep, regs, fe, cl, data, B = 999L, seed = 1L) {
  set.seed(seed)
  rr <- function(v) feols(as.formula(sprintf("%s ~ 1 | %s", v, fe)), data, notes = FALSE)
  y  <- resid(rr(dep)); X <- sapply(regs, function(v) resid(rr(v)))
  g  <- as.integer(factor(data[[cl]][obs(rr(dep))]))       # cluster id on the rows fixest kept
  gl <- split(seq_along(g), g); k <- ncol(X)
  fit <- function(Y) { b <- qr.solve(crossprod(X), crossprod(X, Y)); u <- Y - X %*% b
    meat <- Reduce(`+`, lapply(gl, function(i){ s <- crossprod(X[i, , drop = FALSE], u[i]); tcrossprod(s) }))
    Bi <- chol2inv(chol(crossprod(X))); V <- Bi %*% meat %*% Bi; list(b = b, se = sqrt(diag(V))) }
  m0 <- fit(y); t_obs <- m0$b / m0$se; p <- numeric(k)
  for (j in seq_len(k)) {                                  # impose H0: b_j = 0
    Xr <- X[, -j, drop = FALSE]; yhat <- Xr %*% qr.solve(crossprod(Xr), crossprod(Xr, y)); ur <- y - yhat
    cnt <- 0
    for (b in seq_len(B)) { w <- sample(c(-1, 1), length(gl), TRUE)[g]; mb <- fit(yhat + w * ur)
      cnt <- cnt + (abs(mb$b[j] / mb$se[j]) >= abs(t_obs[j])) }
    p[j] <- (cnt + 1) / (B + 1) }
  data.frame(term = regs, coef = as.numeric(m0$b), crv_se = m0$se, t = as.numeric(t_obs), wcb_p = p)
}

## ============= CROP: build-up to the Blanc & Schlenker (2017) benchmark =====
# Per crop, a build-up: pooled OLS -> two-way FE -> +region-specific quadratic time trends
# (the paper's benchmark: county FE + year FE + county-specific quadratic trends, which absorb
# region-specific technology yield growth). Mixed sample shown as the coverage robustness.
#   9  = NUTS3-only (15 countries)   10 = mixed-level (26; NUTS3 where CropStatHarm has yield,
#   else NUTS2/NUTS0). One crop set only (totals and components are the same tonnage twice).
crop_set <- c("Soft wheat", "Durum wheat", "Spring barley", "Winter barley")  # or c("Total wheat","Total barley")
load_crop <- function(f) { x <- fread(file.path(d, f))[crop %in% crop_set]
  x[, `:=`(precip_mj2 = precip_mj^2, t = year - 2005L, t2 = (year - 2005L)^2)]; x[] }  # rainfall^2, centred trend
cp9  <- load_crop("9.crop_panel_nuts3_estimation.csv")
cp10 <- load_crop("10.crop_panel_mixed_estimation.csv")
f_ols <- ln_yield ~ gdd_mj + heat_mj + frost_mj + precip_mj + precip_mj2                         # no FE
f_fe  <- ln_yield ~ gdd_mj + heat_mj + frost_mj + precip_mj + precip_mj2 | NUTS_ID + year         # two-way FE
f_tr  <- ln_yield ~ gdd_mj + heat_mj + frost_mj + precip_mj + precip_mj2 | NUTS_ID[t, t2] + year  # + region trends

for (k in crop_set) {              # OLS -> +FE -> +region trends (benchmark); Mixed = coverage robustness
  s9  <- common_obs(cp9[crop == k],  list(f_tr), ~cntr)   # anchor the sample to the most saturated spec
  s10 <- common_obs(cp10[crop == k], list(f_tr), ~cntr)
  tab <- etable(feols(f_ols, s9, cluster = ~cntr), feols(f_fe, s9, cluster = ~cntr),
                feols(f_tr,  s9, cluster = ~cntr), feols(f_tr, s10, cluster = ~cntr),
                headers = c("NUTS3 OLS", "NUTS3 FE", "NUTS3 FE+Trends", "Mixed FE+Trends"),
                fitstat = ~ n + r2 + wr2)
  cat("\n===== CROP:", k, "=====\n"); print(tab)
  save_tab(tab, paste0("crop_", gsub(" ", "_", k), "_n3_vs_mixed.txt"))
}

## --- CROP inference robustness (caveats: few clusters #5, snapped regions #1) ---
# The benchmark is estimated on 15 country clusters, where cluster-robust asymptotics are
# unreliable. Three checks on the benchmark spec: (a) two-way clustering (country + year) for
# common time shocks; (b) a restricted wild-cluster bootstrap (wcb(), Rademacher, null imposed) -
# the standard few-cluster inference; (c) dropping the snapped Ionian NUTS3 (EL62*, <=40 km E-OBS
# grid snaps) as a measurement-error defense - expected to move nothing (3-4 tiny islands),
# reported as defense, not a live threat.
wx_crop <- c("gdd_mj", "heat_mj", "frost_mj", "precip_mj", "precip_mj2")
for (k in crop_set) {
  s9   <- common_obs(cp9[crop == k], list(f_tr), ~cntr)
  m_bm <- feols(f_tr, s9, cluster = ~cntr)                    # benchmark, one-way cluster (main)
  m_2w <- feols(f_tr, s9, cluster = ~cntr + year)            # + two-way cluster (country + year)
  s9d  <- s9[!grepl("^EL62", NUTS_ID)]                       # drop snapped Ionian NUTS3
  m_dr <- feols(f_tr, s9d, cluster = ~cntr)
  rtab <- etable(m_bm, m_2w, m_dr,
                 headers = c("Benchmark (1-way)", "Two-way cluster", "Drop snapped EL62*"),
                 fitstat = ~ n + r2)
  ct <- wcb("ln_yield", wx_crop, "NUTS_ID[t, t2] + year", "cntr", s9)
  cat("\n===== CROP robustness:", k, "=====\n"); print(rtab)
  cat("\n-- wild-cluster bootstrap p-values (999 reps, 15 country clusters) --\n"); print(ct)
  writeLines(c(capture.output(print(rtab)), "",
               "wild-cluster bootstrap (restricted, Rademacher, 999 reps, 15 country clusters):",
               capture.output(print(ct))),
             file.path(res, paste0("crop_", gsub(" ", "_", k), "_robustness.txt")))
}

## ======================= ENERGY ===========================================
ep <- fread(file.path(d, "16.energy_panel_estimation.csv"))
ep[, `:=`(t = year - 2007L, t2 = (year - 2007L)^2)]        # centred trend for country-specific slopes
# fuels run separately (different units + weather drivers). Outcome = per-capita (elasticity 1).
# Build-up to the SAME benchmark as crop (Blanc & Schlenker: unit FE + year FE + unit-specific
# quadratic trends, all three). Year FE alone absorbs the spatially-correlated common European
# weather and leaves HDD unidentified (even wrong-signed for gas); the benchmark recovers it.
energy_tab <- function(dt, wx, lab, file) {
  f_ols <- as.formula(paste("ln_energy_pc ~", wx, "+ ln_price"))                        # no FE
  f_yfe <- as.formula(paste("ln_energy_pc ~", wx, "+ ln_price | country_id + year"))     # + year FE
  f_ctr <- as.formula(paste("ln_energy_pc ~", wx, "+ ln_price | country_id[t, t2]"))     # + country trends
  f_bm  <- as.formula(paste("ln_energy_pc ~", wx, "+ ln_price | country_id[t, t2] + year"))  # BENCHMARK
  dt <- common_obs(dt, list(f_yfe, f_ctr, f_bm), ~country_id)   # constant sample across columns
  tab <- etable(feols(f_ols, dt, cluster = ~country_id), feols(f_yfe, dt, cluster = ~country_id),
                feols(f_ctr, dt, cluster = ~country_id), feols(f_bm, dt, cluster = ~country_id),
                headers = c("Pooled OLS", "Year FE", "Country trends", "Benchmark (main)"),
                fitstat = ~ n + r2 + wr2)
  cat("\n===== ENERGY:", lab, "=====\n"); print(tab); save_tab(tab, file)

  # inference robustness (caveats: few clusters #5, snapped region #1) on the benchmark spec
  m_bm2 <- feols(f_bm, dt, cluster = ~country_id)                 # benchmark (main)
  m_2w  <- feols(f_bm, dt, cluster = ~country_id + year)          # + two-way cluster
  m_dr  <- feols(f_bm, dt[country_id != "MLT"], cluster = ~country_id)  # drop Malta (snapped to Sicily)
  rtab  <- etable(m_bm2, m_2w, m_dr,
                  headers = c("Benchmark", "Two-way cluster", "Drop Malta (MLT)"),
                  fitstat = ~ n + r2)
  keep  <- c(trimws(strsplit(wx, "\\+")[[1]]), "ln_price")
  ct    <- wcb("ln_energy_pc", keep, "country_id[t, t2] + year", "country_id", dt)
  cat("\n===== ENERGY robustness:", lab, "=====\n"); print(rtab)
  cat("\n-- wild-cluster bootstrap p-values (999 reps, 29 country clusters) --\n"); print(ct)
  writeLines(c(capture.output(print(rtab)), "",
               "wild-cluster bootstrap (restricted, Rademacher, 999 reps, 29 country clusters):",
               capture.output(print(ct))),
             file.path(res, sub("\\.txt$", "_robustness.txt", file)))
}
energy_tab(ep[fuel == "Electricity"], "hdd_calendar + cdd_jja", "Electricity", "energy_electricity.txt")
energy_tab(ep[fuel == "Natural gas"], "hdd_octmar",             "Natural gas", "energy_gas.txt")
# extra robustness: free population elasticity (ln_energy on ln_population), year FE
save_tab(etable(feols(ln_energy ~ ln_population + hdd_calendar + cdd_jja + ln_price | country_id + year,
                      cluster = ~country_id, data = ep[fuel == "Electricity"])), "energy_pop_elasticity.txt")

cat("\nwrote coefficient tables ->", res, "\n")
