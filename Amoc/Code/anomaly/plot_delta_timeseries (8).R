# ============================================================
# TIME-SERIES overlays of the Europe-mean anomaly THROUGH the run,
# PLUS per-country-mean overlays.
#
# Per variable, two pages:
#   (1) EU-MEAN page (unchanged): all model runs on one axis, raw
#       monthly + LOESS + settled plateau + pooled control 95% band.
#   (2) COUNTRY page: small multiples, one panel per country, the
#       model runs overlaid (raw + LOESS + plateau dot). Country
#       cells are identified by point-in-polygon against the maps
#       world database and area-averaged (cos-lat). Native grids,
#       no interpolation.
#
# Control is DESEASONALISED (own 12-month climatology removed) then
# de-drifted, so the band/ratio are seasonal-free.
# ============================================================

library(ncdf4)
have_maps <- requireNamespace("maps", quietly = TRUE)

## rename-proof: source amoc_common.R from THIS script's own folder
## (keeps working if the folder is moved/renamed; run via `Rscript "<this file>"`).
.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)   # Rscript encodes spaces in --file= as ~+~
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: europe, MAX_PIC_FILES, EFFECT_K, publish_file,
## read_europe_cube, masked_mean_series, pic_dir_for, mean_null, and
## get_control_series (deseasonalised + de-drifted -> same ratio as the CSV).

## ---- options ----
## Countries to break out (names as in the maps 'world' database).
COUNTRIES <- c(
  "Austria","Belgium","Bulgaria","Croatia","Cyprus","Czech Republic","Denmark","Estonia",
  "Finland","France","Germany","Greece","Hungary","Ireland","Italy","Latvia","Lithuania",
  "Luxembourg","Netherlands","Poland","Portugal","Romania","Slovakia","Slovenia",
  "Spain","Sweden",                                  # EU-27
  "Switzerland","Norway","UK")                       # + non-EU requested

## ---- per-grid country label: cell_country() now lives in amoc_common.R ----

## ---- LOESS trajectory + plateau for one series ----
fit_traj <- function(ts) {
  n <- length(ts); tm <- (0:(n - 1)) / 12; pl <- max(12, round(n / 3))
  sp <- max(0.20, min(0.60, 300 / n))
  sm <- tryCatch(predict(loess(ts ~ tm, span = sp, degree = 1), tm), error = function(e) ts)
  list(tm = tm, ts = ts, sm = sm, pl = pl, plate = mean(tail(ts, pl), na.rm = TRUE), xend = max(tm))
}

## ---- runs + paths ----
runs <- list(
  list(model = "EC-Earth3",       proto = "g01", folder = "ECHearth3_anomaly", col = "#E69F00"),
  list(model = "EC-Earth3",       proto = "u03", folder = "ECHearth3_anomaly", col = "#D55E00"),
  list(model = "HadGEM3-GC31-LL", proto = "g01", folder = "LL_anomaly",        col = "#0072B2"),
  list(model = "HadGEM3-GC31-LL", proto = "u03", folder = "LL_anomaly",        col = "#003A6B"),
  list(model = "HadGEM3-GC31-MM", proto = "g01", folder = "MM_anomaly",        col = "#009E73"),
  list(model = "HadGEM3-GC31-MM", proto = "u03", folder = "MM_anomaly",        col = "#00674A"),
  list(model = "IPSL-CM6A-LR",    proto = "u03", folder = "IPSL_anomaly",      col = "#e7298a")   # u03 only (no g01)
)
base <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/anomaly_output"
anomaly_path <- function(r, vn) file.path(base, r$folder,
  if (r$folder %in% c("ECHearth3_anomaly", "IPSL_anomaly")) sprintf("%s_Amon_%s_hos-%s-hos_anomaly.nc", vn, r$model, r$proto)
  else sprintf("%s_anomaly_%s-hos_minus_piControl_1850-1949.nc", vn, r$proto))

vars   <- c("tas", "pr", "tasmax", "tasmin")
outdir <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# RUN
# ============================================================
if (!have_maps) cat("NOTE: 'maps' package not installed -> country pages skipped (EU pages still produced).\n")

CP <- data.frame()   # per-country plateau collector

pdf_final <- file.path(outdir, "delta_timeseries_overlay.pdf")
pdf_tmp   <- tempfile(fileext = ".pdf")     # render locally, publish to OneDrive at the end
pdf(pdf_tmp, width = 11, height = 8.5)

for (vn in vars) {
  is_pr      <- (vn == "pr")
  ylab       <- if (is_pr) "Anomaly (mm/day)" else "Anomaly (K = \u00B0C change)"
  ylab_short <- if (is_pr) "mm/day" else "K"

  ## ---- build per-run data: EU series (+ratio) and country series ----
  RD <- Filter(Negate(is.null), lapply(runs, function(r) {
    p <- anomaly_path(r, vn); if (!file.exists(p)) return(NULL)
    cb <- read_europe_cube(p); ncell <- dim(cb$v)[1] * dim(cb$v)[2]
    eu <- fit_traj(masked_mean_series(cb$v, cb$lat, seq_len(ncell)))
    ctrl_a <- { c <- get_control_series(p); c - mean(c) }
    eu$ratio  <- abs(eu$plate) / sd(mean_null(ctrl_a, eu$pl))
    eu$ctrl_a <- ctrl_a
    ctry <- NULL
    if (have_maps) {
      cc <- cell_country(basename(dirname(p)), cb$lon, cb$lat)
      ctry <- lapply(COUNTRIES, function(C) {
        sel <- which(cc == C); if (!length(sel)) return(NULL)
        ts <- masked_mean_series(cb$v, cb$lat, sel)
        if (!any(is.finite(ts))) NULL
        else { f <- fit_traj(ts); f$n_cells <- length(sel); f }  # n_cells -> CSV
      })
      names(ctry) <- COUNTRIES
    }
    list(r = r, eu = eu, country = ctry)
  }))
  if (!length(RD)) next

  if (have_maps) {                                   # report countries with no cells
    resolved <- vapply(COUNTRIES, function(C) any(vapply(RD, function(z) !is.null(z$country[[C]]), logical(1))), logical(1))
    if (any(!resolved)) cat(sprintf("  [%s] no cells (check name / too small for grid): %s\n",
                                    vn, paste(COUNTRIES[!resolved], collapse = ", ")))
    for (z in RD) for (C in COUNTRIES) {             # collect per-country plateaus
      f <- z$country[[C]]
      if (!is.null(f)) CP <- rbind(CP, data.frame(
        variable = vn, model = z$r$model, protocol = z$r$proto, country = C,
        plateau = round(f$plate, 3), n_cells = f$n_cells,
        run_years = round(length(f$ts) / 12), stringsAsFactors = FALSE))
    }
  }

  ## ---- PAGE 1: EU-mean overlay (unchanged) ----
  par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(4.5, 4.5, 3, 1), mgp = c(2.6, 0.8, 0))
  S <- lapply(RD, function(z) modifyList(z$eu, list(r = z$r)))
  band <- quantile(unlist(lapply(S, `[[`, "ctrl_a")), c(.025, .975), na.rm = TRUE)
  xmax <- max(vapply(S, `[[`, numeric(1), "xend"))
  yr   <- range(c(unlist(lapply(S, `[[`, "ts")), band, 0), na.rm = TRUE)
  plot(NA, xlim = c(0, xmax * 1.02), ylim = yr,
       xlab = "Years from start of hosing run (monthly data)", ylab = ylab,
       main = sprintf("Europe-mean %s through the hosing run \u2014 trajectory & plateau", vn))
  rect(0, band[1], xmax * 1.02, band[2], col = rgb(0.6, 0.6, 0.6, 0.28), border = NA)
  abline(h = 0, lty = 3)
  for (s in S) {
    lines(s$tm, s$ts, col = adjustcolor(s$r$col, alpha.f = 0.18))
    lines(s$tm, s$sm, col = s$r$col, lwd = 2.6)
    segments(s$xend - (s$pl - 1) / 12, s$plate, s$xend, s$plate, col = s$r$col, lwd = 5)
    points(s$xend, s$plate, pch = 19, col = s$r$col, cex = 1.1)
  }
  legend("bottomleft", bty = "n", cex = 0.78, lwd = 2.6,
         col = vapply(S, function(s) s$r$col, character(1)),
         legend = vapply(S, function(s) sprintf("%s %s   plateau %+.2f   ratio %.1f (%s)",
                         s$r$model, s$r$proto, s$plate, s$ratio,
                         if (s$ratio >= EFFECT_K) "effect" else "no"), character(1)))
  mtext("grey = control central-95% pooled across models (deseasonalised; each legend ratio uses its own model's control); thick bar = settled plateau",
        side = 3, line = 0.2, cex = 0.72, col = "grey30")

  ## ---- PAGE 2+: country small-multiples, paginated, SHARED y-scale ----
  if (have_maps) {
    allts <- unlist(lapply(RD, function(z) lapply(z$country, function(f) if (!is.null(f)) f$ts)))
    yrC   <- range(c(allts, 0), na.rm = TRUE)                          # wide raw-data range (as before)
    per_page <- 16; ncols <- 4; nrows <- 4
    chunks   <- split(COUNTRIES, ceiling(seq_along(COUNTRIES) / per_page))
    cols4 <- vapply(runs, function(r) r$col, character(1))
    labs4 <- vapply(runs, function(r) paste(sub("HadGEM3-GC31-", "", r$model), r$proto), character(1))
    at4   <- seq(0.12, 0.88, length.out = length(runs))
    for (chunk in chunks) {
      par(mfrow = c(nrows, ncols), oma = c(1, 0, 4, 0), mar = c(2.3, 2.8, 1.8, 0.6), mgp = c(1.6, 0.5, 0))
      for (ctry in chunk) {
        present <- Filter(function(z) !is.null(z$country[[ctry]]), RD)
        if (!length(present)) { plot.new(); box(); title(ctry, cex.main = 0.85); next }
        xm <- max(vapply(present, function(z) z$country[[ctry]]$xend, numeric(1)))
        plot(NA, xlim = c(0, xm), ylim = yrC, xlab = "", ylab = ylab_short,
             main = ctry, cex.main = 0.85, cex.axis = 0.7, cex.lab = 0.8)
        abline(h = 0, lty = 3, col = "grey60")
        for (z in present) {
          f <- z$country[[ctry]]
          lines(f$tm, f$ts, col = adjustcolor(z$r$col, alpha.f = 0.15))   # raw monthly
          lines(f$tm, f$sm, col = z$r$col, lwd = 1.8)                     # LOESS
          points(f$xend, f$plate, pch = 19, col = z$r$col, cex = 0.6)     # settled level
        }
      }
      mtext(sprintf("Country-mean %s through the hosing run  (shared scale; x = years; thick = LOESS, dot = plateau)", vn),
            outer = TRUE, line = 2.2, font = 2, cex = 0.95)
      for (k in seq_along(runs)) mtext(labs4[k], col = cols4[k], side = 3, line = 0.6, at = at4[k], outer = TRUE, cex = 0.72)
    }
  }
}
dev.off()
publish_file(pdf_tmp, pdf_final)

if (nrow(CP)) {
  CP <- CP[order(CP$variable, CP$country, CP$model), ]
  cp_path <- file.path(dirname(base), "delta_fields", "country_plateaus.csv")
  cp_tmp  <- tempfile(fileext = ".csv")
  write.csv(CP, cp_tmp, row.names = FALSE)
  publish_file(cp_tmp, cp_path)
  cat("Country plateaus ->", cp_path, "\n")
}
cat("PLOTTED -> ", pdf_final, "\n")
