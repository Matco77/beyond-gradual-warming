# Crop long-difference (Burke & Emerick 2016) — a LONG-RUN BOUND on the weather response.
# The panel benchmark (regression.R) identifies the SHORT-RUN response to transient weather. A
# long-difference regresses the long change in mean log yield on the long change in mean climate
# ACROSS NUTS3, so cross-sectional variation in climate *trends* identifies a response that lets
# agents partly adapt (rotations, varieties, inputs). Reported as a BOUND, not a point estimate:
# a 35-year record with a modest trend leaves it imprecise. Energy is NOT done this way — only 29
# cross-sectional units give too little identifying variation (PIPELINE.md §10.1).
#
# Window choice: early 1999-2005 vs late 2016-2022 (crop panel is thin pre-1999; 2023 has 94 obs).
# Keep NUTS3 with >= MINYRS observed years in BOTH windows so each cell's means are well-measured.
# Author: Marco Bova
suppressMessages({library(fixest); library(data.table)})
d   <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
res <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/results")
save_tab <- function(x, f) writeLines(capture.output(print(x)), file.path(res, f))

crop_set <- c("Soft wheat", "Durum wheat", "Spring barley", "Winter barley")
EARLY <- 1999:2005; LATE <- 2016:2022; MINYRS <- 4L
wx <- c("gdd_mj", "heat_mj", "frost_mj", "precip_mj", "precip_mj2")

cp <- fread(file.path(d, "9.crop_panel_nuts3_estimation.csv"))[crop %in% crop_set]
cp[, precip_mj2 := precip_mj^2]                                  # same quadratic as the panel
cp[, win := fifelse(year %in% EARLY, "e", fifelse(year %in% LATE, "l", NA_character_))]
cp <- cp[!is.na(win)]

# per NUTS3 x crop x window: mean of each variable + mean area (long-difference weight)
agg <- cp[, c(.(n = .N, area = mean(area_ha)),
              lapply(.SD, mean)), by = .(NUTS_ID, cntr, crop, win),
          .SDcols = c("ln_yield", wx)]
agg <- agg[n >= MINYRS]
w <- dcast(agg, NUTS_ID + cntr + crop ~ win, value.var = c("ln_yield", wx, "area"))
w <- w[complete.cases(w[, .(ln_yield_e, ln_yield_l)])]          # keep cells present in both windows
for (v in c("ln_yield", wx)) w[, (paste0("D_", v)) := get(paste0(v, "_l")) - get(paste0(v, "_e"))]
w[, area_w := (area_e + area_l) / 2]                            # cropland weight (Burke-Emerick)

f_ld <- as.formula(paste("D_ln_yield ~", paste(paste0("D_", wx), collapse = " + "), "| cntr"))
# panel benchmark for side-by-side (short-run), on the same crop set
cp[, `:=`(t = year - 2005L, t2 = (year - 2005L)^2)]
f_tr <- ln_yield ~ gdd_mj + heat_mj + frost_mj + precip_mj + precip_mj2 | NUTS_ID[t, t2] + year

for (k in crop_set) {
  wk <- w[crop == k]
  m_ld <- feols(f_ld, wk, weights = ~area_w, cluster = ~cntr)   # long-difference (long-run bound)
  m_sr <- feols(f_tr, cp[crop == k], cluster = ~cntr)           # panel benchmark (short-run)
  tab <- etable(m_sr, m_ld,
                headers = c("Panel (short-run)", "Long-diff (long-run bound)"),
                dict = c(D_gdd_mj = "gdd_mj", D_heat_mj = "heat_mj", D_frost_mj = "frost_mj",
                         D_precip_mj = "precip_mj", D_precip_mj2 = "precip_mj2",
                         D_ln_yield = "ln_yield"),
                fitstat = ~ n + r2)
  cat("\n===== CROP long-difference:", k, "  (cells:", nrow(wk), ") =====\n"); print(tab)
  save_tab(tab, paste0("crop_", gsub(" ", "_", k), "_long_difference.txt"))
}
cat("\nwrote long-difference tables ->", res, "\n")
