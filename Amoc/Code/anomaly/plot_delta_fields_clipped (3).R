# ============================================================
# PLOT delta fields  (revised: CLIPPED scale + fingerprint readout)
#
# Two fixes over the first version:
#  (1) the shared scale is now the 98th percentile of |delta|, not
#      the max, so a few extreme cells (domain edge / cold-blob core
#      / reconstruction overshoot) no longer wash out the dominant
#      ~2 K cooling. Values beyond the scale are clipped to the end
#      colours. Scale is still SHARED + SYMMETRIC per variable, so
#      cross-model magnitudes stay comparable.
#  (2) each panel title now reports NW-minus-Med: the AMOC
#      fingerprint number. Negative => NW-amplified cooling = AMOC,
#      not a uniform offset.
#
# Diverging palette centred on 0 (blue = cooling/drying). Same
# Europe window, base-R + fields, coastlines via maps.
# ============================================================

library(ncdf4)
library(fields)
have_maps <- requireNamespace("maps", quietly = TRUE)

## rename-proof: source amoc_common.R from THIS script's own folder
## (keeps working if the folder is moved/renamed; run via `Rscript "<this file>"`).
.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)   # Rscript encodes spaces in --file= as ~+~
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: CLIP_Q, nw_box, med_box, publish_file, panel_grid,
## pal_for, wmean, box_mean.

delta_dir <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/delta_fields"
outdir    <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

panels <- list(
  c(model = "EC-Earth3",       proto = "g01"),
  c(model = "EC-Earth3",       proto = "u03"),
  c(model = "HadGEM3-GC31-LL", proto = "g01"),
  c(model = "HadGEM3-GC31-LL", proto = "u03"),
  c(model = "HadGEM3-GC31-MM", proto = "g01"),
  c(model = "HadGEM3-GC31-MM", proto = "u03"),
  c(model = "IPSL-CM6A-LR",    proto = "u03")   # u03 only (no g01)
)
vars <- c("tas", "pr", "tasmax", "tasmin")

# reads the exported delta_*_lastthird.nc, whose pr is ALREADY in mm/day
# (converted once in amoc_common.R's read_europe_cube before export) --
# so no x86400 here, and the mm/day legend below is correct.
read_delta <- function(path, vn) {
  nc <- nc_open(path); lon <- as.numeric(ncvar_get(nc, "lon"))
  lat <- as.numeric(ncvar_get(nc, "lat")); f <- ncvar_get(nc, vn); nc_close(nc)
  f[f <= -9990] <- NA; list(lon = lon, lat = lat, f = f)
}
delta_path <- function(model, proto, vn)
  file.path(delta_dir, sprintf("delta_%s_%s_%s_lastthird.nc", model, proto, vn))

pdf_final <- file.path(outdir, "delta_fields_panels_clipped.pdf")
pdf_tmp   <- tempfile(fileext = ".pdf")        # render locally, publish to OneDrive at the end
pdf(pdf_tmp, width = 11, height = 8)

## one page = 6 panels on a shared, symmetric, 98th-pct-clipped scale
draw_page <- function(loaded, vn, unit, head_prefix) {
  zmax <- as.numeric(quantile(abs(unlist(lapply(loaded, function(d) as.vector(d$f)))), CLIP_Q, na.rm = TRUE))
  if (!is.finite(zmax) || zmax == 0) zmax <- 1
  pal <- pal_for(vn)
  par(mfrow = panel_grid(length(loaded)), oma = c(0, 0, 3, 0), mar = c(3.2, 3.2, 2.4, 4.2), mgp = c(2, 0.6, 0))
  for (d in loaded) {
    fdisp <- pmin(pmax(d$f, -zmax), zmax)          # clip outliers to scale ends
    image.plot(d$lon, d$lat, fdisp, zlim = c(-zmax, zmax), col = pal,
               xlab = "lon", ylab = "lat", main = d$label,
               legend.lab = unit, legend.line = 2.3, cex.main = 0.92)
    if (have_maps) maps::map("world", add = TRUE, interior = FALSE, lwd = 0.5)
    abline(h = 0, v = 0, col = "grey80", lty = 3)
  }
  dirlab <- if (vn == "pr") "brown = drying, teal = wetting" else "blue = cooling"
  fplab  <- if (vn == "pr") "NW\u2212Med shown for reference" else "NW\u2212Med<0 = AMOC fingerprint"
  mtext(sprintf("%s \u2014 %s [%s]   |   shared scale clipped at %.0fth pct = \u00B1%.2f,  %s,  %s",
                head_prefix, vn, unit, CLIP_Q * 100, zmax, dirlab, fplab),
        outer = TRUE, cex = 0.95, font = 2)
}

for (vn in vars) {

  loaded <- Filter(Negate(is.null), lapply(panels, function(pm) {
    p <- delta_path(pm["model"], pm["proto"], vn); if (!file.exists(p)) return(NULL)
    d <- read_delta(p, vn)
    d$model <- pm["model"]; d$proto <- pm["proto"]
    d$mean <- wmean(d$f, d$lon, d$lat)
    d$grad <- box_mean(d$f, d$lon, d$lat, nw_box) - box_mean(d$f, d$lon, d$lat, med_box)
    d$label <- sprintf("%s  %s   (mean %+.2f,  NW\u2212Med %+.2f)", pm["model"], pm["proto"], d$mean, d$grad)
    d
  }))
  if (!length(loaded)) next

  unit <- if (vn == "pr") "mm/day" else "K"
  draw_page(loaded, vn, unit, "Late-run (plateau) delta")

  ## pr gets a SECOND page: the same delta as a PERCENTAGE of the control
  ## annual-mean precipitation (the literature-standard relative view).
  ## % = 100 * delta / clim; clim = annual mean of the piControl climatology
  ## (equal month weights), same Europe crop, mm/day on both sides.
  if (vn == "pr") {
    pct <- lapply(loaded, function(d) {
      cl <- pic_annual_field(d$model, "pr")$clim
      stopifnot(all(dim(cl) == dim(d$f)))
      d$f <- 100 * d$f / cl
      d$f[cl < 0.1] <- NA   # % undefined where control precip ~0 (guard; none expected in Europe window)
      d$mean <- wmean(d$f, d$lon, d$lat)
      d$grad <- box_mean(d$f, d$lon, d$lat, nw_box) - box_mean(d$f, d$lon, d$lat, med_box)
      d$label <- sprintf("%s  %s   (mean %+.1f%%,  NW\u2212Med %+.1f%%)", d$model, d$proto, d$mean, d$grad)
      d
    })
    draw_page(pct, vn, "%", "Late-run (plateau) delta, RELATIVE (% of control annual-mean precip)")
  }
}

dev.off()
publish_file(pdf_tmp, pdf_final)
cat("PLOTTED -> ", pdf_final, "\n")
