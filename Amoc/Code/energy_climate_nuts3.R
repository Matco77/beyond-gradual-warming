# Energy climate covariates (HDD/CDD) at NUTS3, POPULATION-weighted - consistent with the
# country energy module (energyweather), just at finer spatial resolution.
#   w_pop_{g,r} = pop(cell g inside NUTS3 r) / pop(r), from the 1km Census-2021 grid.
#   HDD_{g,d}=1(T<15)(18-T), CDD_{g,d}=1(T>=24)(T-21), T=tg (slide 19). Mean over 2000-2019.
# Writes Amoc/figures/map_energy_climate_nuts3.png (HDD calendar, HDD Oct-Mar, CDD JJA).
#
# VISUALISATION ONLY: the energy regression stays at country level - Eurostat household energy
# has no subnational series, so NUTS3 HDD/CDD has nothing to regress against. Population-weighting
# barely moves a small NUTS3 (~0.98 corr with area-weighting); it matters most in mountain regions,
# where people sit in the warm valleys and an area-mean over-counts the cold peaks.
# Author: Marco Bova
suppressMessages({library(sf); library(terra); library(ncdf4); library(Matrix); library(data.table)
                  library(ggplot2); library(grid)})
sf_use_s2(FALSE)
d   <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
out <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/figures")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
ncf <- function(v) file.path(d, sprintf("EOBS/%s_ens_mean_0.25deg_reg_v31.0e.nc", v)); YRS <- 2000:2019

## --- population weights: 1km pop -> (E-OBS cell, NUTS3) -----------------------
pop    <- rast(file.path(d, "GISCO_population_grid/ESTAT_OBS-VALUE-T_2021_V2.tiff"))
tmpl   <- rast(ncf("tg"))[[1]]; cellid <- init(tmpl, "cell")
cell_1km <- project(cellid, pop, method = "near")
nuts <- st_read(file.path(d, "ref-nuts-2016-03m.shp/NUTS_RG_03M_2016_4326.shp.zip"), quiet = TRUE)
n3sf <- st_transform(nuts[nuts$LEVL_CODE == 3, "NUTS_ID"], 3035); n3sf$idx <- seq_len(nrow(n3sf))
r3   <- rasterize(vect(n3sf), pop, field = "idx")
z <- data.table(pop = values(pop)[,1], cell = values(cell_1km)[,1], idx = values(r3)[,1])
z <- z[!is.na(pop) & pop > 0 & !is.na(cell) & !is.na(idx)][, nuts := n3sf$NUTS_ID[idx]]
wp <- z[, .(pop = sum(pop)), by = .(cell, nuts)][, w := pop / sum(pop), by = nuts][]
stopifnot(max(abs(wp[, sum(w), by = nuts]$V1 - 1)) < 1e-9)   # weights renormalise per region

## --- fixed weight matrix W [n_nuts x n_cells] --------------------------------
cells <- sort(unique(wp$cell)); nutsv <- sort(unique(wp$nuts))
W <- sparseMatrix(i = match(wp$nuts, nutsv), j = match(wp$cell, cells), x = wp$w,
                  dims = c(length(nutsv), length(cells)))

## --- map crosswalk cells -> ncdf4 [lon,lat] rows -----------------------------
r  <- rast(ncf("tg")); xy <- xyFromCell(r[[1]], cells)
nc0 <- nc_open(ncf("tg")); nc_lon <- ncvar_get(nc0, "longitude"); nc_lat <- ncvar_get(nc0, "latitude")
tvec <- as.Date(ncvar_get(nc0, "time"), origin = "1950-01-01"); nc_close(nc0); nlon <- length(nc_lon)
loni <- match(round(xy[,1],3), round(nc_lon,3)); lati <- match(round(xy[,2],3), round(nc_lat,3))
stopifnot(!anyNA(loni), !anyNA(lati)); ncdf_row <- (lati-1)*nlon + loni
yr <- as.integer(format(tvec,"%Y")); mo <- as.integer(format(tvec,"%m"))
wmean <- function(v) { p <- !is.na(v); v[!p] <- 0                    # pop-weighted NUTS3 mean, NA-safe
  m <- as.matrix(W %*% v) / as.matrix(W %*% (p*1)); m[!is.finite(m)] <- NA; m }

## --- daily HDD/CDD -> NUTS3, 2000-2019 mean ----------------------------------
acc <- data.table(NUTS_ID = nutsv, hdd_calendar = 0, hdd_octmar = 0, cdd_jja = 0, valid = 0L, n = 0L)
nc <- nc_open(ncf("tg"))
for (y in YRS) {
  ti <- which(yr == y); mm <- mo[ti]
  a <- ncvar_get(nc, "tg", start = c(1,1,ti[1]), count = c(-1,-1,length(ti)))
  dim(a) <- c(nlon*length(nc_lat), length(ti)); M <- a[ncdf_row,,drop=FALSE]
  hdd <- wmean(ifelse(M < 15, 18 - M, 0)); cdd <- wmean(ifelse(M >= 24, M - 21, 0))
  acc[, `:=`(hdd_calendar = hdd_calendar + rowSums(hdd, na.rm = TRUE),
             hdd_octmar   = hdd_octmar + rowSums(hdd[, mm %in% c(10,11,12,1,2,3), drop=FALSE], na.rm=TRUE),
             cdd_jja      = cdd_jja + rowSums(cdd[, mm %in% 6:8, drop=FALSE], na.rm=TRUE),
             valid = valid + rowSums(is.finite(hdd)), n = n + 1L)]
  cat(y, "done\n"); flush.console()
}
nc_close(nc)
acc[, `:=`(hdd_calendar = hdd_calendar/n, hdd_octmar = hdd_octmar/n, cdd_jja = cdd_jja/n)]
# regions with no observed weather (no census pop, or population only on E-OBS sea cells) -> NA
for (v in c("hdd_calendar","hdd_octmar","cdd_jja")) acc[valid == 0, (v) := NA_real_]

## --- maps (NUTS3) ------------------------------------------------------------
n0 <- st_transform(nuts[nuts$LEVL_CODE == 0, "NUTS_ID"], 3035)
g  <- merge(st_transform(n3sf, 3035), acc, by = "NUTS_ID")
EUR <- coord_sf(xlim = c(2.5e6,6.1e6), ylim = c(1.4e6,5.4e6), expand = FALSE)
mk <- function(fill, title, opt, dir = 1) ggplot() +
  geom_sf(data = n0, fill = "grey92", color = NA) +
  geom_sf(data = g, aes(fill = .data[[fill]]), color = NA) +
  scale_fill_viridis_c(option = opt, direction = dir, na.value = "grey85", name = "C.days",
                       guide = guide_colorbar(barwidth = 0.6, barheight = 4)) +
  EUR + labs(title = title) + theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 10.5, hjust = 0.01, margin = margin(b = 2)),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        plot.margin = margin(4,4,4,4), plot.background = element_rect(fill = "white", color = NA))
png(file.path(out, "map_energy_climate_nuts3.png"), width = 2100, height = 850, res = 150)
pushViewport(viewport(layout = grid.layout(1, 3)))
print(mk("hdd_calendar","HDD (calendar, base 18/15C) - NUTS3, pop-weighted","G",-1), vp = viewport(layout.pos.col=1))
print(mk("hdd_octmar",  "HDD (Oct-Mar heating) - NUTS3, pop-weighted","G",-1),        vp = viewport(layout.pos.col=2))
print(mk("cdd_jja",     "CDD (Jun-Aug cooling, base 21/24C) - NUTS3, pop-weighted","B"), vp = viewport(layout.pos.col=3))
dev.off()
cat("wrote map_energy_climate_nuts3.png | NUTS3 regions:", acc[!is.na(hdd_calendar), .N], "\n")
