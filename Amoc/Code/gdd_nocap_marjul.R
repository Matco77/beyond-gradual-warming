# uncapped GDD on the FIXED Mar-Jul window: isolates the 28C cap from the moving window
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_engine.R"))
S <- eobs_setup("crop")
gdd_nocap <- function(tg) pmax(tg - 5, 0)          # gdd_window_daily.R:43,57 - no 28C cap
nc <- nc_open(ncf("tg")); out <- list()
for (y in YRS_CROP) {
  ti <- which(.ax$yr == y); mm <- .ax$mo[ti]
  a  <- ncvar_get(nc, "tg", start = c(1,1,ti[1]), count = c(-1,-1,length(ti)))
  dim(a) <- c(.ax$nlon * .ax$nlat, length(ti)); TG <- a[S$ncdf_row, , drop = FALSE]
  M  <- round(.ind_month(S, TG, mm, gdd_nocap), 3)          # same two-step rounding as 7./8.
  keep <- as.integer(colnames(M)) %in% 3:7
  out[[length(out)+1L]] <- data.table(NUTS_ID = S$grp, year = y,
                                      gdd_nocap_mj = round(rowSums(M[, keep, drop = FALSE]), 3))
  cat("."); flush.console()
}
nc_close(nc)
fwrite(rbindlist(out), file.path(d, "11b.gdd_nocap_marjul.csv"))
cat("\nwrote 11b.gdd_nocap_marjul.csv\n")
