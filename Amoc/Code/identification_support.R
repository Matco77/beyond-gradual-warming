# Does the scenario stay inside the variation that identifies beta?
# Author: Marco Bova
#
# Two questions the coefficient tables cannot answer.
#
# (1) HOW MUCH VARIATION IS LEFT after the fixed effects. A benchmark with unit FE, year FE and
#     unit-specific quadratic trends throws away the cross-section on purpose - it is confounded -
#     and identifies beta from year-to-year wobble only. If the wobble is tiny, the coefficient is
#     precise but uninformative. This is what decides how much weight the GDD-window column can
#     carry: there the window is DEFINED by accumulating fixed thermal time, so gdd inside it barely
#     moves, and the fixed effects absorb what little is left.
#
# (2) HOW FAR OUTSIDE IT THE SCENARIO REACHES, measured in units of that same within sd. A response
#     function evaluated many sd beyond its estimation support is an extrapolation of the FUNCTIONAL
#     FORM, not a prediction - and this applies to the energy branch too, not only to crop.
d  <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
RG <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/regression.R")
grab <- function(f, want) { e <- new.env(parent = globalenv())
  for (x in parse(f)) if (is.call(x) && identical(as.character(x[[1]]),"<-") && is.name(x[[2]]) &&
      as.character(x[[2]]) %in% want) eval(x, e); e }
R <- grab(RG, c("common_obs","crop_set","load_crop","f_tr")); FE <- "NUTS_ID[t, t2] + year"
cp9 <- R$load_crop("9.crop_panel_nuts3_estimation.csv")
gw  <- fread(file.path(d,"11.crop_weather_gdd_window.csv"))[, precip2 := precip^2][]
f_gw <- as.formula(paste("ln_yield ~ gdd + heat + frost + precip + precip2 |", FE))

cat("=== IDENTIFYING VARIATION IN gdd: how much survives the fixed effects? ===\n")
cat("within sd = sd of gdd after partialling out", FE, "\n")
cat("share = within variance / raw variance (how much of gdd's variation is used)\n\n")
cat(sprintf("%-14s %26s %26s %8s\n", "", "-------- A: Mar-Jul --------", "------ C: GDD-window -------", ""))
cat(sprintf("%-14s %8s %8s %8s %8s %8s %8s %8s\n", "crop",
            "mean", "within", "share", "mean", "within", "share", "C/A"))
for (k in R$crop_set) {
  s <- R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr)
  z <- R$common_obs(merge(s[, .(NUTS_ID,cntr,year,ln_yield,t,t2)], gw, by = c("NUTS_ID","year")),
                    list(f_gw), ~cntr)
  wA <- resid(feols(as.formula(paste("gdd_mj ~ 1 |", FE)), s, notes = FALSE))
  wC <- resid(feols(as.formula(paste("gdd ~ 1 |", FE)),    z, notes = FALSE))
  cat(sprintf("%-14s %8.0f %8.1f %7.1f%% %8.0f %8.1f %7.1f%% %8.3f\n", k,
      mean(s$gdd_mj), sd(wA), 100*var(wA)/var(s$gdd_mj),
      mean(z$gdd),    sd(wC), 100*var(wC)/var(z$gdd), sd(wC)/sd(wA)))
}

cat("\n=== and how much does the SCENARIO move gdd, under each spec? ===\n")
scA <- fread(file.path(d,"scenario_bins_crop_window.csv.gz"))[window=="Mar-Jul" & model=="IPSL-CM6A-LR"]
scC <- fread(file.path(d,"scenario_bins_crop_gddwin.csv.gz"))[model=="IPSL-CM6A-LR"]
hA  <- fread(file.path(d,"8.eobs_nuts3_crop_weather_window.csv"))[window=="Mar-Jul", .(NUTS_ID,year,h=gdd)]
hC  <- gw[, .(NUTS_ID, year, h = gdd)]
dA <- merge(scA, hA, by=c("NUTS_ID","year"))[, .(d = mean(gdd - h)), by=bin_id]
dC <- merge(scC, hC, by=c("NUTS_ID","year"))[, .(d = mean(gdd - h)), by=bin_id]
m  <- merge(dA, dC, by="bin_id", suffixes=c("_A","_C"))[order(bin_id)]
sA <- sd(resid(feols(as.formula(paste("gdd_mj ~ 1 |", FE)), R$common_obs(cp9[crop=="Soft wheat"], list(R$f_tr), ~cntr), notes=FALSE)))
zz <- R$common_obs(merge(R$common_obs(cp9[crop=="Soft wheat"], list(R$f_tr), ~cntr)[, .(NUTS_ID,cntr,year,ln_yield,t,t2)], gw, by=c("NUTS_ID","year")), list(f_gw), ~cntr)
sC <- sd(resid(feols(as.formula(paste("gdd ~ 1 |", FE)), zz, notes=FALSE)))
cat(sprintf("  %7s %12s %12s   | in units of the within sd (soft wheat)\n","bin Sv","dGDD A","dGDD C"))
for (i in seq_len(nrow(m))) cat(sprintf("  %7.1f %+12.1f %+12.1f   | A %+5.2f sd   C %+5.2f sd\n",
    m$bin_id[i], m$d_A[i], m$d_C[i], m$d_A[i]/sA, m$d_C[i]/sC))

d  <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
RG <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/regression.R")
grab <- function(f,w){e<-new.env(parent=globalenv()); for(x in parse(f)) if(is.call(x)&&identical(as.character(x[[1]]),"<-")&&
  is.name(x[[2]])&&as.character(x[[2]])%in%w) eval(x,e); e}
R <- grab(RG, "common_obs")
ep <- fread(file.path(d,"16.energy_panel_estimation.csv")); ep[, `:=`(t=year-2007L, t2=(year-2007L)^2)]
FE <- "country_id[t, t2] + year"
efit <- function(dt,wx){f<-function(fe)as.formula(paste("ln_energy_pc ~",wx,"+ ln_price |",fe))
  fs<-list(f("country_id + year"),f("country_id[t, t2]"),f("country_id[t, t2] + year"))
  R$common_obs(dt,fs,~country_id)}
sc <- fread(file.path(d,"scenario_bins_energy_country.csv.gz"))[model=="IPSL-CM6A-LR"]
h  <- fread(file.path(d,"13.eobs_country_energy_weather_weighted.csv"))
cat("=== ENERGY: does the scenario stay inside the variation that identifies beta? ===\n")
for (z in list(c("Electricity","hdd_calendar + cdd_jja","hdd_calendar","cdd_jja"),
               c("Natural gas","hdd_octmar","hdd_octmar",NA))) {
  s <- efit(ep[fuel==z[1]], z[2]); vars <- na.omit(c(z[3],z[4]))
  cat(sprintf("\n %s (N=%d, %d countries)\n", z[1], nrow(s), uniqueN(s$country_id)))
  for (v in vars) {
    w  <- resid(feols(as.formula(paste(v,"~ 1 |",FE)), s, notes=FALSE))
    sw <- sd(w)
    dd <- merge(sc, h[, c("cntr","year",v), with=FALSE], by=c("cntr","year"))
    dd[, dlt := get(paste0(v,".x")) - get(paste0(v,".y"))]
    g  <- dd[, .(m=mean(dlt,na.rm=TRUE)), by=bin_id][order(bin_id)]
    cat(sprintf("   %-13s within sd %6.1f (%.1f%% of raw var) | scenario shift: %s\n", v, sw,
        100*var(w)/var(s[[v]]),
        paste(sprintf("%.1fSv %+.1fsd", g$bin_id, g$m/sw), collapse="  ")))
  }
}
