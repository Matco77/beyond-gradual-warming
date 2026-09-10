# ISIMIP-vs-NAHosMIP attribution at matched delta_Sv: how much of the ssp126 impact looks like the
# response to an AMOC weakening of the same size.
# Author: Marco Bova
#
# Reads the two CENTRAL EU tables already on disk - no replay, no beta refit, nothing re-estimated.
# At every (model, bin) that BOTH branches kept:
#     A        = NAHosMIP effect  (AMOC weakening of that size on a piControl background,
#                                  i.e. NO greenhouse forcing)
#     T        = ISIMIP effect    (the ssp126 years whose own AMOC sits in that bin: greenhouse
#                                  warming AND the weakening that comes with it - NOT hosing)
#     residual = T - A            the part of the ssp126 impact that does not look like an AMOC response
#     share    = A / T            the AMOC-analogue share of the ssp126 impact
#
# FIRST-ORDER ONLY. Both assumptions behind the subtraction are violated to some degree:
#  1. ADDITIVITY. A is measured on a piControl background, so T - A treats the AMOC response as the
#     same whether it lands on a pre-industrial or on a warm climate. The indicators are threshold
#     hinges (crop 5/28 C, energy 15/24 C dead-bands) - the whole reason Phase 4 perturbs the DAILY
#     field - so the same delta on a warmer baseline does NOT produce the same indicator change.
#     The residual therefore carries the interaction term; it is not a clean "greenhouse-only" part.
#  2. delta_Sv FOOTING. NAHosMIP's delta_Sv is hos(y) - mean(piControl): unforced, stationary.
#     ISIMIP's is ssp126(y) - mean(historical 1850-2014), and the historical run carries real
#     forcing (plot_impact_common.R, Appendix J.5). Rows are matched on the INTEGER bin, and each
#     branch's own bin-mean delta_sv is carried through to the output so the residual mismatch is
#     visible rather than assumed away.
# Report as indicative. It is a comparison at matched AMOC level, not an identified attribution.
#
# Coverage is structurally small and scenario-dependent: NAHosMIP runs to -9.5 Sv, but ssp126 only
# reaches ~-3.5 Sv (6 model x bin cells) and ssp370 ~-5.5 Sv (7 cells) - a stronger forcing weakens
# the AMOC further, so it overlaps more of the hosing curve. Only the shallow bins can be compared
# either way. Nothing is extrapolated to fill the rest.
#
# A is NAHosMIP, so A does NOT change with the emissions scenario - only T does. Running this for
# both scenarios therefore contrasts the SAME AMOC contribution against two different greenhouse
# signals, which is the informative comparison: see how much of T the fixed A can account for.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/plot_impact_common.R"))

# ISIMIP emissions scenario, from ISIMIP_SCEN (default ssp126) - same convention as the impact
# scripts. The NAHosMIP side is scenario-independent.
ISCEN <- Sys.getenv("ISIMIP_SCEN", "ssp126")
ISFX  <- if (ISCEN == "ssp126") "" else paste0("_", ISCEN)

# below this |T| the ratio A/T is noise, not a share (0.1% in log units)
SHARE_FLOOR <- 1e-3

## ---- pair the two branches on the bins they share ----------------------------------------------
# INNER join on (model, bin_id, component): only bins present on both sides survive, which is
# exactly the level-for-level intersection isimip_bin_fields.R already enforced upstream.
pair <- function(a_file, i_file, extra = NULL) {
  k <- c("model", "bin_id", "component", extra)
  A <- fread(file.path(d, a_file)); I <- fread(file.path(d, sprintf(i_file, ISFX)))
  setnames(A, c("delta_sv", "n_years", "dln"), c("dsv_amoc",   "n_amoc",   "dln_amoc"))
  setnames(I, c("delta_sv", "n_years", "dln"), c("dsv_isimip", "n_isimip", "dln_isimip"))
  merge(A, I, by = k)
}

AT <- rbindlist(list(
  pair("amoc_impact_crop_eu.csv",   "isimip_bin_impact_crop_eu%s.csv")[,   branch := "crop (spec A)"],
  pair("amoc_impact_cropgw_eu.csv", "isimip_bin_impact_cropgw_eu%s.csv")[, branch := "crop (spec C)"],
  pair("amoc_impact_energy_eu.csv", "isimip_bin_impact_energy_eu%s.csv", "fuel")[, branch := "energy"]
), use.names = TRUE, fill = TRUE)
AT[, dln_residual := dln_isimip - dln_amoc]

## ---- totals per (branch, model, bin) and the share ----------------------------------------------
tot <- AT[, .(dsv_amoc = dsv_amoc[1], dsv_isimip = dsv_isimip[1], n_amoc = n_amoc[1],
              n_isimip = n_isimip[1], dln_amoc = sum(dln_amoc), dln_isimip = sum(dln_isimip)),
          by = .(branch, fuel, model, bin_id)]
tot[, dln_residual := dln_isimip - dln_amoc]
# A share is only meaningful when the two effects point the SAME way and T is not ~0. Where the AMOC
# effect opposes the ssp126 effect, "share" would print a plausible-looking number for a case where
# the AMOC is pushing against the total - left NA, the two % columns show what happened instead.
tot[, share := fifelse(abs(dln_isimip) >= SHARE_FLOOR & sign(dln_amoc) == sign(dln_isimip),
                       100 * dln_amoc / dln_isimip, NA_real_)]

## ---- self-check: the merge did not blow up, and the decomposition is exact ---------------------
stopifnot(nrow(AT) > 0,
          uniqueN(AT[, .(branch, fuel, model, bin_id, component)]) == nrow(AT),   # no cartesian dup
          all(abs(AT$dln_amoc  + AT$dln_residual  - AT$dln_isimip)  < 1e-12),
          all(abs(tot$dln_amoc + tot$dln_residual - tot$dln_isimip) < 1e-12))

for (f in c("dln_amoc", "dln_isimip", "dln_residual"))
  set(AT, j = sub("dln", "pct", f), value = pct(AT[[f]]))
for (f in c("dln_amoc", "dln_isimip", "dln_residual"))
  set(tot, j = sub("dln", "pct", f), value = pct(tot[[f]]))

setorder(AT,  branch, fuel, model, -bin_id, component)
setorder(tot, branch, fuel, model, -bin_id)
fwrite(AT,  file.path(d, sprintf("amoc_isimip_attribution%s.csv", ISFX)))
fwrite(tot, file.path(d, sprintf("amoc_isimip_attribution_total%s.csv", ISFX)))

## ---- report ------------------------------------------------------------------------------------
cat(sprintf("\n===== AMOC-analogue share of the %s impact, at matched delta_Sv =====\n", ISCEN))
cat("A = NAHosMIP (AMOC weakening only, piControl background, no greenhouse forcing)\n")
cat(sprintf("T = ISIMIP %s (greenhouse warming + the AMOC weakening it carries; NO hosing)\n", ISCEN))
cat("residual = T - A (carries the interaction, see header) | share = A/T\n")
cat(sprintf("bins compared: %d (model x bin), %.1f to %.1f Sv - %s never reaches deeper\n",
            uniqueN(tot[, .(model, bin_id)]), max(tot$bin_id), min(tot$bin_id), ISCEN))

for (br in unique(tot$branch)) for (fl in unique(tot[branch == br]$fuel)) {
  z <- tot[branch == br & (is.na(fuel) | fuel == fl)][order(model, -bin_id)]
  if (!nrow(z)) next
  cat(sprintf("\n-- %s%s --\n", br, if (is.na(fl)) "" else paste(",", fl)))
  cat(sprintf("  %-14s %6s %8s %8s %9s %9s %10s %8s\n",
              "model", "bin", "dSv(A)", "dSv(T)", "A %", "T %", "residual%", "share%"))
  for (i in seq_len(nrow(z)))
    cat(sprintf("  %-14s %6.1f %8.2f %8.2f %+8.2f%% %+8.2f%% %+9.2f%% %s\n",
                z$model[i], z$bin_id[i], z$dsv_amoc[i], z$dsv_isimip[i],
                z$pct_amoc[i], z$pct_isimip[i], z$pct_residual[i],
                if (is.na(z$share[i])) "      n/a" else sprintf("%7.0f%%", z$share[i])))
}

cat("\n!! The crop (spec A) rows inherit amoc_impact.R's warning: the crop NET is a specification\n")
cat("!! artefact (its sign is beta_gdd < 0, which does not survive spec A -> C). Read the crop\n")
cat("!! COMPONENTS in amoc_isimip_attribution.csv, not these crop totals. Energy has no such issue.\n")
cat("!! share = n/a means the two effects have opposite signs (AMOC pushes against the ssp126\n")
cat("!! total, so no fraction of it is 'due to' the AMOC) or |T| is below the noise floor.\n")
cat(sprintf("\nwrote amoc_isimip_attribution%s.csv (per component) and _total%s.csv\n", ISFX, ISFX))
