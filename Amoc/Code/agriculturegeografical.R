# Agriculture x geography: CropStatHarm (wheat/barley area, production, yield)
# + cropareaRegional (supplementary area source) bound to NUTS 2016 geometry.
# Author: Marco Bova
#
# CropStatHarm is harmonised to NUTS 2016 (Ronchetti et al. 2024, ESSD 16:1623),
# so region codes join 1:1 onto the NUTS 2016 shapefile - no crosswalk needed.

library(dplyr)
library(tidyr)
library(readr)
library(sf)

data_dir <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")

# --- CropStatHarm: area, production, yield -------------------------------
# UNIT encodes the variable: ha = area, t = production, t/ha = yield.
# Pivot those three rows into one row per region x crop x year.
csh <- read_delim(
  file.path(data_dir, "CropStatHarm_sub-national_16234acd-d17d-49de-bcb2-908a366a12be_2025.v01_39.csv"),
  delim = ";",
  col_types = cols_only(
    IDREGION = col_character(), CROP_NAME = col_character(),
    YEAR = col_integer(), VALUE = col_double(), UNIT = col_character()
  )
) |>
  mutate(UNIT = recode(UNIT, "ha" = "area_ha", "t" = "production_t", "t/ha" = "yield_t_ha")) |>
  pivot_wider(
    id_cols = c(IDREGION, CROP_NAME, YEAR),
    names_from = UNIT, values_from = VALUE
  ) |>
  transmute(
    NUTS_ID = IDREGION,
    crop    = CROP_NAME,
    year    = YEAR,
    area_ha, production_t, yield_t_ha
  )

# --- cropareaRegional: supplementary area-only source ---------------------
# Independent regional AREA panel (1980-2022, NUTS1-3) - NOT a spatial/gridded
# crop mask, it reports area at the same administrative levels as CropStatHarm.
# No published IDCROP legend; identity confirmed empirically by correlating
# each code's AREA against CropStatHarm's area_ha over all matching
# region-years (median ratio exactly 1.00, r > 0.9 for all four):
#   1 = Soft wheat   13 = Winter barley   41 = Durum wheat   42 = Total barley
#
# Why bother with a second source for the same crops: coverage differs from
# CropStatHarm in two ways that matter for an area regression -
#  - finer admin level for some countries: e.g. Poland has 17 NUTS2 regions in
#    CropStatHarm but 73 NUTS3 regions here (CropStatHarm's NSI source never
#    reported Polish wheat/barley below NUTS2). Also BE, NL, AT, PT, HR, DE.
#  - longer time series for regions both sources cover, e.g. CZ010 soft wheat
#    runs 1987-2022 here vs 1998-2022 in CropStatHarm.
# So it is used to ADD region-years CropStatHarm lacks, and to cross-check the
# region-years it has. CropStatHarm's value always wins on overlap - it is the
# peer-reviewed harmonised product; this is one of its raw ingredients.
car <- read_delim(
  file.path(data_dir, "cropareaRegional_b2ba07b7-73da-4d7f-96f6-4928b2f2b26e_2022.01_31.csv"),
  delim = ";",
  col_types = cols_only(IDCROP = col_integer(), IDREGION = col_character(),
                         YEAR = col_integer(), AREA = col_double())
) |>
  filter(nchar(IDREGION) %in% c(4, 5)) |>   # drop NUTS1, CropStatHarm never uses it
  transmute(
    NUTS_ID = IDREGION,
    year    = YEAR,
    crop    = recode(as.character(IDCROP), "1" = "Soft wheat", "13" = "Winter barley",
                      "41" = "Durum wheat", "42" = "Total barley"),
    area_alt = AREA
  ) |>
  pivot_wider(names_from = crop, values_from = area_alt) |>
  # complete the other two CropStatHarm crops from the two identities they satisfy
  mutate(
    `Total wheat`   = `Soft wheat` + `Durum wheat`,
    `Spring barley` = if_else(`Total barley` - `Winter barley` >= 0,
                               `Total barley` - `Winter barley`, NA_real_)
  ) |>
  pivot_longer(-c(NUTS_ID, year), names_to = "crop", values_to = "area_alt") |>
  filter(!is.na(area_alt))

# --- combine: CropStatHarm wins on overlap, cropareaRegional fills the rest
crop <- full_join(csh, car, by = c("NUTS_ID", "crop", "year")) |>
  mutate(
    area_source = case_when(
      !is.na(area_ha)  ~ "CropStatHarm",
      !is.na(area_alt) ~ "cropareaRegional",
      TRUE             ~ NA_character_    # neither source reports area here
    ),
    area_ha_alt = if_else(!is.na(area_ha) & !is.na(area_alt), area_alt, NA_real_),
    area_ha     = coalesce(area_ha, area_alt),
    nuts_level  = nchar(NUTS_ID) - 2L   # 2 chars = NUTS0, 4 = NUTS2, 5 = NUTS3
  ) |>
  select(NUTS_ID, nuts_level, crop, year, area_ha, production_t, yield_t_ha,
         area_source, area_ha_alt)

# NUTS 2016 region polygons, all levels, WGS84 (GDAL reads the .shp.zip directly)
nuts <- st_read(
  file.path(data_dir, "ref-nuts-2016-03m.shp/NUTS_RG_03M_2016_4326.shp.zip"),
  quiet = TRUE
) |>
  select(NUTS_ID, LEVL_CODE, CNTR_CODE, NAME_LATN)

# guard the harmonisation assumption: every region in either source must have a 2016 polygon
stopifnot(all(crop$NUTS_ID %in% nuts$NUTS_ID))

# keep only the polygons that appear in the combined data (one row per region)
nuts <- filter(nuts, NUTS_ID %in% crop$NUTS_ID)

cat("panel rows:", nrow(crop), "| regions:", n_distinct(crop$NUTS_ID),
    "| years:", paste(range(crop$year), collapse = "-"), "\n")
print(count(crop, area_source))

# full panel (region x crop x year), both sources, both levels - kept for
# audit/cross-validation (area_source, area_ha_alt), not meant for direct
# regression use since some regions still double-count NUTS2 vs NUTS3 below
write_csv(crop, file.path(data_dir, "6.CropStatHarm_prepared.csv"))

# --- regression panel: one geography per region, no double-counting --------
# nuts_level mixes NUTS0/2/3, and some regions now have both a NUTS2 row (from
# CropStatHarm) and NUTS3 children (from cropareaRegional) for the same
# crop-year - using both would double-count that area. Coverage is uneven (a
# NUTS2 region may have some NUTS3 children reported and others not), so dedup
# per NUTS2-parent, not per country: drop the NUTS2 row only where at least
# one NUTS3 child actually exists for that same crop-year; keep it otherwise.
crop_regression <- crop |>
  mutate(parent2 = if_else(nuts_level == 3, substr(NUTS_ID, 1, 4), NUTS_ID)) |>
  group_by(parent2, crop, year) |>
  filter(n() == 1 | nuts_level == 3) |>
  ungroup() |>
  select(-parent2)

cat("regression panel rows:", nrow(crop_regression), "(dropped",
    nrow(crop) - nrow(crop_regression), "superseded NUTS2 parent rows) | regions:",
    n_distinct(crop_regression$NUTS_ID), "\n")

write_csv(crop_regression, file.path(data_dir, "6.CropStatHarm_regression_panel.csv"))
st_write(filter(nuts, NUTS_ID %in% crop_regression$NUTS_ID),
         file.path(data_dir, "6.NUTS2016_cropregions.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)

# regression-ready sf panel, built on demand:
# panel <- read_csv(file.path(data_dir, "6.CropStatHarm_regression_panel.csv")) |>
#   left_join(st_read(file.path(data_dir, "6.NUTS2016_cropregions.gpkg")), by = "NUTS_ID") |>
#   st_as_sf()
