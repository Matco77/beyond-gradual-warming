#Datasets Preparation for "Thesis"
#Author: Marco Bova

#library


###############
# 1. FAOSTAT
###############

FAOSTAT <- read.csv("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/FAOSTAT_data_en_5-4-2026 (1).csv")
#modify as such that the structure is the following country crop year yield production area


###############
# 2. E-OBS
###############
library(dplyr)
library(ggplot2)

nahos_root <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/NAHosMIP"
ll_dir <- file.path(nahos_root, "HadGEM3-GC31-LL")

out_dir <- file.path(nahos_root, "tables_visual_diagnostics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

inspect_nc_table <- function(file) {
  
  nc <- nc_open(file)
  
  dim_table <- data.frame(
    file = basename(file),
    dimension = names(nc$dim),
    length = sapply(nc$dim, function(x) x$len),
    units = sapply(nc$dim, function(x) x$units)
  )
  
  var_table <- data.frame(
    file = basename(file),
    variable = names(nc$var),
    units = sapply(nc$var, function(x) x$units),
    dimensions = sapply(nc$var, function(x) {
      paste(sapply(x$dim, function(d) d$name), collapse = " x ")
    })
  )
  
  nc_close(nc)
  
  list(
    dimensions = dim_table,
    variables = var_table
  )
}

info_1 <- inspect_nc_table(g01_tas_1)
info_2 <- inspect_nc_table(g01_tas_2)

dimensions_table <- bind_rows(info_1$dimensions, info_2$dimensions)
variables_table <- bind_rows(info_1$variables, info_2$variables)

print(dimensions_table)
print(variables_table)

write_csv(dimensions_table, file.path(out_dir, "g01_tas_dimensions_table.csv"))
write_csv(variables_table, file.path(out_dir, "g01_tas_variables_table.csv"))

fig_dir <- file.path(nahos_root, "figures_visual_diagnostics")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# g01-hos files: HadGEM3-GC31-LL
# -------------------------

g01_tas_1 <- file.path(ll_dir, "tas_Amon_HadGEM3-GC31-LL_g01-hos_r1i1p1f1_gn_205001-205412.nc")
g01_tas_2 <- file.path(ll_dir, "tas_Amon_HadGEM3-GC31-LL_g01-hos_r1i1p1f1_gn_205501-214912.nc")

g01_tasmin_1 <- file.path(ll_dir, "tasmin_Amon_HadGEM3-GC31-LL_g01-hos_r1i1p1f1_gn_205001-205412.nc")
g01_tasmin_2 <- file.path(ll_dir, "tasmin_Amon_HadGEM3-GC31-LL_g01-hos_r1i1p1f1_gn_205501-214912.nc")

g01_tasmax_1 <- file.path(ll_dir, "tasmax_Amon_HadGEM3-GC31-LL_g01-hos_r1i1p1f1_gn_205001-205412.nc")
g01_tasmax_2 <- file.path(ll_dir, "tasmax_Amon_HadGEM3-GC31-LL_g01-hos_r1i1p1f1_gn_205501-214912.nc")

g01_pr_1 <- file.path(ll_dir, "pr_Amon_HadGEM3-GC31-LL_g01-hos_r1i1p1f1_gn_205001-205412.nc")
g01_pr_2 <- file.path(ll_dir, "pr_Amon_HadGEM3-GC31-LL_g01-hos_r1i1p1f1_gn_205501-214912.nc")


# -------------------------
# piControl files: HadGEM3-GC31-LL
# -------------------------

pic_tas_file <- file.path(ll_dir, "tas_Amon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_205001-214912.nc")
pic_tasmin_file <- file.path(ll_dir, "tasmin_Amon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_205001-214912.nc")
pic_tasmax_file <- file.path(ll_dir, "tasmax_Amon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_205001-214912.nc")
pic_pr_file <- file.path(ll_dir, "pr_Amon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_205001-214912.nc")

# -------------------------
# Open g01-hos
# -------------------------

g01_tas <- c(
  rast(paste0("NETCDF:", g01_tas_1, ":tas")),
  rast(paste0("NETCDF:", g01_tas_2, ":tas"))
)

g01_tasmin <- c(
  rast(paste0("NETCDF:", g01_tasmin_1, ":tasmin")),
  rast(paste0("NETCDF:", g01_tasmin_2, ":tasmin"))
)

g01_tasmax <- c(
  rast(paste0("NETCDF:", g01_tasmax_1, ":tasmax")),
  rast(paste0("NETCDF:", g01_tasmax_2, ":tasmax"))
)

g01_pr <- c(
  rast(paste0("NETCDF:", g01_pr_1, ":pr")),
  rast(paste0("NETCDF:", g01_pr_2, ":pr"))
)


# -------------------------
# Open piControl
# -------------------------

pic_tas <- rast(paste0("NETCDF:", pic_tas_file, ":tas"))
pic_tasmin <- rast(paste0("NETCDF:", pic_tasmin_file, ":tasmin"))
pic_tasmax <- rast(paste0("NETCDF:", pic_tasmax_file, ":tasmax"))
pic_pr <- rast(paste0("NETCDF:", pic_pr_file, ":pr"))


###############
# 6. CropStatHarm + NUTS 2016 geometry -> see agriculturegeografical.R
###############
