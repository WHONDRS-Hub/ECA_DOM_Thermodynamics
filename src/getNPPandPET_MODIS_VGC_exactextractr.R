library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(stringr)

# ---- USER SETTINGS ----
ws_path   <- "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Geospatial_variables/Example_Code/shape/v4_RCSFA.shp"
id_col    <- "comid"
modis_dir <- "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Geospatial_variables/Example_Code/data/MODIS 2022-2023"

out_annual_csv <- sprintf("watershed_pet_et_npp_annual_2022_2023_%s.csv", Sys.Date())
out_mean_csv   <- sprintf("watershed_pet_et_npp_samplingwindow_mean_%s.csv", Sys.Date())

get_year <- function(x) as.integer(str_match(basename(x), "doy(\\d{4})")[,2])

# ---- read watersheds ----
ws_sf <- st_read(ws_path, quiet = TRUE)
stopifnot(id_col %in% names(ws_sf))

# ---- list files ----
files_all <- list.files(modis_dir, pattern="\\.tif$", full.names=TRUE, recursive=TRUE)

pet_files <- files_all[str_detect(files_all, "MOD16A3GF\\.061_PET_500m_")]
et_files  <- files_all[str_detect(files_all, "MOD16A3GF\\.061_ET_500m_")]
npp_files <- files_all[str_detect(files_all, "MOD17A3HGF\\.061_Npp_500m_")]

# keep only 2022–2023
keep_years <- 2022:2023
pet_files <- pet_files[get_year(pet_files) %in% keep_years]
et_files  <- et_files [get_year(et_files)  %in% keep_years]
npp_files <- npp_files[get_year(npp_files) %in% keep_years]

# sort
pet_files <- pet_files[order(get_year(pet_files))]
et_files  <- et_files [order(get_year(et_files))]
npp_files <- npp_files[order(get_year(npp_files))]

# ---- read rasters ----
pet_raw <- rast(pet_files)
et_raw  <- rast(et_files)
npp_raw <- rast(npp_files)

# ---- mask fill values (safe) ----
pet_raw[pet_raw >= 65535] <- NA
et_raw [et_raw  >= 65535] <- NA
npp_raw[npp_raw >= 32767] <- NA

# ---- scaling (same as before; verify with global() if you want) ----
pet <- pet_raw * 0.1
et  <- et_raw  * 0.1

# NPP: handle like we did before — DO NOT apply 0.0001 unless your ranges are 0–32766
# quick auto-check: if max > 100, assume unscaled ints and multiply by 0.0001
if (global(npp_raw[[1]], "max", na.rm=TRUE)[1,1] > 100) {
  npp <- npp_raw * 0.0001
} else {
  npp <- npp_raw
}

# name layers by year
pet_years <- sapply(pet_files, get_year); names(pet) <- paste0("PET_", pet_years)
et_years  <- sapply(et_files,  get_year); names(et)  <- paste0("ET_",  et_years)
npp_years <- sapply(npp_files, get_year); names(npp) <- paste0("NPP_", npp_years)

# ---- project watersheds to raster CRS ----
ws_pet <- st_transform(ws_sf, crs(pet))
ws_et  <- st_transform(ws_sf, crs(et))
ws_npp <- st_transform(ws_sf, crs(npp))

# ---- exact area-weighted means for each year ----
pet_df <- as.data.frame(exact_extract(pet, ws_pet, "mean"))
et_df  <- as.data.frame(exact_extract(et,  ws_et,  "mean"))
npp_df <- as.data.frame(exact_extract(npp, ws_npp, "mean"))

names(pet_df) <- names(pet)
names(et_df)  <- names(et)
names(npp_df) <- names(npp)

wide <- bind_cols(
  st_drop_geometry(ws_sf) |> dplyr::select(all_of(id_col)),
  pet_df, et_df, npp_df
)

# ---- annual long table (optional) ----
annual <- wide |>
  pivot_longer(cols = -all_of(id_col),
               names_to = "var_year",
               values_to = "value") |>
  separate(var_year, into = c("var", "year"), sep = "_", convert = TRUE) |>
  pivot_wider(names_from = var, values_from = value) |>
  arrange(!!sym(id_col), year)

write.csv(annual, out_annual_csv, row.names = FALSE)

# ---- sampling-window mean (weighted) ----
w2022 <- 9/19   # Apr–Dec
w2023 <- 10/19  # Jan–Oct

mean_tbl <- annual |>
  filter(year %in% c(2022, 2023)) |>
  group_by(!!sym(id_col)) |>
  summarise(
    PET_mean = w2022 * PET[year == 2022] + w2023 * PET[year == 2023],
    ET_mean  = w2022 * ET [year == 2022] + w2023 * ET [year == 2023],
    NPP_mean = w2022 * NPP[year == 2022] + w2023 * NPP[year == 2023],
    .groups = "drop"
  )

write.csv(mean_tbl, out_mean_csv, row.names = FALSE)

cat("Wrote annual: ", out_annual_csv, "\n")
cat("Wrote sampling-window mean: ", out_mean_csv, "\n")