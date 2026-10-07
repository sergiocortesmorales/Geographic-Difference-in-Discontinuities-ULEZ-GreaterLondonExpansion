# Install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, tidyverse, sf)

# Read PPD-ONSPD merge
pp_onspd_south <- readRDS(here("data", "pp_onspd_south.rds"))

# 1. Boundary and distance ----

# Greater London boundary (edge of the London-wide ULEZ), dissolved, British National Grid
gla_boundary <- st_read(here("data", "Greater_London_Authority_(GLA).shp")) %>%
  st_make_valid() %>%
  st_union() %>%
  st_sf() %>%
  st_transform(27700)

gla_line <- st_boundary(gla_boundary)
plot(st_geometry(gla_boundary))

# Sales to spatial points (BNG)
pp_onspd_south_sf <- pp_onspd_south %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(27700)
rm(pp_onspd_south)

# 10km buffer around the boundary
gla_buffer_10km <- st_buffer(gla_boundary, dist = 10000)

# Keep sales within the buffer
near_boundary <- st_intersects(pp_onspd_south_sf, gla_buffer_10km, sparse = FALSE)[, 1]
pp_onspd_south_sf <- pp_onspd_south_sf[near_boundary, ]
rm(near_boundary, gla_buffer_10km)

# Distance to the boundary (m)
pp_onspd_south_sf$dist_to_boundary <- as.numeric(st_distance(pp_onspd_south_sf, gla_line))

# Signed distance: positive inside Greater London, negative outside
pp_onspd_south_sf$signed_dist <- ifelse(
  st_intersects(pp_onspd_south_sf, gla_boundary, sparse = FALSE)[, 1],
  pp_onspd_south_sf$dist_to_boundary,
  -pp_onspd_south_sf$dist_to_boundary
)
# Unsigned distance no longer needed
pp_onspd_south_sf$dist_to_boundary <- NULL

summary(pp_onspd_south_sf$signed_dist)
cat("Properties inside GLA:", sum(pp_onspd_south_sf$signed_dist > 0), "\n")
cat("Properties outside GLA:", sum(pp_onspd_south_sf$signed_dist < 0), "\n")

# Date, treatment and log price variables
pp_onspd_south_sf <- pp_onspd_south_sf %>%
  mutate(
    # Date
    date_full    = as.Date(date),
    year         = as.numeric(format(date_full, "%Y")),
    month        = format(date_full, "%Y-%m"),
    quarter      = paste0(format(date_full, "%Y"), "-Q", ceiling(as.numeric(format(date_full, "%m")) / 3)),
    day          = date_full,
    # Treatment: inside Greater London, sold from 29 Aug 2023
    inside_ulez  = factor(
      ifelse(signed_dist > 0, "Inside ULEZ", "Outside ULEZ"),
      levels = c("Outside ULEZ", "Inside ULEZ")
    ),
    post_ulez    = factor(
      ifelse(date_full >= as.Date("2023-08-29"), "Post", "Pre"),
      levels = c("Pre", "Post")
    ),
    # Log price
    log_price    = log(price)
  )

# 2. Defra PCM air quality ----
# Modelled annual mean NO2 and PM2.5 (ug/m3), 1km grid, BNG
# Source: https://uk-air.defra.gov.uk/data/pcm-data
# csvs edited by hand first: metadata rows removed, 4 columns kept
# 2022 = pre-ULEZ, 2024 = first full post-ULEZ year

# NO2
pcm_no2_2022 <- read.csv(here("rawdata", "mapno22022.csv")) %>%
  rename(no2_2022 = 4) %>%
  mutate(no2_2022 = as.numeric(no2_2022))

pcm_no2_2024 <- read.csv(here("rawdata", "mapno22024.csv")) %>%
  rename(no2_2024 = 4) %>%
  mutate(no2_2024 = as.numeric(no2_2024))

pcm_no2 <- pcm_no2_2022 %>%
  select(gridcode, x, y, no2_2022) %>%
  left_join(pcm_no2_2024 %>% select(gridcode, no2_2024), by = "gridcode")

rm(pcm_no2_2022, pcm_no2_2024)

# PM2.5
pcm_pm25_2022 <- read.csv(here("rawdata", "mappm252022g.csv")) %>%
  rename(pm25_2022 = 4) %>%
  mutate(pm25_2022 = as.numeric(pm25_2022))

pcm_pm25_2024 <- read.csv(here("rawdata", "mappm252024g.csv")) %>%
  rename(pm25_2024 = 4) %>%
  mutate(pm25_2024 = as.numeric(pm25_2024))

pcm_pm25 <- pcm_pm25_2022 %>%
  select(gridcode, x, y, pm25_2022) %>%
  left_join(pcm_pm25_2024 %>% select(gridcode, pm25_2024), by = "gridcode")

rm(pcm_pm25_2022, pcm_pm25_2024)

# Combine into one grid
pcm_grid <- pcm_no2 %>%
  left_join(pcm_pm25 %>% select(gridcode, pm25_2022, pm25_2024), by = "gridcode")

rm(pcm_no2, pcm_pm25)

cat("PCM grid cells:", nrow(pcm_grid), "\n")
cat("With both NO2 and PM2.5:", sum(complete.cases(pcm_grid)), "\n")

# Grid cell centres to points
pcm_grid_sf <- pcm_grid %>%
  st_as_sf(coords = c("x", "y"), crs = 27700)

# Nearest grid cell for each sale
nearest_pcm_idx <- st_nearest_feature(pp_onspd_south_sf, pcm_grid_sf)

# Add pollution values and 2022-2024 change
pp_onspd_south_sf <- pp_onspd_south_sf %>%
  mutate(
    pcm_cell_id = pcm_grid$gridcode[nearest_pcm_idx],
    no2_2022  = pcm_grid$no2_2022[nearest_pcm_idx],
    no2_2024  = pcm_grid$no2_2024[nearest_pcm_idx],
    pm25_2022 = pcm_grid$pm25_2022[nearest_pcm_idx],
    pm25_2024 = pcm_grid$pm25_2024[nearest_pcm_idx],
    # 2024 minus 2022, negative = improvement
    delta_no2  = no2_2024 - no2_2022,
    delta_pm25 = pm25_2024 - pm25_2022
  )

rm(pcm_grid, pcm_grid_sf, nearest_pcm_idx)

# Merge checks
cat("\n=== PCM merge diagnostics ===\n")
cat("Missing NO2:", sum(is.na(pp_onspd_south_sf$no2_2022)), "of", nrow(pp_onspd_south_sf), "\n")
cat("Missing PM2.5:", sum(is.na(pp_onspd_south_sf$pm25_2022)), "of", nrow(pp_onspd_south_sf), "\n")
cat("\nNO2 2022 summary (ug/m3):\n"); print(summary(pp_onspd_south_sf$no2_2022))
cat("NO2 2024 summary (ug/m3):\n"); print(summary(pp_onspd_south_sf$no2_2024))
cat("NO2 change summary:\n"); print(summary(pp_onspd_south_sf$delta_no2))
cat("\nPM2.5 2022 summary (ug/m3):\n"); print(summary(pp_onspd_south_sf$pm25_2022))
cat("PM2.5 2024 summary (ug/m3):\n"); print(summary(pp_onspd_south_sf$pm25_2024))
cat("PM2.5 change summary:\n"); print(summary(pp_onspd_south_sf$delta_pm25))

# 3. NaPTAN rail and metro stations ----

stops_raw <- read.csv(here("rawdata", "Stops.csv"))

# Active rail and metro access points with valid coordinates
stations <- stops_raw %>%
  filter(Status == "active",
         StopType %in% c("RLY", "MET", "RSE", "TMU", "PLT")) %>%
  mutate(Easting = as.numeric(Easting), Northing = as.numeric(Northing)) %>%
  filter(!is.na(Easting), !is.na(Northing), Easting > 0, Northing > 0) %>%
  select(ATCOCode, CommonName, StopType, Easting, Northing)

rm(stops_raw)

stations_sf <- stations %>%
  st_as_sf(coords = c("Easting", "Northing"), crs = 27700)

saveRDS(stations_sf, here("data", "stations_sf.rds"))
cat("\nStation points (active rail/metro):", nrow(stations_sf), "\n")

# Distance to nearest station (m)
nearest_station_idx <- st_nearest_feature(pp_onspd_south_sf, stations_sf)
pp_onspd_south_sf$dist_to_station <- as.numeric(
  st_distance(pp_onspd_south_sf, stations_sf[nearest_station_idx, ], by_element = TRUE)
)

rm(stations, stations_sf, nearest_station_idx)

cat("Distance to nearest station summary (m):\n"); print(summary(pp_onspd_south_sf$dist_to_station))
cat("Within 800m of a station:",  sum(pp_onspd_south_sf$dist_to_station <= 800),  "\n")
cat("Within 1600m of a station:", sum(pp_onspd_south_sf$dist_to_station <= 1600), "\n")

# Save sales and boundary objects
saveRDS(pp_onspd_south_sf, here("data", "pp_onspd_south_sf.rds"))
saveRDS(gla_boundary, here("data", "gla_boundary.rds"))
saveRDS(gla_line, here("data", "gla_line.rds"))

# Clean memory
rm(pp_onspd_south_sf, gla_boundary, gla_line)