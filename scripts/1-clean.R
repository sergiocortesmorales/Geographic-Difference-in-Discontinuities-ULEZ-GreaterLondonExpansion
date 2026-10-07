# Install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, tidyverse, readODS)

# 1. Land Registry Price Paid Data ----
pp_full <- read.csv(here("rawdata", "pp-complete.csv"), header = FALSE)

# Sales from 2013 onwards
pp_2013to2026 <- pp_full %>%
  filter(V3 >="2013")
rm(pp_full)
pp_2013to2026 <- pp_2013to2026 %>%
  rename(
    transaction_id = V1,
    price = V2,
    date = V3,
    pcds = V4,
    property_type = V5,    # D=Detached, S=Semi, T=Terraced, F=Flat, O=Other
    new_build = V6,        # Y/N
    tenure = V7,           # F=Freehold, L=Leasehold
    paon = V8,             # house number or name
    saon = V9,             # flat or unit number
    street = V10,
    locality = V11,
    town = V12,
    district = V13,
    county = V14,
    ppd_cat = V15,         # A=standard, B=additional (repossessions, transfers)
    record_status = V16    # A=Addition, C=Change, D=Deletion
  )
saveRDS(pp_2013to2026, here("data", "pp_2013to2026.rds"))
# pp_2013to2026 <- readRDS(here("data", "pp_2013to2026.rds"))

# Restrict to 2013-2025
pp_2013to2025 <- pp_2013to2026 %>% 
  filter(date < "2026")
saveRDS(pp_2013to2025, here("data", "pp_2013to2025.rds"))
# pp_2013to2025 <- readRDS(here("data", "pp_2013to2025.rds"))

# 2. ONS Postcode Directory ----
onspd_full <- read.csv(here("rawdata", "ONSPD_NOV_2025_UK.csv"))
# Coordinates, postcode dates, area codes and IMD
onspd <- onspd_full %>%
  select(pcds, lat, long, gridind, dointr, doterm, lsoa11cd, lsoa21cd, msoa21cd, lad25cd, rgn25cd, imd20ind)
rm(onspd_full)

# 3. DfT Journey Time Statistics, JTS0501 (2017) ----
# Travel times to employment centres by mode, 2011 LSOAs
# Source: https://www.gov.uk/government/statistical-data-sets/journey-time-statistics-data-tables-jts
# ods edited by hand first: variable names in row 1, 2017 sheet only

jts0501_full <- read_ods(here("rawdata", "jts0501.ods"), sheet = 1)

# Average minimum time (mins) to the nearest centre with 100+, 500+, 5000+ jobs
jts0501 <- jts0501_full %>%
  select(
    lsoa11cd       = LSOA_code,
    # public transport
    jts_pt_100emp_mins  = `100EmpPTt`,
    jts_pt_500emp_mins  = `500EmpPTt`,
    jts_pt_5000emp_mins = `5000EmpPTt`,
    # car
    jts_car_100emp_mins  = `100EmpCart`,
    jts_car_500emp_mins  = `500EmpCart`,
    jts_car_5000emp_mins = `5000EmpCart`,
    # cycle
    jts_cyc_100emp_mins  = `100EmpCyct`,
    jts_cyc_500emp_mins  = `500EmpCyct`,
    jts_cyc_5000emp_mins = `5000EmpCyct`
  ) %>%
  mutate(across(-lsoa11cd, ~ as.numeric(.x)))

cat("JTS0501 rows:", nrow(jts0501), "\n")
cat("Non-missing rows:", sum(complete.cases(jts0501)), "of", nrow(jts0501), "\n")

rm(jts0501_full)

# 4. Merge sales with postcode data ----
pp_onspd <- pp_2013to2025 %>%
  left_join(onspd, by = "pcds")

# Drop category B sales (not at full market value)
pp_onspd <- pp_onspd %>%
  filter(
    ppd_cat != "B",
  )

# Sales with no coordinates
cat("Total rows:", nrow(pp_onspd), "\n")
cat("Missing lat:", sum(is.na(pp_onspd$lat)), "\n")
cat("% missing:", round(100 * mean(is.na(pp_onspd$lat)), 2), "%\n")

# Missing coordinates by year, property type and county
pp_onspd %>%
  filter(is.na(lat)) %>%
  mutate(year = substr(date, 1, 4)) %>%
  count(year) # no clustering
pp_onspd %>%
  filter(is.na(lat)) %>%
  count(property_type) %>%
  arrange(desc(n)) # no clustering
pp_onspd %>%
  filter(is.na(lat)) %>%
  count(county) %>%
  arrange(desc(n)) # London first, as expected

# Check unmatched postcodes for formatting issues
unmatched <- pp_onspd %>%
  filter(is.na(lat)) %>%
  pull(pcds) %>%
  unique()
cat("Unique unmatched postcodes:", length(unmatched), "\n")
head(unmatched, 30)
cat("Empty strings:", sum(unmatched == ""), "\n")
cat("Extra spaces:", sum(str_detect(unmatched, "  ")), "\n")
cat("Leading/trailing space:", sum(unmatched != str_trim(unmatched)), "\n")
cat("Unusual length:", sum(nchar(unmatched) < 5 | nchar(unmatched) > 8), "\n")
# Look up unmatched postcodes in ONSPD (terminated ones are kept on file)
onspd %>%
  filter(pcds %in% unmatched[unmatched != ""]) %>%
  select(pcds, gridind, lat, long, doterm) # 0 rows
# Same check with spaces removed
sample_unmatched <- head(unmatched[unmatched != ""], 10)
onspd %>%
  filter(str_remove_all(pcds, " ") %in% str_remove_all(sample_unmatched, " ")) %>%
  select(pcds, gridind, lat, long, doterm) # 0 rows
# Same check against the raw ONSPD file
onspd_check <- read.csv(here("rawdata", "ONSPD_NOV_2025_UK.csv")) %>%
  filter(pcds %in% unmatched[unmatched != ""]) %>%
  select(pcds, gridind, lat, long, doterm)
cat("Unmatched postcodes found in raw ONSPD:", nrow(onspd_check), "\n")
print(onspd_check) # 0 rows, postcodes not in ONSPD
rm(onspd_check)
cat("% missing:", round(100 * mean(is.na(pp_onspd$lat)), 2), "%\n")
# About 0.02% of sales, dropped
pp_onspd <- pp_onspd %>%
  filter(!is.na(lat))

# 5. Southern England regions ----
pp_onspd_south <- pp_onspd %>%
  filter(rgn25cd %in% c(
    "E12000007",  # London
    "E12000008",  # South East
    "E12000006",  # East of England
    "E12000009"   # South West
  ))

# 6. Add JTS travel times (joined on 2011 LSOA) ----
pp_onspd_south <- pp_onspd_south %>%
  left_join(jts0501, by = "lsoa11cd")

jts_merge_rate <- mean(!is.na(pp_onspd_south$jts_pt_500emp_mins))
cat("JTS merge rate:", round(100 * jts_merge_rate, 1), "%\n")
cat("Missing JTS:", sum(is.na(pp_onspd_south$jts_pt_500emp_mins)),
    "of", nrow(pp_onspd_south), "\n")

# Save full merge and clear memory
saveRDS(pp_onspd, here("data", "pp_onspd.rds"))
rm(onspd, pp_onspd, pp_2013to2025, pp_2013to2026, sample_unmatched, unmatched, jts0501)
# Save southern England subset
saveRDS(pp_onspd_south, here("data", "pp_onspd_south.rds"))