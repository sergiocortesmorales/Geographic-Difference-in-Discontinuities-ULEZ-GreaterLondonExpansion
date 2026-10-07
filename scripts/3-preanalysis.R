# Install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, tidyverse, sf, rdrobust, patchwork, fixest)

# Read files
pp_onspd_south_sf <- readRDS(here("data", "pp_onspd_south_sf.rds"))
gla_boundary <- readRDS(here("data", "gla_boundary.rds"))

# Sales within 2km of the boundary, geometry dropped
pp_onspd_south_2km <- pp_onspd_south_sf %>%
  st_drop_geometry() %>%
  filter(abs(signed_dist) <= 2000)

# Plot settings
ulez_year_pos <- 2023.66   # 29 Aug 2023 as a decimal year
font_family   <- "serif"
theme_set(theme_minimal(base_family = font_family))

theme_eventstudy <- theme_minimal(base_size = 11, base_family = font_family) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", colour = "grey85"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

theme_descriptive <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

col_in_out <- c("Outside ULEZ" = "darkorange", "Inside ULEZ" = "darkgreen")

# Save a figure to output/ and print it
dir.create(here("output"), showWarnings = FALSE)
save_fig <- function(plot, name, width = 9, height = 5) {
  ggsave(here("output", paste0(name, ".png")), plot = plot, width = width, height = height, dpi = 400)
  print(plot)
  invisible(NULL)
}

# 1. Transaction volumes, inside vs outside (annual mean of quarterly counts) ----

local({
  vol_by_year <- pp_onspd_south_2km %>%
    count(quarter, inside_ulez) %>%
    mutate(year = as.integer(substr(quarter, 1, 4))) %>%
    group_by(year, inside_ulez) %>%
    summarise(n = mean(n), .groups = "drop")
  
  save_fig(
    (ggplot(vol_by_year, aes(year, n, colour = inside_ulez, group = inside_ulez)) +
       geom_line(linewidth = 0.8) +
       geom_vline(xintercept = ulez_year_pos, linetype = "dashed", colour = "red") +
       scale_colour_manual(values = col_in_out) +
       scale_x_continuous(breaks = 2013:2025) +
       labs(title = "Levels", x = NULL, y = "Mean quarterly transactions", colour = NULL) +
       theme_descriptive) +
      (vol_by_year %>%
         group_by(inside_ulez) %>%
         mutate(index = n / n[year == 2019] * 100) %>%
         ggplot(aes(year, index, colour = inside_ulez, group = inside_ulez)) +
         geom_line(linewidth = 0.8) +
         geom_vline(xintercept = ulez_year_pos, linetype = "dashed", colour = "red") +
         geom_hline(yintercept = 100, linetype = "dotted", colour = "grey50") +
         scale_colour_manual(values = col_in_out) +
         scale_x_continuous(breaks = 2013:2025) +
         labs(title = "Indexed (2019 = 100)", x = NULL, y = "Volume index", colour = NULL) +
         theme_descriptive) +
      plot_layout(guides = "collect"),
    "01_volumes", width = 10, height = 4
  )
})

# 2. Median prices, inside vs outside (annual mean of quarterly medians) ----

local({
  price_by_year <- pp_onspd_south_2km %>%
    group_by(quarter, inside_ulez) %>%
    summarise(median_price = median(price), .groups = "drop") %>%
    mutate(year = as.integer(substr(quarter, 1, 4))) %>%
    group_by(year, inside_ulez) %>%
    summarise(median_price = mean(median_price), .groups = "drop")
  
  save_fig(
    (ggplot(price_by_year, aes(year, median_price / 1000, colour = inside_ulez, group = inside_ulez)) +
       geom_line(linewidth = 0.8) +
       geom_vline(xintercept = ulez_year_pos, linetype = "dashed", colour = "red") +
       scale_colour_manual(values = col_in_out) +
       scale_x_continuous(breaks = 2013:2025) +
       labs(title = "Levels", x = NULL, y = "Median price (£000s)", colour = NULL) +
       theme_descriptive) +
      (price_by_year %>%
         group_by(inside_ulez) %>%
         mutate(index = median_price / median_price[year == 2019] * 100) %>%
         ggplot(aes(year, index, colour = inside_ulez, group = inside_ulez)) +
         geom_line(linewidth = 0.8) +
         geom_vline(xintercept = ulez_year_pos, linetype = "dashed", colour = "red") +
         geom_hline(yintercept = 100, linetype = "dotted", colour = "grey50") +
         scale_colour_manual(values = col_in_out) +
         scale_x_continuous(breaks = 2013:2025) +
         labs(title = "Indexed (2019 = 100)", x = NULL, y = "Price index", colour = NULL) +
         theme_descriptive) +
      plot_layout(guides = "collect"),
    "02_prices", width = 10, height = 4
  )
})

# 3. RD plot of log price, pre vs post ----

local({
  prepost_p <- 1
  
  # rdplot bins and fitted line for one period
  get_rdplot_layers <- function(period_label) {
    d <- filter(pp_onspd_south_2km, post_ulez == period_label)
    
    rp <- rdplot(y = d$log_price, x = d$signed_dist, c = 0, p = prepost_p,
                 nbins = c(20, 20), masspoints = "adjust", hide = TRUE)
    # Split the fitted line at the cutoff
    cut_row <- which(rp$vars_poly$rdplot_x >= 0)[1]
    
    pf <- factor(period_label, levels = c("Pre", "Post"))
    list(
      bins = tibble(x = rp$vars_bins$rdplot_mean_x, y = rp$vars_bins$rdplot_mean_y, post_ulez = pf),
      poly = tibble(x = rp$vars_poly$rdplot_x, y = rp$vars_poly$rdplot_y,
                    side = c(rep("Outside ULEZ", cut_row),
                             rep("Inside ULEZ", nrow(rp$vars_poly) - cut_row)),
                    post_ulez = pf)
    )
  }
  
  rdplot_layers <- map(c("Pre", "Post"), get_rdplot_layers)
  
  save_fig(
    ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
      geom_point(data = map_dfr(rdplot_layers, "bins"), aes(x / 1000, y),
                 alpha = 0.6, colour = "grey30", size = 1.2) +
      geom_line(data = map_dfr(rdplot_layers, "poly"), aes(x / 1000, y, colour = side), linewidth = 0.8) +
      facet_wrap(~ post_ulez) +
      scale_colour_manual(values = col_in_out) +
      labs(x = "Distance to boundary (km, positive = inside ULEZ)", y = "Mean log price", colour = NULL) +
      theme(legend.position = "bottom"),
    "03_rd_prepost", width = 9, height = 4.5
  )
})

# 4. Event study (RD): boundary gap in log price by period ----

local({
  event_study_h     <- 2000
  event_study_start <- 2013
  
  event_study_data <- pp_onspd_south_sf %>%
    st_drop_geometry() %>%
    filter(abs(signed_dist) <= 10000) %>%
    mutate(event_period = case_when(
      year < 2023                        ~ as.character(year),
      year == 2023 & post_ulez == "Pre"  ~ "2023-Pre",
      year == 2023 & post_ulez == "Post" ~ "2023-Post",
      year > 2023                        ~ as.character(year)
    ))
  
  # RD per period (h = 2km), NA if under 50 sales on either side
  event_study_results <- map_dfr(sort(unique(event_study_data$event_period)), function(pd) {
    d <- filter(event_study_data, event_period == pd)
    na_row <- tibble(event_period = pd, estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
    if (sum(d$signed_dist > 0) < 50 | sum(d$signed_dist < 0) < 50) return(na_row)
    
    rd_fit <- tryCatch(
      rdrobust(y = d$log_price, x = d$signed_dist, c = 0, h = event_study_h,
               cluster = d$lsoa21cd, masspoints = "adjust"),
      error = function(e) NULL
    )
    if (is.null(rd_fit)) return(na_row)
    
    tibble(event_period = pd,
           estimate = rd_fit$coef["Bias-Corrected", ],   # bias-corrected point, robust CI
           ci_lower = rd_fit$ci["Robust", "CI Lower"],
           ci_upper = rd_fit$ci["Robust", "CI Upper"])
  }) %>%
    left_join(tibble(event_period = c(as.character(2013:2022), "2023-Pre", "2023-Post", as.character(2024:2025)),
                     x_pos = c(2013:2022, 2023.0, 2023.5, 2024:2025)),
              by = "event_period") %>%
    filter(x_pos >= event_study_start)
  
  # Mean pre-period gap (dotted line)
  pre_gamma_mean <- mean(event_study_results$estimate[event_study_results$x_pos < 2023.25], na.rm = TRUE)
  
  # y range from all periods except 2023-Post
  es_ylim <- with(filter(event_study_results, event_period != "2023-Post"),
                  range(c(ci_lower, ci_upper), na.rm = TRUE))
  es_ylim <- es_ylim + c(-1, 1) * 0.04 * diff(es_ylim)
  
  save_fig(
    ggplot(event_study_results, aes(x_pos, estimate)) +
      geom_hline(yintercept = pre_gamma_mean, linetype = "dotted", colour = "grey50") +
      geom_vline(xintercept = 2023.25, linetype = "dashed", colour = "red") +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), colour = "black", width = 0.3, linewidth = 0.6) +
      geom_point(colour = "black", size = 2.2) +
      scale_x_continuous(breaks = event_study_results$x_pos, labels = event_study_results$event_period) +
      coord_cartesian(ylim = es_ylim) +
      labs(x = NULL, y = "RD gap in log price") +
      theme_eventstudy,
    "04_event_study_rd", width = 9, height = 4.5
  )
})

# 5. Event study (DiD): inside x year, postcode and year FE, ref 2022 ----

local({
  did_event_reg <- feols(
    log_price ~ i(year, inside_ulez_dummy, ref = "2022") | pcds + year,
    data = pp_onspd_south_2km %>%
      mutate(inside_ulez_dummy = as.numeric(inside_ulez == "Inside ULEZ"),
             year = as.factor(year)),
    cluster = "lsoa21cd"
  )
  print(summary(did_event_reg))
  
  # Coefficients with 95% CI, 2022 added as 0
  ct <- coeftable(did_event_reg)
  did_event_coefs <- tibble(term = rownames(ct), estimate = ct[, 1], se = ct[, 2]) %>%
    filter(str_detect(term, "inside_ulez_dummy")) %>%
    transmute(year = as.integer(str_extract(term, "[0-9]{4}")),
              estimate, conf.low = estimate - 1.96 * se, conf.high = estimate + 1.96 * se) %>%
    bind_rows(tibble(year = 2022, estimate = 0, conf.low = 0, conf.high = 0)) %>%
    arrange(year)
  
  save_fig(
    ggplot(did_event_coefs, aes(year, estimate)) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
      geom_vline(xintercept = 2022.5, linetype = "dashed", colour = "red") +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), colour = "black", width = 0.3, linewidth = 0.6) +
      geom_point(colour = "black", size = 2.2) +
      annotate("text", x = 2022, y = 0, label = "ref", vjust = -0.9, size = 3, colour = "grey40", family = font_family) +
      scale_x_continuous(breaks = sort(unique(did_event_coefs$year))) +
      labs(x = "Year", y = "Inside × Year coefficient (log price)") +
      theme_eventstudy,
    "05_event_study_did", width = 8, height = 4.5
  )
})

# 6. Transaction density near the boundary, pre vs post ----

local({
  # Post minus pre share of sales by 250m bin
  density_change <- pp_onspd_south_2km %>%
    mutate(dist_bin = floor(signed_dist / 250) * 250 + 125) %>%
    count(post_ulez, dist_bin) %>%
    group_by(post_ulez) %>%
    mutate(share = n / sum(n)) %>%
    ungroup() %>%
    select(post_ulez, dist_bin, share) %>%
    pivot_wider(names_from = post_ulez, values_from = share, values_fill = 0) %>%
    mutate(share_diff = Post - Pre)
  
  save_fig(
    (ggplot(filter(pp_onspd_south_2km, post_ulez == "Pre"), aes(signed_dist / 1000)) +
       geom_histogram(binwidth = 0.25, boundary = 0, fill = "grey20", colour = "white") +
       geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
       labs(title = "Pre (before 29 Aug 2023)", x = "Distance to boundary (km)", y = "Count")) +
      (ggplot(filter(pp_onspd_south_2km, post_ulez == "Post"), aes(signed_dist / 1000)) +
         geom_histogram(binwidth = 0.25, boundary = 0, fill = "grey20", colour = "white") +
         geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
         labs(title = "Post (from 29 Aug 2023)", x = "Distance to boundary (km)", y = "Count")) +
      (ggplot(density_change, aes(dist_bin / 1000, share_diff)) +
         geom_col(fill = "grey20", colour = "white") +
         geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
         geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
         labs(title = "Change in transaction share (post minus pre)",
              x = "Distance to boundary (km)", y = "Δ share of transactions")) +
      plot_layout(ncol = 3),
    "06_density", width = 12, height = 4
  )
})

# 7. Map of sales within 2km of the boundary ----

save_fig(
  ggplot() +
    geom_sf(data = gla_boundary, fill = NA, colour = "black", linewidth = 0.5) +
    geom_sf(data = pp_onspd_south_sf %>% filter(abs(signed_dist) <= 2000),
            aes(colour = inside_ulez), size = 0.2) +
    scale_colour_manual(values = col_in_out) +
    labs(colour = NULL),
  "07_map_transactions", width = 7, height = 6
)

# 8. Maps of NO2 and PM2.5 change, 2022 to 2024 ----
# purple = reduction, orange = increase
local({
  pollution_map_sample <- pp_onspd_south_sf %>%
    filter(abs(signed_dist) <= 2000, !is.na(delta_no2), !is.na(delta_pm25)) %>%
    mutate(pct_no2  = 100 * delta_no2  / no2_2022,
           pct_pm25 = 100 * delta_pm25 / pm25_2022)
  
  # Shaded inside and outside areas, clipped to the map extent
  frame_poly   <- st_as_sfc(st_bbox(pollution_map_sample))
  inside_area  <- st_intersection(st_geometry(gla_boundary), frame_poly)
  outside_area <- st_difference(frame_poly, st_geometry(gla_boundary))
  
  pollution_map <- function(pct_var, title) {
    ggplot() +
      geom_sf(data = outside_area, fill = "darkorange", colour = NA, alpha = 0.15) +
      geom_sf(data = inside_area,  fill = "darkgreen",  colour = NA, alpha = 0.15) +
      geom_sf(data = gla_boundary, fill = NA, colour = "black", linewidth = 0.5) +
      geom_sf(data = pollution_map_sample, aes(colour = .data[[pct_var]]), size = 0.2) +
      scale_colour_distiller(palette = "PuOr", direction = -1,
                             limits = c(-20, 20), oob = scales::squish,
                             name = "Change (%)") +
      coord_sf(xlim = as.numeric(st_bbox(frame_poly)[c("xmin", "xmax")]),
               ylim = as.numeric(st_bbox(frame_poly)[c("ymin", "ymax")]), expand = FALSE)
  }
  
  save_fig(pollution_map("pct_no2",  "NO2 change (2022 to 2024)"),  "08_map_no2_change",  width = 7, height = 6)
  save_fig(pollution_map("pct_pm25", "PM2.5 change (2022 to 2024)"), "08_map_pm25_change", width = 7, height = 6)
})
# 9. Zoomed pollution maps at one boundary section ----
local({
  zoom_half     <- 3000
  zoom_coverage <- 2000   # max distance from the boundary to plot
  # Frame centred on the 2nd busiest 2km cell within 250m of the boundary
  xy   <- st_coordinates(pp_onspd_south_sf)
  cell <- as_tibble(xy)[abs(pp_onspd_south_sf$signed_dist) <= 250, ] %>%
    transmute(gx = round(X / 2000) * 2000, gy = round(Y / 2000) * 2000) %>%
    count(gx, gy, sort = TRUE) %>%
    slice(2)
  bbox <- st_bbox(c(xmin = cell$gx - zoom_half, xmax = cell$gx + zoom_half,
                    ymin = cell$gy - zoom_half, ymax = cell$gy + zoom_half),
                  crs = st_crs(pp_onspd_south_sf))
  in_frame <- xy[, 1] >= bbox["xmin"] & xy[, 1] <= bbox["xmax"] &
    xy[, 2] >= bbox["ymin"] & xy[, 2] <= bbox["ymax"] &
    abs(pp_onspd_south_sf$signed_dist) <= zoom_coverage
  zoom_data  <- pp_onspd_south_sf[in_frame, ] %>%
    mutate(pct_no2  = 100 * delta_no2  / no2_2022,
           pct_pm25 = 100 * delta_pm25 / pm25_2022)
  zoom_line  <- st_crop(st_boundary(gla_boundary), bbox)
  zoom_coord <- coord_sf(xlim = as.numeric(bbox[c("xmin", "xmax")]),
                         ylim = as.numeric(bbox[c("ymin", "ymax")]), expand = FALSE)
  # Shaded inside and outside areas
  frame_poly   <- st_as_sfc(bbox)
  inside_zoom  <- st_intersection(st_geometry(gla_boundary), frame_poly)
  outside_zoom <- st_difference(frame_poly, st_geometry(gla_boundary))
  pollution_zoom <- function(pct_var, title) {
    ggplot() +
      geom_sf(data = outside_zoom, fill = "darkorange", colour = NA, alpha = 0.15) +
      geom_sf(data = inside_zoom,  fill = "darkgreen",  colour = NA, alpha = 0.15) +
      geom_sf(data = zoom_line, colour = "black", linewidth = 0.5) +
      geom_sf(data = zoom_data, aes(fill = .data[[pct_var]]), shape = 21,
              colour = "grey20", stroke = 0.2, size = 1.7) +
      scale_fill_distiller(palette = "PuOr", direction = -1,
                           limits = c(-20, 20), oob = scales::squish,
                           name = "Change (%)") +
      zoom_coord
  }
  save_fig(pollution_zoom("pct_no2",  "NO2 change (2km coverage)"),  "08b_map_no2_change_zoom",  width = 6, height = 6)
  save_fig(pollution_zoom("pct_pm25", "PM2.5 change (2km coverage)"), "08b_map_pm25_change_zoom", width = 6, height = 6)
})