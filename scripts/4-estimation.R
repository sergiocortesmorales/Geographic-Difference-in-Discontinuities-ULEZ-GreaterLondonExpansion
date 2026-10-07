# Install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, tidyverse, sf, rdrobust, fixest, flextable, officer)

# Read files
pp_onspd_south_sf <- readRDS(here("data", "pp_onspd_south_sf.rds"))

dir.create(here("output"), showWarnings = FALSE)

# Plot font
font_family <- "serif"

# 1. Repeat-sales panel ----

# Property id from postcode and address
pp_onspd_south_sf$property_id <- paste(
  pp_onspd_south_sf$pcds, pp_onspd_south_sf$paon,
  pp_onspd_south_sf$saon, pp_onspd_south_sf$street, sep = "|"
)

# Last pre-ULEZ sale per property
pp_onspd_south_sf_pre <- pp_onspd_south_sf %>%
  st_drop_geometry() %>%
  filter(post_ulez == "Pre") %>%
  arrange(property_id, desc(date_full)) %>%
  distinct(property_id, .keep_all = TRUE) %>%
  select(property_id, signed_dist, log_price_pre = log_price, date_pre = date_full,
         property_type, tenure, new_build, lsoa21cd, pcm_cell_id, dist_to_station,
         delta_no2, delta_pm25)

# First post-ULEZ sale per property
pp_onspd_south_sf_post <- pp_onspd_south_sf %>%
  st_drop_geometry() %>%
  filter(post_ulez == "Post") %>%
  arrange(property_id, date_full) %>%
  distinct(property_id, .keep_all = TRUE) %>%
  select(property_id, log_price_post = log_price, date_post = date_full)

# Properties sold in both periods
pp_onspd_south_sf_repeatsales <- pp_onspd_south_sf_pre %>%
  inner_join(pp_onspd_south_sf_post, by = "property_id") %>%
  mutate(
    delta_log_price = log_price_post - log_price_pre,
    years_between   = as.numeric(date_post - date_pre) / 365.25,
    year_pre  = as.numeric(format(date_pre, "%Y")),
    year_post = as.numeric(format(date_post, "%Y")),
    year_pair = paste(year_pre, year_post, sep = "_"),
    is_flat      = as.numeric(property_type == "F"),
    is_terraced  = as.numeric(property_type == "T"),
    is_semi      = as.numeric(property_type == "S"),
    is_leasehold = as.numeric(tenure == "L"),
    is_new_build = as.numeric(new_build == "Y")
  )

rm(pp_onspd_south_sf_pre, pp_onspd_south_sf_post)

# Year-pair dummies (first dropped) and property covariates
year_pair_dummies <- model.matrix(~ factor(year_pair) - 1, data = pp_onspd_south_sf_repeatsales)
year_pair_dummies <- year_pair_dummies[, -1]

property_covariates <- with(pp_onspd_south_sf_repeatsales,
                            cbind(is_flat, is_terraced, is_semi, is_leasehold, is_new_build))

all_covariates <- cbind(property_covariates, year_pair_dummies)

cat("Repeat sales properties:", nrow(pp_onspd_south_sf_repeatsales), "\n")
cat("Within 2km:", sum(abs(pp_onspd_south_sf_repeatsales$signed_dist) <= 2000), "\n")
cat("Within 1km:", sum(abs(pp_onspd_south_sf_repeatsales$signed_dist) <= 1000), "\n")

# 2. Helpers ----

# Bias-corrected estimate with robust SE, CI and effective N
extract_rd <- function(rd_obj, label = "") {
  tibble(
    specification   = label,
    estimate_bc     = rd_obj$coef["Bias-Corrected", ],
    robust_se       = rd_obj$se["Robust", ],
    robust_z        = rd_obj$z["Robust", ],
    robust_pval     = rd_obj$pv["Robust", ],
    ci_lower        = rd_obj$ci["Robust", "CI Lower"],
    ci_upper        = rd_obj$ci["Robust", "CI Upper"],
    bandwidth       = rd_obj$bws["h", "left"],
    n_eff_left      = rd_obj$N_h[1],
    n_eff_right     = rd_obj$N_h[2],
    n_eff_total     = sum(rd_obj$N_h)
  )
}

# Drop covariates that are constant within a subsample
drop_constant_cols <- function(covs_matrix) {
  varying <- apply(covs_matrix, 2, function(col) length(unique(col)) > 1)
  covs_matrix[, varying, drop = FALSE]
}

# 3. Main estimates: 4 bandwidths x 2 polynomial orders x 4 specifications ----

bandwidth_configs <- list(
  list(h = 500,  label = "Bandwidth = 500m"),
  list(h = 1000, label = "Bandwidth = 1000m"),
  list(h = 1500, label = "Bandwidth = 1500m"),
  list(h = 2000, label = "Bandwidth = 2000m")
)

# Donut samples
donut_100_idx <- abs(pp_onspd_south_sf_repeatsales$signed_dist) > 100
donut_200_idx <- abs(pp_onspd_south_sf_repeatsales$signed_dist) > 200

run_four_specs <- function(poly_order, h_val) {
  h_args <- list(h = h_val)
  
  rd_fe <- do.call(rdrobust, c(list(
    y = pp_onspd_south_sf_repeatsales$delta_log_price,
    x = pp_onspd_south_sf_repeatsales$signed_dist,
    c = 0, p = poly_order,
    covs = year_pair_dummies,
    cluster = pp_onspd_south_sf_repeatsales$lsoa21cd,
    masspoints = "adjust"), h_args))
  
  rd_prop <- do.call(rdrobust, c(list(
    y = pp_onspd_south_sf_repeatsales$delta_log_price,
    x = pp_onspd_south_sf_repeatsales$signed_dist,
    c = 0, p = poly_order,
    covs = all_covariates,
    cluster = pp_onspd_south_sf_repeatsales$lsoa21cd,
    masspoints = "adjust"), h_args))
  
  rd_d100 <- do.call(rdrobust, c(list(
    y = pp_onspd_south_sf_repeatsales$delta_log_price[donut_100_idx],
    x = pp_onspd_south_sf_repeatsales$signed_dist[donut_100_idx],
    c = 0, p = poly_order,
    covs = drop_constant_cols(all_covariates[donut_100_idx, ]),
    cluster = pp_onspd_south_sf_repeatsales$lsoa21cd[donut_100_idx],
    masspoints = "adjust"), h_args))
  
  rd_d200 <- do.call(rdrobust, c(list(
    y = pp_onspd_south_sf_repeatsales$delta_log_price[donut_200_idx],
    x = pp_onspd_south_sf_repeatsales$signed_dist[donut_200_idx],
    c = 0, p = poly_order,
    covs = drop_constant_cols(all_covariates[donut_200_idx, ]),
    cluster = pp_onspd_south_sf_repeatsales$lsoa21cd[donut_200_idx],
    masspoints = "adjust"), h_args))
  
  bind_rows(
    extract_rd(rd_fe,   "FE only"),
    extract_rd(rd_prop, "All covariates"),
    extract_rd(rd_d100, "Donut 100m"),
    extract_rd(rd_d200, "Donut 200m")
  )
}

all_results <- map_dfr(bandwidth_configs, function(cfg) {
  map_dfr(c(1, 2), function(p) {
    cat("Running:", cfg$label, "| p =", p, "\n")
    run_four_specs(poly_order = p, h_val = cfg$h) %>%
      mutate(bw_label = cfg$label, poly_order = p)
  })
})

cat("\n=== ALL RESULTS ===\n")
print(all_results, n = Inf)

# Preferred: h = 1000m, linear, all covariates
rd_preferred <- with(pp_onspd_south_sf_repeatsales,
                     rdrobust(y = delta_log_price, x = signed_dist, c = 0, h = 1000,
                              covs = all_covariates, cluster = lsoa21cd, masspoints = "adjust"))

# 4. Result tables: linear (main) and quadratic (appendix) ----

build_panel_table <- function(results_df, poly, table_notes,
                              bw_labels = c("Bandwidth = 500m", "Bandwidth = 1000m",
                                            "Bandwidth = 1500m", "Bandwidth = 2000m")) {
  df <- results_df %>% filter(poly_order == poly)
  
  col_labels <- c("FE only", "All covariates", "Donut 100m", "Donut 200m")
  
  rows <- list()
  
  for (bw_lab in bw_labels) {
    panel_df <- df %>% filter(bw_label == bw_lab)
    
    est_stars <- panel_df %>%
      mutate(
        stars = case_when(robust_pval < 0.01 ~ "***",
                          robust_pval < 0.05 ~ "**",
                          robust_pval < 0.10 ~ "*",
                          TRUE ~ ""),
        est_str = paste0(round(estimate_bc, 4), stars),
        se_str  = paste0("(", round(robust_se, 4), ")"),
        n_str   = format(n_eff_total, big.mark = ",")
      )
    
    has_prop     <- c("", "\u2713", "\u2713", "\u2713")
    has_yearpair <- c("\u2713", "\u2713", "\u2713", "\u2713")
    
    rows <- c(rows, list(
      c(bw_lab, rep("", 4)),
      c("ULEZ effect",         est_stars$est_str),
      c("",                    est_stars$se_str),
      c("Property covariates", has_prop),
      c("Year-pair FE",        has_yearpair),
      c("N",                   est_stars$n_str)
    ))
  }
  
  table_df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(table_df) <- c(" ", col_labels)
  
  panel_rows <- which(table_df[, 1] %in% bw_labels)
  
  ft <- flextable(table_df) %>%
    align(j = 2:5, align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "all") %>%
    bold(i = panel_rows, part = "body") %>%
    italic(i = panel_rows, part = "body") %>%
    border_remove() %>%
    hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) %>%
    hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
    add_footer_lines(table_notes) %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    autofit()
  
  for (pr in panel_rows) {
    ft <- hline(ft, i = pr, part = "body",
                border = fp_border(color = "grey70", width = 0.5))
  }
  
  ft
}

ft_table1 <- build_panel_table(
  all_results, poly = 1,
  table_notes = "Notes: Dependent variable is first-differenced log price (post minus pre). Bias-corrected point estimates with robust SE in parentheses (CCT 2014). LSOA-clustered. Donut columns exclude properties within the stated radius and include all covariates. *** p<0.01, ** p<0.05, * p<0.10.",
  bw_labels = c("Bandwidth = 500m", "Bandwidth = 1000m", "Bandwidth = 1500m")
)

ft_table_appendix <- build_panel_table(
  all_results, poly = 2,
  table_notes = "Notes: Dependent variable is first-differenced log price (post minus pre). Bias-corrected point estimates with robust SE in parentheses (CCT 2014). LSOA-clustered. Donut columns exclude properties within the stated radius and include all covariates. *** p<0.01, ** p<0.05, * p<0.10."
)

# Full sales data no longer needed
rm(pp_onspd_south_sf)

# 5. First-differenced RD plot ----

rd_main_sub <- filter(pp_onspd_south_sf_repeatsales, abs(signed_dist) <= 2000)

rp_main <- rdplot(y = rd_main_sub$delta_log_price, x = rd_main_sub$signed_dist,
                  c = 0, p = 1, nbins = c(20, 20), masspoints = "adjust", hide = TRUE)

# Split the fitted line at the cutoff
cut_row <- which(rp_main$vars_poly$rdplot_x >= 0)[1]

rd_main_bins <- tibble(x = rp_main$vars_bins$rdplot_mean_x,
                       y = rp_main$vars_bins$rdplot_mean_y)

rd_main_poly <- tibble(x = rp_main$vars_poly$rdplot_x,
                       y = rp_main$vars_poly$rdplot_y,
                       side = c(rep("Outside ULEZ", cut_row),
                                rep("Inside ULEZ", nrow(rp_main$vars_poly) - cut_row)))

rd_main_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  geom_point(data = rd_main_bins, aes(x / 1000, y), alpha = 0.6, colour = "grey30", size = 1.2) +
  geom_line(data = rd_main_poly, aes(x / 1000, y, colour = side), linewidth = 0.8) +
  scale_colour_manual(values = c("Outside ULEZ" = "darkorange", "Inside ULEZ" = "darkgreen")) +
  labs(x = "Distance to boundary (km, positive = inside ULEZ)",
       y = "Delta log price (post - pre)", colour = NULL) +
  theme_minimal(base_family = font_family) +
  theme(legend.position = "bottom")
rd_main_plot

# 6. Bandwidth sensitivity ----

bandwidths <- seq(500, 5000, by = 250)

bandwidth_results <- map_dfr(bandwidths, function(bw) {
  rd_fit <- tryCatch(
    with(pp_onspd_south_sf_repeatsales,
         rdrobust(y = delta_log_price, x = signed_dist, c = 0,
                  h = bw, covs = all_covariates,
                  cluster = lsoa21cd, masspoints = "adjust")),
    error = function(e) NULL
  )
  if (is.null(rd_fit)) {
    return(tibble(bandwidth_m = bw, estimate = NA_real_, ci_lower = NA_real_,
                  ci_upper = NA_real_, n_effective = NA_integer_))
  }
  tibble(
    bandwidth_m = bw,
    estimate    = rd_fit$coef["Bias-Corrected", ],
    ci_lower    = rd_fit$ci["Robust", "CI Lower"],
    ci_upper    = rd_fit$ci["Robust", "CI Upper"],
    n_effective = sum(rd_fit$N_h)
  )
})

bandwidth_sensitivity_plot <- ggplot(bandwidth_results, aes(x = bandwidth_m / 1000, y = estimate)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 1, fill = "grey80") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_point(size = 2) +
  geom_line() +
  geom_vline(xintercept = c(0.5, 1, 1.5), linetype = "dotted", colour = "black") +
  labs(x = "Bandwidth (km)", y = "Treatment effect (Delta log price)") +
  theme_minimal(base_family = font_family)
bandwidth_sensitivity_plot

# 7. Placebo boundaries (h = 1000m) ----

placebo_shifts <- c(-2000, -1000, -500, 500, 1000, 2000)

placebo_results <- map_dfr(placebo_shifts, function(shift) {
  rd_fit <- tryCatch(
    with(pp_onspd_south_sf_repeatsales,
         rdrobust(y = delta_log_price, x = signed_dist, c = shift,
                  h = 1000, covs = all_covariates,
                  cluster = lsoa21cd, masspoints = "adjust")),
    error = function(e) NULL
  )
  if (is.null(rd_fit)) {
    return(tibble(placebo_shift_m = shift, estimate = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_))
  }
  tibble(
    placebo_shift_m = shift,
    estimate        = rd_fit$coef["Bias-Corrected", ],
    ci_lower        = rd_fit$ci["Robust", "CI Lower"],
    ci_upper        = rd_fit$ci["Robust", "CI Upper"]
  )
})

# Add the actual boundary at 0
placebo_results <- bind_rows(
  tibble(placebo_shift_m = 0,
         estimate = rd_preferred$coef["Bias-Corrected", ],
         ci_lower = rd_preferred$ci["Robust", "CI Lower"],
         ci_upper = rd_preferred$ci["Robust", "CI Upper"]),
  placebo_results
) %>% arrange(placebo_shift_m)

placebo_plot <- ggplot(placebo_results, aes(x = placebo_shift_m / 1000, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
  geom_point(aes(colour = placebo_shift_m == 0), size = 3) +
  scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey30"), guide = "none") +
  labs(x = "Boundary shift (km, negative = moved outward)",
       y = "Treatment effect (Delta log price)") +
  theme_minimal(base_family = font_family)
placebo_plot

# 8. Balance checks ----

# Sale timing and property mix at the boundary (h = 1000m)
balance_vars <- c("year_pre", "year_post", "years_between",
                  "is_flat", "is_terraced", "is_semi", "is_leasehold", "is_new_build")

balance_results <- map_dfr(balance_vars, function(v) {
  rd_bal <- rdrobust(
    y = pp_onspd_south_sf_repeatsales[[v]],
    x = pp_onspd_south_sf_repeatsales$signed_dist,
    c = 0, h = 1000, cluster = pp_onspd_south_sf_repeatsales$lsoa21cd,
    masspoints = "adjust")
  extract_rd(rd_bal, v)
})

# Distance to station, time-invariant (h = 1000m)
balance_station <- extract_rd(
  rdrobust(y = pp_onspd_south_sf_repeatsales$dist_to_station,
           x = pp_onspd_south_sf_repeatsales$signed_dist,
           c = 0, h = 1000, cluster = pp_onspd_south_sf_repeatsales$lsoa21cd,
           masspoints = "adjust"),
  "dist_to_station")

cat("\n=== BALANCE: composition and timing ===\n"); print(balance_results, n = Inf)
cat("\n=== BALANCE: station proximity (cross-sectional diagnostic) ===\n"); print(balance_station, n = Inf)

# Balance table
balance_label <- c(
  year_pre = "Year of pre-sale", year_post = "Year of post-sale",
  years_between = "Years between sales", is_flat = "Flat (share)",
  is_terraced = "Terraced (share)", is_semi = "Semi-detached (share)",
  is_leasehold = "Leasehold (share)", is_new_build = "New build (share)",
  dist_to_station = "Distance to nearest station (m)"
)

format_balance <- function(df) {
  df %>% mutate(
    stars = case_when(robust_pval < 0.01 ~ "***", robust_pval < 0.05 ~ "**",
                      robust_pval < 0.10 ~ "*", TRUE ~ ""),
    label   = balance_label[specification],
    est_str = paste0(round(estimate_bc, 4), stars),
    se_str  = paste0("(", round(robust_se, 4), ")"),
    n_str   = format(n_eff_total, big.mark = ",")
  )
}

bal_comp <- format_balance(balance_results)
bal_station <- format_balance(balance_station)

table_balance_rows <- c(
  list(c("Panel A: Composition and timing", "", "", "")),
  lapply(seq_len(nrow(bal_comp)), function(i)
    c(bal_comp$label[i], bal_comp$est_str[i], bal_comp$se_str[i], bal_comp$n_str[i])),
  list(c("Panel B: Station proximity (diagnostic)", "", "", "")),
  lapply(seq_len(nrow(bal_station)), function(i)
    c(bal_station$label[i], bal_station$est_str[i], bal_station$se_str[i], bal_station$n_str[i]))
)

table_balance_df <- as.data.frame(do.call(rbind, table_balance_rows), stringsAsFactors = FALSE)
names(table_balance_df) <- c(" ", "Discontinuity", "Robust SE", "N")

panel_a_row <- 1
panel_b_row <- 1 + nrow(bal_comp) + 1

ft_balance <- flextable(table_balance_df) %>%
  align(j = 2:4, align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  bold(i = c(panel_a_row, panel_b_row), part = "body") %>%
  italic(i = c(panel_a_row, panel_b_row), part = "body") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) %>%
  hline(i = panel_b_row - 1, part = "body", border = fp_border(color = "grey70", width = 0.5)) %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
  add_footer_lines("Notes: RD discontinuity in each covariate at the boundary, forced h = 1000m, bias-corrected estimate with robust SE in parentheses (CCT 2014). LSOA-clustered. Panel A covariates should be continuous (no differential selection into the repeat-sales sample). Panel B is a cross-sectional RD diagnostic: a discontinuity in time-invariant station proximity is expected and is differenced out under the diff-in-disc design. *** p<0.01, ** p<0.05, * p<0.10.") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit()

# 9. Transport channel ----

# Near vs far from a station (1km)
station_threshold <- 1000
near_station_idx <- pp_onspd_south_sf_repeatsales$dist_to_station <= station_threshold
far_station_idx  <- pp_onspd_south_sf_repeatsales$dist_to_station >  station_threshold

cat("Near station (<=1000m):", sum(near_station_idx), "obs | Far station (>1000m):", sum(far_station_idx), "obs\n")

transport_bandwidths <- c(1000, 1500, 2000)

run_station_rd <- function(group_idx, bw) {
  rdrobust(
    y = pp_onspd_south_sf_repeatsales$delta_log_price[group_idx],
    x = pp_onspd_south_sf_repeatsales$signed_dist[group_idx],
    c = 0, h = bw,
    covs = drop_constant_cols(all_covariates[group_idx, ]),
    cluster = pp_onspd_south_sf_repeatsales$lsoa21cd[group_idx],
    masspoints = "adjust")
}

transport_results <- map_dfr(transport_bandwidths, function(bw) {
  bind_rows(
    extract_rd(run_station_rd(near_station_idx, bw), paste0("Near (h=", bw/1000, "km)")),
    extract_rd(run_station_rd(far_station_idx,  bw), paste0("Far (h=",  bw/1000, "km)"))
  ) %>%
    mutate(station_group = c("Near station (<=1000m)", "Far station (>1000m)"),
           bw_m = bw)
})

cat("\n=== TRANSPORT CHANNEL RESULTS ===\n")
print(transport_results, n = Inf)

# Near/far plot
transport_plot <- ggplot(transport_results %>%
                           mutate(bw_label = paste0(bw_m / 1000, "km"),
                                  station_group = factor(station_group,
                                                         levels = c("Near station (<=1000m)", "Far station (>1000m)"))),
                         aes(x = bw_label, y = estimate_bc)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15) +
  geom_point(size = 3, colour = "steelblue") +
  facet_wrap(~ station_group) +
  labs(x = "Bandwidth", y = "Treatment effect (Delta log price)") +
  theme_minimal(base_family = font_family)
transport_plot

# Near/far table
tr_est <- transport_results %>%
  mutate(
    stars = case_when(robust_pval < 0.01 ~ "***",
                      robust_pval < 0.05 ~ "**",
                      robust_pval < 0.10 ~ "*",
                      TRUE ~ ""),
    est_str = paste0(round(estimate_bc, 4), stars),
    se_str  = paste0("(", round(robust_se, 4), ")"),
    n_str   = format(n_eff_total, big.mark = ","),
    bw_label = paste0("h = ", bw_m / 1000, "km")
  )

tr_near <- tr_est %>% filter(station_group == "Near station (<=1000m)")
tr_far  <- tr_est %>% filter(station_group == "Far station (>1000m)")
col_labels_tr <- tr_near$bw_label

table_tr_rows <- list(
  c("Panel A: Near station (<=1000m)", rep("", length(col_labels_tr))),
  c("ULEZ effect", tr_near$est_str),
  c("",            tr_near$se_str),
  c("N effective", tr_near$n_str),
  c("Panel B: Far station (>1000m)", rep("", length(col_labels_tr))),
  c("ULEZ effect", tr_far$est_str),
  c("",            tr_far$se_str),
  c("N effective", tr_far$n_str)
)

table_tr_df <- as.data.frame(do.call(rbind, table_tr_rows), stringsAsFactors = FALSE)
names(table_tr_df) <- c(" ", col_labels_tr)

ft_transport <- flextable(table_tr_df) %>%
  align(j = 2:ncol(table_tr_df), align = "center", part = "all") %>%
  bold(i = c(1, 5), part = "body") %>%
  italic(i = c(1, 5), part = "body") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) %>%
  hline(i = 4, part = "body", border = fp_border(color = "grey70", width = 0.5)) %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
  add_footer_lines("Notes: Dependent variable is first-differenced log price (post minus pre), estimated separately for properties within and beyond 1000m walking distance of the nearest rail or metro station. A larger effect near stations would indicate a transport (mode-substitution) channel. Property and year-pair covariates included. LSOA-clustered. h = 1000m is preferred for the price outcome; the Green Belt land-use confound widens with bandwidth. Bias-corrected estimates with robust SE in parentheses (CCT 2014). *** p<0.01, ** p<0.05, * p<0.10.") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit()

# Continuous interactions: treatment x moderator (station distance, NO2, PM2.5)
interaction_h <- 1000

pp_onspd_south_sf_repeatsales <- pp_onspd_south_sf_repeatsales %>%
  mutate(
    inside_ulez = as.numeric(signed_dist >= 0),
    run_dist_km = signed_dist / 1000,
    station_km  = dist_to_station / 1000
  )

run_interaction <- function(moderator_col, cluster_col, label) {
  d <- pp_onspd_south_sf_repeatsales
  d$moderator <- d[[moderator_col]]
  d$cluster   <- d[[cluster_col]]
  # Window, centred moderator, triangular kernel weights
  d <- d[abs(d$signed_dist) <= interaction_h & !is.na(d$moderator), ]
  d$moderator_c <- d$moderator - mean(d$moderator)
  d$kernel_w    <- 1 - abs(d$signed_dist) / interaction_h
  
  fit <- feols(delta_log_price ~ inside_ulez * run_dist_km * moderator_c | year_pair,
               data = d, weights = ~kernel_w, cluster = ~cluster)
  
  ct <- coeftable(fit)
  # inside_ulez:moderator_c term
  delta_row <- rownames(ct)[grepl("inside_ulez", rownames(ct)) &
                              grepl("moderator_c", rownames(ct)) &
                              !grepl("run_dist_km", rownames(ct))]
  
  tibble(
    moderator  = label,
    tau        = ct["inside_ulez", "Estimate"],
    tau_se     = ct["inside_ulez", "Std. Error"],
    tau_pval   = ct["inside_ulez", "Pr(>|t|)"],
    delta      = ct[delta_row, "Estimate"],
    delta_se   = ct[delta_row, "Std. Error"],
    delta_pval = ct[delta_row, "Pr(>|t|)"],
    n_obs      = fit$nobs
  )
}

interaction_results <- bind_rows(
  run_interaction("station_km", "lsoa21cd",    "Station distance (per km)"),
  run_interaction("delta_no2",  "pcm_cell_id", "NO2 change (per ug/m3)"),
  run_interaction("delta_pm25", "pcm_cell_id", "PM2.5 change (per ug/m3)")
)

cat("\n=== CONTINUOUS INTERACTION RESULTS ===\n")
print(interaction_results, n = Inf)

# Interaction table
int_est <- interaction_results %>%
  mutate(
    tau_stars    = case_when(tau_pval < 0.01 ~ "***", tau_pval < 0.05 ~ "**",
                             tau_pval < 0.10 ~ "*", TRUE ~ ""),
    delta_stars  = case_when(delta_pval < 0.01 ~ "***", delta_pval < 0.05 ~ "**",
                             delta_pval < 0.10 ~ "*", TRUE ~ ""),
    tau_str      = paste0(round(tau, 4), tau_stars),
    tau_se_str   = paste0("(", round(tau_se, 4), ")"),
    delta_str    = paste0(round(delta, 4), delta_stars),
    delta_se_str = paste0("(", round(delta_se, 4), ")"),
    n_str        = format(n_obs, big.mark = ",")
  )

table_int_df <- data.frame(
  ` ` = c("Jump at mean moderator", "", "Interaction with moderator", "", "N"),
  check.names = FALSE, stringsAsFactors = FALSE
)
for (i in seq_len(nrow(int_est))) {
  table_int_df[[int_est$moderator[i]]] <- c(
    int_est$tau_str[i], int_est$tau_se_str[i],
    int_est$delta_str[i], int_est$delta_se_str[i],
    int_est$n_str[i]
  )
}

ft_interaction <- flextable(table_int_df) %>%
  align(j = 2:ncol(table_int_df), align = "center", part = "all") %>%
  bold(i = c(1, 3), part = "body") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) %>%
  hline(i = 4, part = "body", border = fp_border(color = "grey70", width = 0.5)) %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
  add_footer_lines("Notes: Local-linear WLS, triangular kernel, h = 1000m, year-pair FE, full running-variable interaction. Moderator centred at its within-window mean. Dependent variable is first-differenced log price. Station clustered by LSOA; pollution clustered by PCM 1km cell. Conventional SE in parentheses. *** p<0.01, ** p<0.05, * p<0.10.") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit()

# 10. Air quality channel: RD on NO2 and PM2.5 change, 2022 to 2024 ----

no2_valid_idx  <- !is.na(pp_onspd_south_sf_repeatsales$delta_no2)
pm25_valid_idx <- !is.na(pp_onspd_south_sf_repeatsales$delta_pm25)

channel_bandwidths <- c(1000, 1500, 2000)

# PCM cells within the bandwidth
count_cells <- function(valid_idx, bw) {
  d <- pp_onspd_south_sf_repeatsales[valid_idx, ]
  length(unique(d$pcm_cell_id[abs(d$signed_dist) <= bw]))
}

channel_results <- map_dfr(channel_bandwidths, function(bw) {
  rd_no2 <- rdrobust(
    y = pp_onspd_south_sf_repeatsales$delta_no2[no2_valid_idx],
    x = pp_onspd_south_sf_repeatsales$signed_dist[no2_valid_idx],
    c = 0, h = bw,
    cluster = pp_onspd_south_sf_repeatsales$pcm_cell_id[no2_valid_idx],
    masspoints = "adjust")
  
  rd_pm25 <- rdrobust(
    y = pp_onspd_south_sf_repeatsales$delta_pm25[pm25_valid_idx],
    x = pp_onspd_south_sf_repeatsales$signed_dist[pm25_valid_idx],
    c = 0, h = bw,
    cluster = pp_onspd_south_sf_repeatsales$pcm_cell_id[pm25_valid_idx],
    masspoints = "adjust")
  
  bind_rows(
    extract_rd(rd_no2,  paste0("NO2 (h=", bw/1000, "km)")),
    extract_rd(rd_pm25, paste0("PM2.5 (h=", bw/1000, "km)"))
  ) %>%
    mutate(pollutant = c("NO2", "PM2.5"),
           bw_m = bw,
           n_cells = c(count_cells(no2_valid_idx, bw), count_cells(pm25_valid_idx, bw)))
})

cat("\n=== TRANSMISSION CHANNEL RESULTS ===\n")
print(channel_results, n = Inf)

# Channel plot
channel_plot <- ggplot(channel_results %>% mutate(bw_label = paste0(bw_m/1000, "km")),
                       aes(x = bw_label, y = estimate_bc)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15) +
  geom_point(size = 3, colour = "darkred") +
  facet_wrap(~ pollutant, scales = "free_y") +
  labs(x = "Bandwidth", y = "Discontinuity in pollution change (ug/m3)") +
  theme_minimal(base_family = font_family)
channel_plot

# Channel table
ch_est <- channel_results %>%
  mutate(
    stars = case_when(robust_pval < 0.01 ~ "***",
                      robust_pval < 0.05 ~ "**",
                      robust_pval < 0.10 ~ "*",
                      TRUE ~ ""),
    est_str = paste0(round(estimate_bc, 4), stars),
    se_str  = paste0("(", round(robust_se, 4), ")"),
    n_str   = format(n_eff_total, big.mark = ","),
    cells_str = format(n_cells, big.mark = ","),
    bw_label = paste0("h = ", bw_m/1000, "km")
  )

ch_no2  <- ch_est %>% filter(pollutant == "NO2")
ch_pm25 <- ch_est %>% filter(pollutant == "PM2.5")
col_labels_ch <- ch_no2$bw_label

table_ch_rows <- list(
  c("Panel A: NO2 (ug/m3)", rep("", length(col_labels_ch))),
  c("ULEZ effect", ch_no2$est_str),
  c("",            ch_no2$se_str),
  c("N effective", ch_no2$n_str),
  c("N cells",     ch_no2$cells_str),
  c("Panel B: PM2.5 (ug/m3)", rep("", length(col_labels_ch))),
  c("ULEZ effect", ch_pm25$est_str),
  c("",            ch_pm25$se_str),
  c("N effective", ch_pm25$n_str),
  c("N cells",     ch_pm25$cells_str)
)

table_ch_df <- as.data.frame(do.call(rbind, table_ch_rows), stringsAsFactors = FALSE)
names(table_ch_df) <- c(" ", col_labels_ch)

ft_channel <- flextable(table_ch_df) %>%
  align(j = 2:ncol(table_ch_df), align = "center", part = "all") %>%
  bold(i = c(1, 6), part = "body") %>%
  italic(i = c(1, 6), part = "body") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) %>%
  hline(i = 5, part = "body", border = fp_border(color = "grey70", width = 0.5)) %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
  add_footer_lines("Notes: Dependent variable is change in Defra PCM modelled annual mean concentration (2022 to 2024). Negative = improvement inside ULEZ. Inference clustered at the PCM 1km grid cell, the level at which the pollution outcome varies; N cells is the number of distinct cells within the bandwidth. Wider bandwidths defensible: the Green Belt land-use confound affects prices but not atmospheric concentrations. Bias-corrected estimates with robust SE (CCT 2014). *** p<0.01, ** p<0.05, * p<0.10.") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit()

# 11. Export tables to Word ----

doc_all <- read_docx() %>%
  body_add_par("Regression results", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft_table1) %>%
  body_add_break() %>%
  body_add_flextable(ft_channel) %>%
  body_add_break() %>%
  body_add_flextable(ft_transport) %>%
  body_add_break() %>%
  body_add_flextable(ft_interaction) %>%
  body_add_break() %>%
  body_add_par("Appendix", style = "heading 1") %>%
  body_add_flextable(ft_table_appendix) %>%
  body_add_break() %>%
  body_add_flextable(ft_balance)

print(doc_all, target = here("output", "regression_tables.docx"))
cat("\nSaved: output/regression_tables.docx\n")

# 12. Save plots ----

ggsave(here("output", "rd_main_plot.png"), rd_main_plot, width = 8, height = 5)
ggsave(here("output", "bandwidth_sensitivity_plot.png"), bandwidth_sensitivity_plot, width = 8, height = 5)
ggsave(here("output", "placebo_plot.png"), placebo_plot, width = 8, height = 5)
ggsave(here("output", "transport_plot.png"), transport_plot, width = 8, height = 5)
ggsave(here("output", "channel_plot.png"), channel_plot, width = 8, height = 5)