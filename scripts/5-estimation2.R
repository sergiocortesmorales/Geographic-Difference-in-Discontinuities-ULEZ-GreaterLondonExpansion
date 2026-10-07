# Install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, tidyverse, sf, rdrobust, flextable, officer)

# Read files
pp_onspd_south_sf <- readRDS(here("data", "pp_onspd_south_sf.rds"))

dir.create(here("output"), showWarnings = FALSE)

# Plot font
font_family <- "serif"

# 1. Repeated cross-section (all sales) ----

pp_onspd_south_allsales <- pp_onspd_south_sf %>%
  st_drop_geometry() %>%
  filter(post_ulez %in% c("Pre", "Post")) %>%
  mutate(
    property_id  = paste(pcds, paon, saon, street, sep = "|"),
    sale_year    = as.numeric(format(date_full, "%Y")),
    is_flat      = as.numeric(property_type == "F"),
    is_terraced  = as.numeric(property_type == "T"),
    is_semi      = as.numeric(property_type == "S"),
    is_leasehold = as.numeric(tenure == "L"),
    is_new_build = as.numeric(new_build == "Y")
  )

rm(pp_onspd_south_sf)

# Period indicators
post_idx <- pp_onspd_south_allsales$post_ulez == "Post"
pre_idx  <- pp_onspd_south_allsales$post_ulez == "Pre"

cat("All sales:", nrow(pp_onspd_south_allsales),
    "| Post:", sum(post_idx), "| Pre:", sum(pre_idx), "\n")
cat("Within 1km  -> Post:", sum(post_idx & abs(pp_onspd_south_allsales$signed_dist) <= 1000),
    "| Pre:", sum(pre_idx & abs(pp_onspd_south_allsales$signed_dist) <= 1000), "\n")

# 2. Helpers ----

# Drop covariates that are constant within a subsample
drop_constant_cols <- function(covs_matrix) {
  varying <- apply(covs_matrix, 2, function(col) length(unique(col)) > 1)
  covs_matrix[, varying, drop = FALSE]
}

# Covariate matrix for a subsample: none, property, or property + sale-year FE
build_covs <- function(mask, mode) {
  if (mode == "none") return(NULL)
  d <- pp_onspd_south_allsales[mask, ]
  prop <- cbind(is_flat = d$is_flat, is_terraced = d$is_terraced, is_semi = d$is_semi,
                is_leasehold = d$is_leasehold, is_new_build = d$is_new_build)
  m <- if (mode == "property") prop else cbind(prop, model.matrix(~ factor(d$sale_year))[, -1, drop = FALSE])
  m <- drop_constant_cols(m)
  if (ncol(m) == 0) NULL else m
}

# rdrobust on log price for one subsample
fit_side <- function(mask, h_val, poly_order, mode, cutoff = 0) {
  d <- pp_onspd_south_allsales[mask, ]
  covs <- build_covs(mask, mode)
  args <- list(y = d$log_price, x = d$signed_dist, c = cutoff, p = poly_order,
               cluster = d$lsoa21cd, masspoints = "adjust")
  if (!is.null(covs)) args$covs <- covs
  if (!is.null(h_val)) args$h <- h_val
  do.call(rdrobust, args)
}

# Post and pre RDs on the same window
run_didc <- function(base_mask, h_val, poly_order, mode, cutoff = 0) {
  post_mask <- base_mask & post_idx
  pre_mask  <- base_mask & pre_idx
  rd_post <- fit_side(post_mask, h_val, poly_order, mode, cutoff)
  rd_pre  <- fit_side(pre_mask,  h_val, poly_order, mode, cutoff)
  list(post = rd_post, pre = rd_pre)
}

# Post minus pre estimate, SE assumes the two RDs are independent
extract_didc <- function(rd_post, rd_pre, label = "") {
  est <- rd_post$coef["Bias-Corrected", ] - rd_pre$coef["Bias-Corrected", ]
  se  <- sqrt(rd_post$se["Robust", ]^2 + rd_pre$se["Robust", ]^2)
  z   <- est / se
  tibble(
    specification = label,
    estimate_bc   = est,
    robust_se     = se,
    robust_z      = z,
    robust_pval   = 2 * pnorm(-abs(z)),
    ci_lower      = est - qnorm(0.975) * se,
    ci_upper      = est + qnorm(0.975) * se,
    bandwidth     = rd_post$bws["h", "left"],
    n_eff_post    = sum(rd_post$N_h),
    n_eff_pre     = sum(rd_pre$N_h),
    n_eff_total   = sum(rd_post$N_h) + sum(rd_pre$N_h)
  )
}

# Single RD summary (station balance and pollution outcomes)
extract_rd <- function(rd_obj, label = "") {
  tibble(
    specification = label,
    estimate_bc   = rd_obj$coef["Bias-Corrected", ],
    robust_se     = rd_obj$se["Robust", ],
    robust_z      = rd_obj$z["Robust", ],
    robust_pval   = rd_obj$pv["Robust", ],
    ci_lower      = rd_obj$ci["Robust", "CI Lower"],
    ci_upper      = rd_obj$ci["Robust", "CI Upper"],
    bandwidth     = rd_obj$bws["h", "left"],
    n_eff_total   = sum(rd_obj$N_h)
  )
}

# 3. Main estimates: 3 bandwidths x 2 polynomial orders x 5 specifications ----

bandwidth_configs <- list(
  list(h = 500,  label = "h = 500m"),
  list(h = 1000, label = "h = 1000m"),
  list(h = 1500, label = "h = 1500m")
)

# Full sample and donut masks
all_mask        <- rep(TRUE, nrow(pp_onspd_south_allsales))
donut_100_mask  <- abs(pp_onspd_south_allsales$signed_dist) > 100
donut_200_mask  <- abs(pp_onspd_south_allsales$signed_dist) > 200

run_five_specs <- function(poly_order, h_val) {
  s_base <- run_didc(all_mask,        h_val, poly_order, "none")
  s_prop <- run_didc(all_mask,        h_val, poly_order, "property")
  s_all  <- run_didc(all_mask,        h_val, poly_order, "all")
  s_d100 <- run_didc(donut_100_mask,  h_val, poly_order, "all")
  s_d200 <- run_didc(donut_200_mask,  h_val, poly_order, "all")
  bind_rows(
    extract_didc(s_base$post, s_base$pre, "Baseline"),
    extract_didc(s_prop$post, s_prop$pre, "Property covs"),
    extract_didc(s_all$post,  s_all$pre,  "Prop + Year FE"),
    extract_didc(s_d100$post, s_d100$pre, "Donut 100m"),
    extract_didc(s_d200$post, s_d200$pre, "Donut 200m")
  )
}

all_results <- map_dfr(bandwidth_configs, function(cfg) {
  map_dfr(c(1, 2), function(p) {
    cat("Running:", cfg$label, "| p =", p, "\n")
    run_five_specs(poly_order = p, h_val = cfg$h) %>%
      mutate(bw_label = cfg$label, poly_order = p)
  })
})

cat("\n=== ALL RESULTS ===\n")
print(all_results, n = Inf)

# Preferred: h = 1000m, linear, property + year FE
didc_preferred_obj <- run_didc(all_mask, 1000, 1, "all")
didc_preferred <- extract_didc(didc_preferred_obj$post, didc_preferred_obj$pre, "Full sample")

# 4. Result tables: linear (main) and quadratic (appendix) ----

build_panel_table <- function(results_df, poly, table_notes) {
  df <- results_df %>% filter(poly_order == poly)
  
  col_labels <- c("Baseline", "Property covs", "Prop + Year FE", "Donut 100m", "Donut 200m")
  bw_labels  <- c("h = 500m", "h = 1000m", "h = 1500m")
  
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
        bw_str  = format(round(bandwidth, 0), big.mark = ","),
        n_str   = format(n_eff_total, big.mark = ",")
      )
    
    has_prop     <- c("", "\u2713", "\u2713", "\u2713", "\u2713")
    has_yearfe   <- c("", "", "\u2713", "\u2713", "\u2713")
    
    rows <- c(rows, list(
      c(bw_lab, rep("", 5)),
      c("ULEZ effect",         est_stars$est_str),
      c("",                    est_stars$se_str),
      c("Property covariates", has_prop),
      c("Year FE",             has_yearfe),
      c("Bandwidth (m)",       est_stars$bw_str),
      c("N effective",         est_stars$n_str)
    ))
  }
  
  table_df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(table_df) <- c(" ", col_labels)
  
  panel_rows <- which(table_df[, 1] %in% bw_labels)
  
  ft <- flextable(table_df) %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    align(j = 2:6, align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "all") %>%
    bold(i = panel_rows, part = "body") %>%
    italic(i = panel_rows, part = "body") %>%
    border_remove() %>%
    hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) %>%
    hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
    autofit() %>%
    add_footer_lines(table_notes)
  
  for (pr in panel_rows) {
    ft <- hline(ft, i = pr, part = "body",
                border = fp_border(color = "grey70", width = 0.5))
  }
  
  ft
}

main_notes <- "Notes: Dependent variable is log price. Estimate is the post-period boundary discontinuity minus the pre-period boundary discontinuity, each a bias-corrected RD on log price levels (CCT 2014), clustered at LSOA. Combined robust SE in parentheses treats the two period RDs as independent (conservative). Both period RDs use the same window. Year FE are sale-year dummies within each period. Donut columns exclude properties within the stated radius. N effective is the pre plus post effective sample. *** p<0.01, ** p<0.05, * p<0.10."

ft_table1 <- build_panel_table(all_results, poly = 1, table_notes = main_notes)
ft_table_appendix <- build_panel_table(all_results, poly = 2, table_notes = main_notes)

# 5. RD plot by period ----

plot_sub <- filter(pp_onspd_south_allsales, abs(signed_dist) <= 2000)

make_period_bins <- function(period_label) {
  d <- plot_sub[plot_sub$post_ulez == period_label, ]
  rp <- rdplot(y = d$log_price, x = d$signed_dist, c = 0, p = 1,
               nbins = c(20, 20), masspoints = "adjust", hide = TRUE)
  # Split the fitted line at the cutoff
  cut_row <- which(rp$vars_poly$rdplot_x >= 0)[1]
  list(
    bins = tibble(period = period_label,
                  x = rp$vars_bins$rdplot_mean_x, y = rp$vars_bins$rdplot_mean_y),
    poly = tibble(period = period_label,
                  x = rp$vars_poly$rdplot_x, y = rp$vars_poly$rdplot_y,
                  side = c(rep("Outside ULEZ", cut_row),
                           rep("Inside ULEZ", nrow(rp$vars_poly) - cut_row)))
  )
}

pre_bins  <- make_period_bins("Pre")
post_bins <- make_period_bins("Post")

period_levels <- c("Pre", "Post")
period_labels <- c("Pre-ULEZ", "Post-ULEZ")

rd_main_bins <- bind_rows(pre_bins$bins, post_bins$bins) %>%
  mutate(period = factor(period, levels = period_levels, labels = period_labels))
rd_main_poly <- bind_rows(pre_bins$poly, post_bins$poly) %>%
  mutate(period = factor(period, levels = period_levels, labels = period_labels))

rd_main_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  geom_point(data = rd_main_bins, aes(x / 1000, y), alpha = 0.6, colour = "grey30", size = 1.2) +
  geom_line(data = rd_main_poly, aes(x / 1000, y, colour = side), linewidth = 0.8) +
  scale_colour_manual(values = c("Outside ULEZ" = "darkorange", "Inside ULEZ" = "darkgreen")) +
  facet_wrap(~ period, scales = "free_y") +
  labs(x = "Distance to boundary (km, positive = inside ULEZ)",
       y = "Log price", colour = NULL) +
  theme_minimal(base_family = font_family) +
  theme(legend.position = "bottom")
rd_main_plot

# 6. Bandwidth sensitivity ----

bandwidths <- seq(500, 5000, by = 250)

bandwidth_results <- map_dfr(bandwidths, function(bw) {
  obj <- tryCatch(run_didc(all_mask, bw, 1, "all"), error = function(e) NULL)
  if (is.null(obj)) {
    return(tibble(bandwidth_m = bw, estimate = NA_real_, ci_lower = NA_real_,
                  ci_upper = NA_real_, n_effective = NA_integer_))
  }
  r <- extract_didc(obj$post, obj$pre)
  tibble(bandwidth_m = bw, estimate = r$estimate_bc, ci_lower = r$ci_lower,
         ci_upper = r$ci_upper, n_effective = r$n_eff_total)
})

bandwidth_sensitivity_plot <- ggplot(bandwidth_results, aes(x = bandwidth_m / 1000, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
  geom_point(size = 2) +
  geom_line() +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dotted", colour = "darkgreen", alpha = 0.5) +
  annotate("text", x = 1.05, y = max(bandwidth_results$ci_upper, na.rm = TRUE),
           label = "Preferred (1km)", hjust = 0, size = 3, colour = "red", family = font_family) +
  annotate("text", x = 2.05, y = max(bandwidth_results$ci_upper, na.rm = TRUE),
           label = "Green Belt threshold", hjust = 0, size = 3, colour = "darkgreen", family = font_family) +
  labs(x = "Bandwidth (km)", y = "Treatment effect (post - pre discontinuity)") +
  theme_minimal(base_family = font_family)
bandwidth_sensitivity_plot

# 7. Placebo boundaries (h = 1000m) ----

placebo_shifts <- c(-2000, -1000, -500, 500, 1000, 2000)

placebo_results <- map_dfr(placebo_shifts, function(shift) {
  obj <- tryCatch(run_didc(all_mask, 1000, 1, "all", cutoff = shift), error = function(e) NULL)
  if (is.null(obj)) {
    return(tibble(placebo_shift_m = shift, estimate = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_))
  }
  r <- extract_didc(obj$post, obj$pre)
  tibble(placebo_shift_m = shift, estimate = r$estimate_bc,
         ci_lower = r$ci_lower, ci_upper = r$ci_upper)
})

# Add the actual boundary at 0
placebo_results <- bind_rows(
  tibble(placebo_shift_m = 0,
         estimate = didc_preferred$estimate_bc,
         ci_lower = didc_preferred$ci_lower,
         ci_upper = didc_preferred$ci_upper),
  placebo_results
) %>% arrange(placebo_shift_m)

placebo_plot <- ggplot(placebo_results, aes(x = placebo_shift_m / 1000, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
  geom_point(aes(colour = placebo_shift_m == 0), size = 3) +
  scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey30"), guide = "none") +
  labs(x = "Boundary shift (km, negative = moved outward)",
       y = "Treatment effect (post - pre discontinuity)") +
  theme_minimal(base_family = font_family)
placebo_plot

# 8. Balance checks ----

# Post minus pre discontinuity in property mix (h = 1000m)
balance_vars <- c("is_flat", "is_terraced", "is_semi", "is_leasehold", "is_new_build")

post_df <- pp_onspd_south_allsales[post_idx, ]
pre_df  <- pp_onspd_south_allsales[pre_idx, ]

balance_results <- map_dfr(balance_vars, function(v) {
  rd_post <- rdrobust(y = post_df[[v]], x = post_df$signed_dist, c = 0, h = 1000,
                      cluster = post_df$lsoa21cd, masspoints = "adjust")
  rd_pre  <- rdrobust(y = pre_df[[v]],  x = pre_df$signed_dist,  c = 0, h = 1000,
                      cluster = pre_df$lsoa21cd,  masspoints = "adjust")
  extract_didc(rd_post, rd_pre, v)
})

# Distance to station, one row per property, time-invariant
station_sample <- distinct(pp_onspd_south_allsales, property_id, .keep_all = TRUE)

balance_station <- extract_rd(
  rdrobust(y = station_sample$dist_to_station, x = station_sample$signed_dist,
           c = 0, cluster = station_sample$lsoa21cd, masspoints = "adjust"),
  "dist_to_station")

cat("\n=== BALANCE: composition stability (post - pre discontinuity) ===\n"); print(balance_results, n = Inf)
cat("\n=== BALANCE: station proximity (cross-sectional diagnostic) ===\n"); print(balance_station, n = Inf)

# Balance table
balance_label <- c(
  is_flat = "Flat (share)", is_terraced = "Terraced (share)",
  is_semi = "Semi-detached (share)", is_leasehold = "Leasehold (share)",
  is_new_build = "New build (share)",
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

bal_comp    <- format_balance(balance_results)
bal_station <- format_balance(balance_station)

table_balance_rows <- c(
  list(c("Panel A: Composition stability (post - pre)", "", "", "")),
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
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  align(j = 2:4, align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  bold(i = c(panel_a_row, panel_b_row), part = "body") %>%
  italic(i = c(panel_a_row, panel_b_row), part = "body") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) %>%
  hline(i = panel_b_row - 1, part = "body", border = fp_border(color = "grey70", width = 0.5)) %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
  autofit() %>%
  add_footer_lines("Notes: Panel A is the post-period minus pre-period RD discontinuity in each covariate (forced h = 1000m), bias-corrected with combined robust SE in parentheses (CCT 2014), LSOA-clustered. Values near zero mean the composition discontinuity is time-invariant and differences out under the design. Panel B is a single cross-sectional RD on time-invariant station proximity (MSE-optimal bandwidth, one row per property); a discontinuity is expected and differences out. *** p<0.01, ** p<0.05, * p<0.10.")

# 9. Heterogeneity: station distance and property type (h = 1000m) ----

# 9a. Near vs far from a station
station_threshold <- 800   # ~10 min walk
near_mask <- pp_onspd_south_allsales$dist_to_station <= station_threshold
far_mask  <- pp_onspd_south_allsales$dist_to_station >  station_threshold

cat("Near station (<=800m):", sum(near_mask), "sales | Far station (>800m):", sum(far_mask), "sales\n")

obj_near <- run_didc(near_mask, 1000, 1, "all")
obj_far  <- run_didc(far_mask,  1000, 1, "all")

# 9b. Flats vs non-flats
flat_mask    <- pp_onspd_south_allsales$is_flat == 1
nonflat_mask <- pp_onspd_south_allsales$is_flat == 0

cat("Flats:", sum(flat_mask), "sales | Non-flats:", sum(nonflat_mask), "sales\n")

obj_flat    <- run_didc(flat_mask,    1000, 1, "all")
obj_nonflat <- run_didc(nonflat_mask, 1000, 1, "all")

heterogeneity_results <- bind_rows(
  didc_preferred,
  extract_didc(obj_near$post,    obj_near$pre,    "Near station (<=800m)"),
  extract_didc(obj_far$post,     obj_far$pre,     "Far station (>800m)"),
  extract_didc(obj_flat$post,    obj_flat$pre,    "Flats"),
  extract_didc(obj_nonflat$post, obj_nonflat$pre, "Non-flats (D/S/T)")
)

cat("\n=== HETEROGENEITY RESULTS ===\n")
print(heterogeneity_results, n = Inf)

# Heterogeneity plot
heterogeneity_results_plot <- heterogeneity_results %>%
  mutate(
    dimension = c("Baseline", "Station", "Station", "Property", "Property"),
    specification = factor(specification,
                           levels = c("Full sample", "Near station (<=800m)", "Far station (>800m)", "Flats", "Non-flats (D/S/T)"))
  )

heterogeneity_plot <- ggplot(heterogeneity_results_plot,
                             aes(x = specification, y = estimate_bc)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15) +
  geom_point(aes(colour = dimension), size = 3) +
  scale_colour_manual(values = c("Baseline" = "red", "Station" = "steelblue",
                                 "Property" = "darkorange"), name = "Split") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  labs(x = NULL, y = "Treatment effect (post - pre discontinuity)") +
  theme_minimal(base_family = font_family)
heterogeneity_plot

# Heterogeneity table
het_est <- heterogeneity_results %>%
  mutate(
    stars = case_when(robust_pval < 0.01 ~ "***",
                      robust_pval < 0.05 ~ "**",
                      robust_pval < 0.10 ~ "*",
                      TRUE ~ ""),
    est_str = paste0(round(estimate_bc, 4), stars),
    se_str  = paste0("(", round(robust_se, 4), ")"),
    n_str   = format(n_eff_total, big.mark = ",")
  )

table_het_df <- data.frame(
  ` ` = c("ULEZ effect", "", "N effective"),
  check.names = FALSE, stringsAsFactors = FALSE
)
for (i in seq_len(nrow(het_est))) {
  table_het_df[[het_est$specification[i]]] <- c(
    het_est$est_str[i], het_est$se_str[i], het_est$n_str[i]
  )
}

ft_heterogeneity <- flextable(table_het_df) %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  align(j = 2:ncol(table_het_df), align = "center", part = "all") %>%
  bold(i = 1, part = "body") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) %>%
  hline(i = 2, part = "body", border = fp_border(color = "grey70", width = 0.5)) %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
  autofit() %>%
  add_footer_lines("Notes: Post-period minus pre-period boundary discontinuity in log price, bias-corrected with combined robust SE in parentheses (CCT 2014), LSOA-clustered. All splits use forced h = 1000m. Station split at the 800m walking catchment to the nearest rail or metro station. Property split: flats and non-flats use year FE only. N effective is the pre plus post effective sample. *** p<0.01, ** p<0.05, * p<0.10.")

# 10. Air quality channel: RD on NO2 and PM2.5 change, 2022 to 2024 ----

# One row per property, clustered by PCM cell
channel_sample <- distinct(pp_onspd_south_allsales, property_id, .keep_all = TRUE)

no2_valid_idx  <- !is.na(channel_sample$delta_no2)
pm25_valid_idx <- !is.na(channel_sample$delta_pm25)

channel_bandwidths <- c(1000, 1500, 2000)

# PCM cells within the bandwidth
count_cells <- function(valid_idx, bw) {
  d <- channel_sample[valid_idx, ]
  length(unique(d$pcm_cell_id[abs(d$signed_dist) <= bw]))
}

channel_results <- map_dfr(channel_bandwidths, function(bw) {
  rd_no2 <- rdrobust(
    y = channel_sample$delta_no2[no2_valid_idx],
    x = channel_sample$signed_dist[no2_valid_idx],
    c = 0, h = bw,
    cluster = channel_sample$pcm_cell_id[no2_valid_idx],
    masspoints = "adjust")
  
  rd_pm25 <- rdrobust(
    y = channel_sample$delta_pm25[pm25_valid_idx],
    x = channel_sample$signed_dist[pm25_valid_idx],
    c = 0, h = bw,
    cluster = channel_sample$pcm_cell_id[pm25_valid_idx],
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
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  align(j = 2:ncol(table_ch_df), align = "center", part = "all") %>%
  bold(i = c(1, 6), part = "body") %>%
  italic(i = c(1, 6), part = "body") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) %>%
  hline(i = 5, part = "body", border = fp_border(color = "grey70", width = 0.5)) %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
  autofit() %>%
  add_footer_lines("Notes: Dependent variable is the change in Defra PCM modelled annual mean concentration (2022 to 2024), one row per property. Negative = improvement inside ULEZ. Inference clustered at the PCM 1km grid cell, the level the pollution outcome varies on; N cells is the number of distinct cells within the bandwidth. Wider bandwidths defensible: the Green Belt land-use confound affects prices but not atmospheric concentrations. Bias-corrected estimates with robust SE (CCT 2014). *** p<0.01, ** p<0.05, * p<0.10.")

# 11. Export tables to Word ----

doc_all <- read_docx() %>%
  body_add_par("Regression results (repeated cross-section)", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft_table1) %>%
  body_add_break() %>%
  body_add_flextable(ft_heterogeneity) %>%
  body_add_break() %>%
  body_add_flextable(ft_channel) %>%
  body_add_break() %>%
  body_add_par("Appendix", style = "heading 1") %>%
  body_add_flextable(ft_table_appendix) %>%
  body_add_break() %>%
  body_add_flextable(ft_balance)

print(doc_all, target = here("output", "regression_tables_crosssection.docx"))
cat("\nSaved: output/regression_tables_crosssection.docx\n")

# 12. Save plots ----

ggsave(here("output", "rd_main_plot_cs.png"), rd_main_plot, width = 9, height = 5)
ggsave(here("output", "bandwidth_sensitivity_plot_cs.png"), bandwidth_sensitivity_plot, width = 8, height = 5)
ggsave(here("output", "placebo_plot_cs.png"), placebo_plot, width = 8, height = 5)
ggsave(here("output", "heterogeneity_plot_cs.png"), heterogeneity_plot, width = 8, height = 5)
ggsave(here("output", "channel_plot_cs.png"), channel_plot, width = 8, height = 5)
cat("All plots saved to output/\n")