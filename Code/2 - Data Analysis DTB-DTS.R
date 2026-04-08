##########################################################################
#####                    Data Analysis – DB Check                    #####
##########################################################################

# ------------------------------------------------------------------------
# Overview:
# This script analyses the direct-to-buffer and direct-to-swab control
# experiments run to validate stock DNA concentrations (cell-free trout
# and cellular mouse) and to establish a maximum-recoverable baseline
# for transfer efficiency calculations.
#
# The workflow follows the same structure as scripts 1–3:
#   1. Filter unknowns, remove blanks and undetermined wells
#   2. Calculate Quantity_Total (ng in 250 µL eluate = Quantity × 50)
#   3. Attach sample metadata via lookup (concentration, route, DNA type)
#   4. Average PCR duplicate pairs within each replicate
#   5. Pivot to Cell_DNA / Cell_free_DNA columns and sum to Total_DNA
#   6. QC: standard curve plots (slope, R², efficiency)
#   7. Recovery plots (ng recovered and % recovery vs nominal input)
#   8. Export summary CSVs for supplementary information
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Section 1: Sample metadata lookup tables
# ------------------------------------------------------------------------

# Naming convention used on the plate:
#   D       = Direct deposit to hand
#   B / S   = Buffer / Swab extraction
#   F / C / M = Free (cell-free trout) / Cellular (mouse) / Mixed (both species)
#   1 / 2 / 3 = Concentration level (high / mid / low)
#
# Note: SDF3 / SDC3 appear in the plate layout in place of DSF3 / DSC3
#       (plate naming inconsistency; treated as equivalent below).
#
# Nominal input per 20 µL aliquot (ng):
#   Single species   1 = 35,   2 = 25,   3 = 5
#   Mixed per sp.    1 = 17.5, 2 = 12.5, 3 = 2.5

nominal_ng <- c(
  DBF1 = 35,   DBF2 = 25,   DBF3 = 5,
  DSF1 = 35,   DSF2 = 25,   DSF3 = 5,   SDF3 = 5,
  DBC1 = 35,   DBC2 = 25,   DBC3 = 5,
  DSC1 = 35,   DSC2 = 25,   DSC3 = 5,   SDC3 = 5,
  DBM1 = 17.5, DBM2 = 12.5, DBM3 = 2.5,
  DSM1 = 17.5, DSM2 = 12.5, DSM3 = 2.5
)

extraction_route <- c(
  DBF1 = "Buffer", DBF2 = "Buffer", DBF3 = "Buffer",
  DSF1 = "Swab",   DSF2 = "Swab",   DSF3 = "Swab",   SDF3 = "Swab",
  DBC1 = "Buffer", DBC2 = "Buffer", DBC3 = "Buffer",
  DSC1 = "Swab",   DSC2 = "Swab",   DSC3 = "Swab",   SDC3 = "Swab",
  DBM1 = "Buffer", DBM2 = "Buffer", DBM3 = "Buffer",
  DSM1 = "Swab",   DSM2 = "Swab",   DSM3 = "Swab"
)

dna_type <- c(
  DBF1 = "Cell-free (trout)", DBF2 = "Cell-free (trout)", DBF3 = "Cell-free (trout)",
  DSF1 = "Cell-free (trout)", DSF2 = "Cell-free (trout)", DSF3 = "Cell-free (trout)",
  SDF3 = "Cell-free (trout)",
  DBC1 = "Cellular (mouse)",  DBC2 = "Cellular (mouse)",  DBC3 = "Cellular (mouse)",
  DSC1 = "Cellular (mouse)",  DSC2 = "Cellular (mouse)",  DSC3 = "Cellular (mouse)",
  SDC3 = "Cellular (mouse)",
  DBM1 = "Mixed", DBM2 = "Mixed", DBM3 = "Mixed",
  DSM1 = "Mixed", DSM2 = "Mixed", DSM3 = "Mixed"
)

# Attach metadata and calculate Quantity_Total (same × 50 factor as most studies)
attach_metadata <- function(df) {
  df %>%
    mutate(
      Quantity       = as.numeric(Quantity),
      Quantity_Total = Quantity * 50,
      Nominal_ng     = nominal_ng[`Sample Name`],
      Route          = extraction_route[`Sample Name`],
      DNA_type       = dna_type[`Sample Name`]
    ) %>%
    filter(!is.na(Nominal_ng))   # drops any unrecognised sample codes
}

trout_data <- attach_metadata(trout_unk)
mouse_data <- attach_metadata(mouse_unk)

# ------------------------------------------------------------------------
# Section 2: Outlier removal
# ------------------------------------------------------------------------

# Two SDC3 replicates (wells A9 and B9 on the mouse plate) are removed.
# Their Quantity values (0.3385 and 0.3848 ng / 5 uL) are 3-4x higher
# than the other four replicates (0.095-0.112 ng / 5 uL). The anomalously
# low Ct values (~24.0 vs ~26.0) are consistent with a pipetting artefact.
# This matches the exclusion made manually in the processed data spreadsheet.

mouse_data_clean <- mouse_data %>%
  filter(!(Well %in% c("A9", "B9")))

trout_data_clean <- trout_data   # no exclusions in the trout plate

# ------------------------------------------------------------------------
# Section 3: Average PCR duplicates, then rename targets
# ------------------------------------------------------------------------

# Each sample has 6 wells = 3 experimental replicates run in PCR duplicate.
# Consecutive pairs (wells 1+2, 3+4, 5+6) belong to the same replicate.
# For each pair: take the mean if both values are non-zero, otherwise keep
# the non-zero value.

average_replicates <- function(df) {
  df %>%
    arrange(`Sample Name`, `Target Name`, Well) %>%
    group_by(`Sample Name`, `Target Name`, Route, DNA_type, Nominal_ng) %>%
    mutate(PCR_pair = ceiling(row_number() / 2)) %>%
    group_by(`Sample Name`, `Target Name`, Route, DNA_type, Nominal_ng, PCR_pair) %>%
    summarise(
      Replicate_Quantity = if (all(Quantity_Total == 0)) {
        0
      } else {
        mean(Quantity_Total[Quantity_Total != 0], na.rm = TRUE)
      },
      .groups = "drop"
    )
}

trout_avg <- average_replicates(trout_data_clean)
mouse_avg <- average_replicates(mouse_data_clean)

# Rename Target Name to match the convention used in scripts 2-3
trout_avg <- trout_avg %>%
  mutate(`Target Name` = recode(`Target Name`, "TROUT 1" = "Cell_free_DNA"))
mouse_avg <- mouse_avg %>%
  mutate(`Target Name` = recode(`Target Name`, "Mouse" = "Cell_DNA"))


# ------------------------------------------------------------------------
# Section 4: Pivot to Cell_DNA / Cell_free_DNA columns (same as scripts 2–3)
# ------------------------------------------------------------------------
pivot_db <- function(df) {
  df %>%
    pivot_wider(names_from = `Target Name`, values_from = Replicate_Quantity) %>%
    # Each single-species file will be missing one column after pivot_wider;
    # create it as zero before mutating so the function works for both files
    # and for the mixed samples where both columns are present.
    { if (!"Cell_DNA"      %in% names(.)) mutate(., Cell_DNA      = 0) else . } %>%
    { if (!"Cell_free_DNA" %in% names(.)) mutate(., Cell_free_DNA = 0) else . } %>%
    mutate(
      Cell_DNA      = replace_na(Cell_DNA,      0),
      Cell_free_DNA = replace_na(Cell_free_DNA, 0),
      Total_DNA     = Cell_DNA + Cell_free_DNA
    )
}

trout_wide <- pivot_db(trout_avg)
mouse_wide <- pivot_db(mouse_avg)

# Combine for summary and plotting
all_wide <- bind_rows(
  mutate(trout_wide, Species = "Trout"),
  mutate(mouse_wide, Species = "Mouse")
)

# ------------------------------------------------------------------------
# Section 5: Summary table
# ------------------------------------------------------------------------

# Each row in all_wide is one replicate (PCR-duplicate-averaged).
# PCR_pair is dropped here as it has served its purpose; mean and SD
# are then computed across the 3 replicates per sample group.
summary_db <- all_wide %>%
  dplyr::select(-PCR_pair) %>%
  dplyr::group_by(`Sample Name`, Species, DNA_type, Route, Nominal_ng) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_ng = mean(Total_DNA, na.rm = TRUE),
    sd_ng = sd(Total_DNA, na.rm = TRUE),
    cv_pct = 100 * sd(Total_DNA, na.rm = TRUE) / mean(Total_DNA, na.rm = TRUE),
    recovery_pct = 100 * mean(Total_DNA, na.rm = TRUE) / dplyr::first(Nominal_ng),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Species, DNA_type, Route, dplyr::desc(Nominal_ng))

write.csv(summary_db, "./Results/DTBS_check_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------
# Section 6: Recovery plots – single-species data
# ------------------------------------------------------------------------

plot_recovery_single <- function(df_wide) {
  df <- df_wide %>%
    filter(DNA_type != "Mixed") %>%
    mutate(
      Panel = case_when(
        Species == "Trout" ~ "Cell-free",
        Species == "Mouse" ~ "Cellular"
      ),
      Panel = factor(Panel, levels = c("Cell-free", "Cellular")),
      Route = factor(Route, levels = c("Buffer", "Swab")),
      Nominal_label = factor(paste0(Nominal_ng, " ng"),
                             levels = paste0(sort(unique(Nominal_ng), decreasing = TRUE), " ng"))
    )
  
  ggplot(df, aes(x = Nominal_label, y = Total_DNA, fill = Route)) +
    geom_boxplot(position = position_dodge(width = 0.7), width = 0.55,
                 outlier.shape = NA, colour = "black", linewidth = 0.35) +
    geom_jitter(aes(group = Route),
                position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.08),
                size = 1.3, alpha = 0.7, colour = "black") +
    scale_fill_manual(values = c("Buffer" = "#C6DBEF", "Swab" = "#6BAED6"),
                      name = "Extraction route") +
    facet_wrap(~ Panel, scales = "free_x") +
    labs(x = "\nNominal input (ng)", y = "DNA recovered (ng)\n") +
    theme_bw() +
    theme(
      text = element_text(family = "Arial", size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12, colour = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

plot_single <- plot_recovery_single(all_wide)
plot_single

ggsave("./Results/DTBS_check_recovery_single_species.png", plot_single,
       width = 8, height = 5, dpi = 600)

# ------------------------------------------------------------------------
# Section 7: Mixed samples comparison plot
# ------------------------------------------------------------------------

plot_recovery_mixed <- function(df_wide) {
  df <- df_wide %>%
    filter(DNA_type == "Mixed") %>%
    mutate(
      Panel = case_when(
        Species == "Trout" ~ "Cell-free DNA recovered",
        Species == "Mouse" ~ "Cellular DNA recovered"
      ),
      Panel = factor(Panel, levels = c("Cell-free DNA recovered", "Cellular DNA recovered")),
      Route = factor(Route, levels = c("Buffer", "Swab")),
      Nominal_label = factor(paste0(Nominal_ng, " ng"),
                             levels = paste0(sort(unique(Nominal_ng), decreasing = TRUE), " ng"))
    )
  
  ggplot(df, aes(x = Nominal_label, y = Total_DNA, fill = Route)) +
    geom_boxplot(position = position_dodge(width = 0.7), width = 0.55,
                 outlier.shape = NA, colour = "black", linewidth = 0.35) +
    geom_jitter(aes(group = Route),
                position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.08),
                size = 1.3, alpha = 0.7, colour = "black") +
    scale_fill_manual(values = c("Buffer" = "#C6DBEF", "Swab" = "#6BAED6"),
                      name = "Extraction route") +
    facet_wrap(~ Panel, scales = "free_x") +
    labs(x = "\nNominal input per species (ng)", y = "DNA recovered (ng)\n") +
    theme_bw() +
    theme(
      text = element_text(family = "Arial", size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12, colour = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

plot_mixed <- plot_recovery_mixed(all_wide)
plot_mixed

ggsave("./Results/DTBS_check_recovery_mixed.png", plot_mixed,
       width = 8, height = 5, dpi = 600, units = "in")

# ------------------------------------------------------------------------
# Section 8: Recovery plots – % recovery (boxplots)
# ------------------------------------------------------------------------

# Calculates % recovery per replicate (Total_DNA / Nominal_ng * 100),
# so the boxplots show the spread across replicates.
# The horizontal dashed line marks 100% recovery.

plot_recovery_pct_single <- function(df_wide) {
  df <- df_wide %>%
    filter(DNA_type != "Mixed") %>%
    mutate(
      Panel = case_when(
        Species == "Trout" ~ "Cell-free",
        Species == "Mouse" ~ "Cellular"
      ),
      Panel = factor(Panel, levels = c("Cell-free", "Cellular")),
      Route = factor(Route, levels = c("Buffer", "Swab")),
      Nominal_label = factor(paste0(Nominal_ng, " ng"),
                             levels = paste0(sort(unique(Nominal_ng), decreasing = TRUE), " ng")),
      recovery_pct = Total_DNA / Nominal_ng * 100
    )
  
  ggplot(df, aes(x = Nominal_label, y = recovery_pct, fill = Route)) +
    geom_boxplot(position = position_dodge(0.7), width = 0.55,
                 outlier.shape = NA, colour = "black", linewidth = 0.35) +
    geom_jitter(aes(group = Route),
                position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.08),
                size = 1.3, alpha = 0.7, colour = "black") +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "grey40", linewidth = 0.45) +
    scale_fill_manual(values = c("Buffer" = "#C6DBEF", "Swab" = "#6BAED6"),
                      name = "Extraction route") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    facet_wrap(~ Panel, scales = "free_x") +
    labs(x = "\nNominal input (ng)", y = "Recovery (% of nominal input)\n") +
    theme_bw() +
    theme(
      text = element_text(family = "Arial", size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

plot_recovery_pct_mixed <- function(df_wide) {
  df <- df_wide %>%
    filter(DNA_type == "Mixed") %>%
    mutate(
      Panel = case_when(
        Species == "Trout" ~ "Cell-free DNA recovered",
        Species == "Mouse" ~ "Cellular DNA recovered"
      ),
      Panel = factor(Panel, levels = c("Cell-free DNA recovered", "Cellular DNA recovered")),
      Route = factor(Route, levels = c("Buffer", "Swab")),
      Nominal_label = factor(paste0(Nominal_ng, " ng"),
                             levels = paste0(sort(unique(Nominal_ng), decreasing = TRUE), " ng")),
      recovery_pct = Total_DNA / Nominal_ng * 100
    )
  
  ggplot(df, aes(x = Nominal_label, y = recovery_pct, fill = Route)) +
    geom_boxplot(position = position_dodge(0.7), width = 0.55,
                 outlier.shape = NA, colour = "black", linewidth = 0.35) +
    geom_jitter(aes(group = Route),
                position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.08),
                size = 1.3, alpha = 0.7, colour = "black") +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "grey40", linewidth = 0.45) +
    scale_fill_manual(values = c("Buffer" = "#C6DBEF", "Swab" = "#6BAED6"),
                      name = "Extraction route") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    facet_wrap(~ Panel, scales = "free_x") +
    labs(x = "\nNominal input per species (ng)", y = "Recovery (% of nominal input)\n") +
    theme_bw() +
    theme(
      text = element_text(family = "Arial", size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

plot_pct_single <- plot_recovery_pct_single(all_wide)
plot_pct_single

plot_pct_mixed <- plot_recovery_pct_mixed(all_wide)
plot_pct_mixed

ggsave("./Results/DTBS_check_recovery_pct_single.png", plot_pct_single,
       width = 8, height = 5, dpi = 600, units = "in")
ggsave("./Results/DTBS_check_recovery_pct_mixed.png", plot_pct_mixed,
       width = 8, height = 5, dpi = 600, units = "in")

# ------------------------------------------------------------------------
# Section 8: statistical analysis of DTB/DTS control data
# ------------------------------------------------------------------------

# Prepare replicate-level data
stats_df <- all_wide %>%
  mutate(
    Route = factor(Route, levels = c("Buffer", "Swab")),
    Recovery_pct = Total_DNA / Nominal_ng * 100
  ) %>%
  arrange(Species, DNA_type, desc(Nominal_ng), Route)

# ------------------------------------------------------------------------
# 8.1 Descriptive statistics by route
# ------------------------------------------------------------------------

db_descriptive_stats <- stats_df %>%
  group_by(Species, DNA_type, Nominal_ng, Route) %>%
  summarise(
    n = n(),
    mean_ng = mean(Total_DNA, na.rm = TRUE),
    sd_ng = sd(Total_DNA, na.rm = TRUE),
    median_ng = median(Total_DNA, na.rm = TRUE),
    min_ng = min(Total_DNA, na.rm = TRUE),
    max_ng = max(Total_DNA, na.rm = TRUE),
    mean_recovery_pct = mean(Recovery_pct, na.rm = TRUE),
    sd_recovery_pct = sd(Recovery_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(DNA_type, Species, desc(Nominal_ng), Route)

db_descriptive_stats

write.csv(
  db_descriptive_stats,
  "./Results/DTBS_check_descriptive_stats_by_route.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------
# 8.2 Statistical comparison between Buffer and Swab
# ------------------------------------------------------------------------

analyse_dtbs_condition <- function(df) {
  
  buffer_vals <- df %>%
    filter(Route == "Buffer") %>%
    pull(Total_DNA)
  
  swab_vals <- df %>%
    filter(Route == "Swab") %>%
    pull(Total_DNA)
  
  shapiro_buffer_p <- if (length(buffer_vals) >= 3) shapiro.test(buffer_vals)$p.value else NA_real_
  shapiro_swab_p   <- if (length(swab_vals) >= 3) shapiro.test(swab_vals)$p.value else NA_real_
  
  both_normal <- !is.na(shapiro_buffer_p) &&
    !is.na(shapiro_swab_p) &&
    shapiro_buffer_p > 0.05 &&
    shapiro_swab_p > 0.05
  
  if (both_normal) {
    variance_test <- "F test"
    variance_out <- var.test(buffer_vals, swab_vals)
    variance_p <- variance_out$p.value
    
    if (variance_p > 0.05) {
      location_test <- "Student t-test"
      test_out <- t.test(buffer_vals, swab_vals, var.equal = TRUE)
    } else {
      location_test <- "Welch t-test"
      test_out <- t.test(buffer_vals, swab_vals, var.equal = FALSE)
    }
    
  } else {
    variance_test <- NA_character_
    variance_p <- NA_real_
    location_test <- "Wilcoxon rank-sum test"
    test_out <- wilcox.test(buffer_vals, swab_vals, exact = FALSE)
  }
  
  tibble(
    n_buffer = length(buffer_vals),
    n_swab = length(swab_vals),
    mean_buffer = mean(buffer_vals, na.rm = TRUE),
    mean_swab = mean(swab_vals, na.rm = TRUE),
    sd_buffer = sd(buffer_vals, na.rm = TRUE),
    sd_swab = sd(swab_vals, na.rm = TRUE),
    shapiro_buffer_p = shapiro_buffer_p,
    shapiro_swab_p = shapiro_swab_p,
    variance_test = variance_test,
    variance_p = variance_p,
    location_test = location_test,
    location_p = test_out$p.value,
    note = "Exploratory analysis; n = 3 per route"
  )
}

db_stats_results <- stats_df %>%
  group_by(Species, DNA_type, Nominal_ng) %>%
  group_modify(~ analyse_dtbs_condition(.x)) %>%
  ungroup() %>%
  mutate(sd_ratio_buffer_to_swab = sd_buffer / sd_swab) %>%
  arrange(DNA_type, Species, desc(Nominal_ng))

db_stats_results

write.csv(
  db_stats_results,
  "./Results/DTBS_check_exploratory_stats.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------
# 8.3 Holm correction
# ------------------------------------------------------------------------

db_stats_results_adj <- db_stats_results %>%
  mutate(
    variance_p_holm = p.adjust(variance_p, method = "holm"),
    location_p_holm = p.adjust(location_p, method = "holm")
  )

db_stats_results_adj

write.csv(
  db_stats_results_adj,
  "./Results/DTBS_check_exploratory_stats_adjusted.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------
# 8.4 Compact table for manuscript / SI
# ------------------------------------------------------------------------

db_stats_compact <- db_stats_results_adj %>%
  transmute(
    Species,
    DNA_type,
    Nominal_ng,
    n_buffer,
    n_swab,
    Buffer_mean_sd = sprintf("%.2f ± %.2f", mean_buffer, sd_buffer),
    Swab_mean_sd   = sprintf("%.2f ± %.2f", mean_swab, sd_swab),
    Shapiro_p_Buffer = round(shapiro_buffer_p, 4),
    Shapiro_p_Swab   = round(shapiro_swab_p, 4),
    Variance_test = variance_test,
    Variance_p = round(variance_p, 4),
    Variance_p_Holm = round(variance_p_holm, 4),
    Location_test = location_test,
    Location_p = round(location_p, 4),
    Location_p_Holm = round(location_p_holm, 4),
    Note = note
  )

db_stats_compact

write.csv(
  db_stats_compact,
  "./Results/DTBS_check_exploratory_stats_compact.csv",
  row.names = FALSE
)