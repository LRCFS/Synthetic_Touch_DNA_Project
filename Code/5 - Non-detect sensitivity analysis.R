##########################################################################
#####           5 - Non-detect sensitivity analysis                  #####
##########################################################################

if (!exists("Data_with_non_detects")) {
  stop("Run Global.R after adding 'Data_with_non_detects <- Data' BEFORE the Ct filter in script 1.")
}

out_dir <- "Results/Non_detects"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------
# 1. Prepare qPCR data (same corrected names and volume factors as script 1)
# ------------------------------------------------------------------------
wells <- Data_with_non_detects %>%
  filter(Study %in% c("Poetsch", "Fonnelop", "Goray", "Meakin",
                      "Daly", "Bowman", "Thomasma", "Lim"),
         Task == "UNKNOWN", !is.na(`Sample Name corrected`)) %>%
  mutate(
    sample_id = as.character(`Sample Name corrected`),
    target = case_when(
      grepl("mouse", `Target Name`, ignore.case = TRUE) ~ "Cell",
      grepl("trout", `Target Name`, ignore.case = TRUE) ~ "Cell_free",
      TRUE ~ NA_character_
    ),
    ct_value = suppressWarnings(as.numeric(as.character(`Cт`))),
    # Detection uses the measured Ct and quantity, not a filled-in zero.
    detected_well = !is.na(ct_value) & !is.na(Quantity_Total) &
      Quantity_Total > 0
  ) %>%
  filter(!is.na(target))

# ------------------------------------------------------------------------
# 2. Average technical replicates per target, retaining non-detects
# ------------------------------------------------------------------------
# Same averaging rule as custom_average(): ignore undetected wells when at
# least one technical replicate is detectable. Both undetected -> 0.
targets <- wells %>%
  group_by(Study, sample_id, target) %>%
  summarise(
    n_wells = n(),
    n_detected_wells = sum(detected_well),
    amount_ng = if (any(detected_well)) {
      mean(Quantity_Total[detected_well])
    } else {
      0
    },
    .groups = "drop"
  )

# ------------------------------------------------------------------------
# 3. Combine mouse and trout into one row per experimental sample
# ------------------------------------------------------------------------

samples <- targets %>%
  pivot_wider(
    names_from = target,
    values_from = c(n_wells, n_detected_wells, amount_ng)
  ) %>%
  mutate(
    both_targets_present = !is.na(n_wells_Cell) & !is.na(n_wells_Cell_free),
    detected = coalesce(n_detected_wells_Cell, 0L) +
      coalesce(n_detected_wells_Cell_free, 0L) > 0,
    total_ng = if_else(
      both_targets_present,
      amount_ng_Cell + amount_ng_Cell_free,
      NA_real_
    )
  ) %>%
  separate(
    sample_id,
    into = c("Code", "Operator", "Condition", "Hand", "Repeat"),
    sep = "_",
    remove = FALSE,
    extra = "merge",
    fill = "right"
  ) %>%
  mutate(
    # Apply study-specific experimental conditions
    Condition = case_when(
      Study == "Bowman" & grepl("_3s$", sample_id) ~ "3s",
      Study == "Bowman" & grepl("_15s$", sample_id) ~ "15s",
      Study %in% c("Poetsch", "Fonnelop", "Goray") &
        Condition == "L" ~ "Low",
      Study %in% c("Poetsch", "Fonnelop", "Goray") &
        Condition == "A" ~ "Medium",
      Study %in% c("Poetsch", "Fonnelop", "Goray") &
        Condition == "H" ~ "High",
      TRUE ~ Condition
    ),
    
    # Standardise hand labels
    Hand = recode(
      Hand,
      "L" = "Left",
      "R" = "Right",
      "left" = "Left",
      "right" = "Right"
    )
  )

# Replace missing DNA quantities with 0
# and recalculate total DNA
samples <- samples %>%
  mutate(
    amount_ng_Cell = replace_na(amount_ng_Cell, 0),
    amount_ng_Cell_free = replace_na(amount_ng_Cell_free, 0),
    total_ng = amount_ng_Cell + amount_ng_Cell_free
  )

# ------------------------------------------------------------------------
# 4. Compare original detectable-only mean with zero-substitution scenario
# ------------------------------------------------------------------------
# Analysis is SAMPLE-level, not individual PCR-well-level.
# Incomplete target assays are excluded from BOTH quantitative summaries.
summarise_recovery <- function(df) {
  df %>%
    filter(both_targets_present) %>%
    summarise(
      n_samples = n(),
      n_detected = sum(detected),
      n_non_detect = sum(!detected),
      pct_non_detect = if (n() > 0) 100 * mean(!detected) else NA_real_,
      mean_detected_only = if (any(detected)) mean(total_ng[detected]) else NA_real_,
      sd_detected_only = if (sum(detected) > 1) sd(total_ng[detected]) else NA_real_,
      mean_including_zero = if (n() > 0) mean(total_ng) else NA_real_,
      sd_including_zero = if (n() > 1) sd(total_ng) else NA_real_,
      .groups = "drop"
    )
}

by_study <- samples %>% group_by(Study) %>% summarise_recovery()
by_condition <- samples %>% group_by(Study, Condition) %>% summarise_recovery()
by_condition_hand <- samples %>%
  group_by(Study, Condition, Hand) %>% summarise_recovery()

# ------------------------------------------------------------------------
# 5. Goray: left-right PAIRED Wilcoxon test under the two assumptions
# ------------------------------------------------------------------------
# Match pairs by Operator + concentration + repeat (same as script 4).
goray_pairs <- samples %>%
  filter(Study == "Goray", both_targets_present,
         Hand %in% c("Left", "Right")) %>%
  select(Operator, Condition, Repeat, Hand, total_ng, detected) %>%
  pivot_wider(names_from = Hand, values_from = c(total_ng, detected))

# Guard against duplicate pair identifiers, which invalidate pivoting.
if (anyDuplicated(samples %>%
                  filter(Study == "Goray", both_targets_present) %>%
                  select(Operator, Condition, Repeat, Hand))) {
  stop("Duplicate Goray pairing identifiers. Check sample names before testing.")
}

original_pairs <- goray_pairs %>%
  filter(!is.na(total_ng_Left), !is.na(total_ng_Right),
         detected_Left, detected_Right)
zero_pairs <- goray_pairs %>%
  filter(!is.na(total_ng_Left), !is.na(total_ng_Right))

paired_test <- function(df) {
  if (nrow(df) < 2 || all(df$total_ng_Left == df$total_ng_Right)) {
    return(NA_real_)
  }
  wilcox.test(df$total_ng_Left, df$total_ng_Right,
              paired = TRUE, exact = FALSE)$p.value
}

goray_wilcoxon <- data.frame(
  analysis = c("Both samples detected", "Include sample non-detects as zero"),
  n_pairs = c(nrow(original_pairs), nrow(zero_pairs)),
  p_value = c(paired_test(original_pairs), paired_test(zero_pairs))
)

# ------------------------------------------------------------------------
# 6. Export results to a separate folder
# ------------------------------------------------------------------------
write.csv(by_study, file.path(out_dir, "01_by_study.csv"), row.names = FALSE)
write.csv(by_condition, file.path(out_dir, "02_by_condition.csv"), row.names = FALSE)
write.csv(goray_wilcoxon, file.path(out_dir, "03_Goray_Wilcoxon.csv"), row.names = FALSE)

print(by_study)
print(by_condition %>% filter(Study == "Goray"))
print(goray_wilcoxon)