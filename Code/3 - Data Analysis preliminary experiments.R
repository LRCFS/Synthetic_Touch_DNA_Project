##########################################################################
#####                      Data analysis washes                      #####
##########################################################################

# ------------------------------------------------------------------------
# Overview:
# This R script analyses preliminary DNA quantity data from two washing-related
# experiments: hand washing and hand swabbing. It processes the cleaned datasets
# for each experiment, extracts the relevant experimental variables from sample
# names, averages replicate measurements, and combines cell-associated and
# cell-free DNA quantities to calculate total DNA per sample.
#
# The workflow includes:
# 1. Selecting the relevant columns from each dataset
# 2. Splitting sample names into experimental factors such as operator,
#    condition, hand, volume, and repeat
# 3. Averaging replicate measurements within each grouping level
# 4. Renaming target categories into cell DNA and cell-free DNA
# 5. Calculating total DNA quantity per sample
# 6. Assessing data distribution and applying group comparisons
# 7. Producing summary tables and boxplots for visualisation
#
# Output: Summary .csv files and comparative plots for the hand washing and
# hand swab datasets.
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Section 1: Hand washing
# ------------------------------------------------------------------------
# Select the columns of interest
Data_preliminary_HWashing <- Data_preliminary_HWashing %>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity_Total")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Data_preliminary_HWashing_intermediate <- Data_preliminary_HWashing %>%
  separate(`Sample Name`, into = c("Operator","Condition", "Hand", "Repeat"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Operator)) %>%
  select(Operator, Condition, Hand, Repeat, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Data_preliminary_HWashing_averaged <- Data_preliminary_HWashing_intermediate %>%
  group_by(Operator, Condition, Hand, `Target Name`, Repeat) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Rename Target Name values
Data_preliminary_HWashing_renamed <- Data_preliminary_HWashing_averaged %>%
  mutate(`Target Name` = recode(`Target Name`,
                                "MOUSE" = "Cell_DNA",
                                "TROUT" = "Cell_free_DNA"))

# Calculate total DNA quantity per row per hand
Data_preliminary_HWashing_with_total <- Data_preliminary_HWashing_renamed %>%
  pivot_wider(names_from = `Target Name`, values_from = Total_Quantity) %>%
  mutate(
    Cell_DNA = replace_na(Cell_DNA, 0),
    Cell_free_DNA = replace_na(Cell_free_DNA, 0),
    Total_DNA = Cell_DNA + Cell_free_DNA
  )

# Check normality (per Operator × Condition)
# If all p-values > 0.05 → data are approximately normal → use t-test
# If any < 0.05 → not normal → use Wilcoxon test
Data_preliminary_HWashing_with_total %>%
  group_by(Operator, Condition) %>%
  summarise(
    p_shapiro = shapiro.test(Total_DNA)$p.value,
    .groups = "drop"
  )

# If normal (t-test):
t_test_results <- Data_preliminary_HWashing_with_total %>%
  group_by(Operator) %>%
  summarise(
    t_test = list(t.test(Total_DNA ~ Condition)),
    results = list(tidy(t_test[[1]]))
  ) %>%
  tidyr::unnest(results)
t_test_results

# Extract p-values per operator
p_values <- t_test_results %>%
  select(Operator, p.value)

# Summarise data
summary_df <- Data_preliminary_HWashing_with_total %>%
  group_by(Operator, Condition) %>%
  summarise(
    n = n(),
    mean_DNA = mean(Total_DNA, na.rm = TRUE),
    se_DNA = sd(Total_DNA, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )
summary_df

# Export to CSV
write.csv(summary_df, "./Results/summary_HWashing_TotalDNA.csv", row.names = FALSE)

# Boxplot + t-test significance labels
# Rename operators
Data_preliminary_HWashing_with_total <- Data_preliminary_HWashing_with_total %>%
  mutate(Operator = recode(Operator,
                           "Op1" = "Operator 1",
                           "Op2" = "Operator 2"))

p <- ggplot(Data_preliminary_HWashing_with_total, aes(x = Condition, y = Total_DNA, fill = Condition)) +
  geom_boxplot(alpha = 0.9,outlier.shape = NA,colour = "black") +
  scale_fill_manual(values = c("unwashed" = "grey80", "washed" = "#C6DBEF"),  # same blue for both
    name = "Condition") +
  facet_wrap(~Operator) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    comparisons = list(c("washed", "unwashed"))) +
  labs(
    title = "",
    x = "Condition",
    y = "Total DNA (ng)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )
p

# ------------------------------------------------------------------------
# Section 2: Hand swab
# ------------------------------------------------------------------------
# Select the columns of interest
Data_preliminary_HSwab <- Data_preliminary_HSwab %>%
  select(`Sample Name`=`Sample Name corrected`,`Target Name`, "Quantity_Total")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Data_preliminary_HSwab_intermediate <- Data_preliminary_HSwab %>%
  separate(`Sample Name`, into = c("Volume", "Operator","Hand", "Repeat"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Volume)) %>%
  select(Volume, Operator, Hand, Repeat, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Data_preliminary_HSwab_averaged <- Data_preliminary_HSwab_intermediate %>%
  group_by(Volume, Operator, Hand, Repeat, `Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Rename Target Name values
Data_preliminary_HSwab_renamed <- Data_preliminary_HSwab_averaged %>%
  mutate(`Target Name` = recode(`Target Name`,
                                "mouse" = "Cell_DNA",
                                "TROUT 2" = "Cell_free_DNA"))

# Calculate total DNA quantity per row per hand
Data_preliminary_HSwab_with_total <- Data_preliminary_HSwab_renamed %>%
  pivot_wider(names_from = `Target Name`, values_from = Total_Quantity) %>%
  mutate(
    Cell_DNA = replace_na(Cell_DNA, 0),
    Cell_free_DNA = replace_na(Cell_free_DNA, 0),
    Total_DNA = Cell_DNA + Cell_free_DNA
  )

# Calculate ratio between cfDNA and cDNA
Data_preliminary_HSwab_with_total$Ratio <- Data_preliminary_HSwab_with_total$Cell_DNA/Data_preliminary_HSwab_with_total$Cell_free_DNA
mean(Data_preliminary_HSwab_with_total$Ratio[
  is.finite(Data_preliminary_HSwab_with_total$Ratio)
])

# Check normality (per Operator × Hand)
# If all p-values > 0.05 → data are approximately normal → use t-test
# If any < 0.05 → not normal → use Wilcoxon test
Data_preliminary_HSwab_with_total %>%
  group_by(Operator, Hand) %>%
  summarise(
    p_shapiro = shapiro.test(Total_DNA)$p.value,
    .groups = "drop"
  )

# If not normal (Wilcoxon test):
wilcox_results <- Data_preliminary_HSwab_with_total %>%
  group_by(Operator) %>%
  summarise(
    wilcox_test = list(wilcox.test(Total_DNA ~ Hand)),
    results = list(tidy(wilcox_test[[1]]))
  ) %>%
  tidyr::unnest(results)
wilcox_results

# Extract p-values per operator
p_values <- wilcox_results %>%
  select(Operator, p.value)

# Summarise data
summary_df <- Data_preliminary_HSwab_with_total %>%
  group_by(Operator, Hand) %>%
  summarise(
    n = n(),
    mean_DNA = mean(Total_DNA, na.rm = TRUE),
    sd_DNA = sd(Total_DNA, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
summary_df

# Export to CSV
write.csv(summary_df, "./Results/summary_HSwab_with_total.csv", row.names = FALSE)


# Boxplot + t-test significance labels
# Rename operators
Data_preliminary_HSwab_with_total <- Data_preliminary_HSwab_with_total %>%
  mutate(Operator = recode(Operator,
                           "Op1" = "Operator 1",
                           "Op2" = "Operator 2",
                           "Op3" = "Operator 3"))

p <- ggplot(Data_preliminary_HSwab_with_total, aes(x = Hand, y = Total_DNA, fill = Hand)) +
  geom_boxplot(alpha = 0.9,outlier.shape = NA,colour = "black") +
  scale_fill_manual(values = c("Left" = "grey80", "Right" = "#C6DBEF"),  # same blue for both
                    name = "Hand") +
  facet_wrap(~Operator) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    comparisons = list(c("Left", "Right"))) +
  labs(
    title = "",
    x = "Hand",
    y = "Total DNA (ng)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )
p

