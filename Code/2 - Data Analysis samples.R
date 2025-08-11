##########################################################################
#####                      Data analysis washes                      #####
##########################################################################

# ------------------------------------------------------------------------
# Overview:
# This R script analyses DNA quantity data from multiple published studies 
# and from repeated experimental data collected in the current project.
# For each study, the workflow:
#   1. Imports DNA quantification data from the original publication
#   2. Cleans and reformats the repeated experimental data (proxy experiments)
#   3. Calculates total DNA quantities
#   4. Produces visual comparisons between published results and repeated data
#   5. Exports processed datasets for supplementary information
#   6. Combines plots into publication-ready figures
#
# Included studies and their corresponding sections:
#   Section 1 – Poetsch et al.
#   Section 2 – Fonnelop et al.
#   Section 3 – Goray et al.
#   Section 4 – Meakin et al.
#   Section 5 – Daly et al.
#   Section 6 – Bowman et al.
#   Section 7 – Thomasma and Foran
#   Section 8 – Lim et al.
#
# Output: Cleaned datasets (.csv) and high-resolution combined figures (.png)
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Section 1: Poetsch et al.
# ------------------------------------------------------------------------
#### Data from the article ####
# Load the Excel file
Poetsch_data <- read.table("./Data/Approximate_Poetsch_dna_quantities.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Reorder factor levels for Age_Group
Poetsch_data$Age_Group <- factor(Poetsch_data$Age_Group, levels = c(
  "1-4 y", "5-10 y", "11-14 y", "15-20 y", "21-24 y", "25-30 y",
  "31-40 y", "41-50 y", "51-60 y", "61-70 y", "71-80 y", "> 80 y"
))

# Plot
plot_Poetsch_data <- ggplot(Poetsch_data, aes(x = Age_Group, y = DNA_qt)) +
  geom_bar(stat = "identity", fill = "#C6DBEF", colour = "black", width = 0.7) +
  geom_errorbar(aes(ymin = DNA_qt, ymax = DNA_qt + Error), width = 0.3) +
  scale_y_continuous(
    limits = c(0, 5),
    breaks = seq(0, 5, by = 0.5),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  labs(x = "Age Group",y = "DNA Quantity (ng)") +
  theme(text = element_text(family = "Arial", size = 14),  # Set overall font
        axis.title = element_text(size = 14),              # Axis title font size
        axis.text = element_text(size = 12),               # Axis tick labels
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

# Show plot
plot_Poetsch_data

#### Repeated data with Proxy ####
# Remove rows with unused data
Poetsch_repeat <- Poetsch_repeat %>%
  filter(`Sample Name` != "LH" & `Sample Name` != "LH_1" & `Sample Name` != "LH_2"
         & `Sample Name` != "RH" & `Sample Name` != "RH_1" & `Sample Name` != "RH_2"
         & `Sample Name` != "C")

# Select the columns of interest
Poetsch_repeat <- Poetsch_repeat %>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Poetsch_data_intermediate <- Poetsch_repeat %>%
  separate(`Sample Name`, into = c("Study", "Operator","Concentration", "Hand", "Repeat"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study)) %>%
  select(Study, Operator, Concentration, Hand,Repeat, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Poetsch_data_averaged <- Poetsch_data_intermediate %>%
  group_by(Study, Operator, Concentration, Hand, `Target Name`, Repeat) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Rename Target Name values
Poetsch_data_renamed <- Poetsch_data_averaged %>%
  mutate(`Target Name` = recode(`Target Name`,
                                "Mouse" = "Cell_DNA",
                                "TROUT 2" = "Cell_free_DNA"))

# Rename Concentration values
Poetsch_data_renamed <- Poetsch_data_renamed %>%
  mutate(Concentration = recode(Concentration,
                                "L" = "Low",
                                "A" = "Average",
                                "H" = "High",))

# Calculate total DNA quantity per row per hand
Poetsch_data_with_total <- Poetsch_data_renamed %>%
  pivot_wider(names_from = `Target Name`, values_from = Total_Quantity) %>%
  mutate(
    Cell_DNA = replace_na(Cell_DNA, 0),
    Cell_free_DNA = replace_na(Cell_free_DNA, 0),
    Total_DNA = Cell_DNA + Cell_free_DNA
  )

# force the x-axis (Concentration) to appear in the order Low, Average, High instead of alphabetical (A, H, L)
Poetsch_data_with_total <- Poetsch_data_with_total %>%
  mutate(Concentration = factor(Concentration, levels = c("Low", "Average", "High")))

# Plot repeated data as a boxplot
plot_Poetsch_repeat <- Poetsch_data_with_total %>%
  ggplot(aes(x = Concentration, y = Total_DNA)) +
  geom_boxplot(fill = "#C6DBEF", colour = "black", width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  labs(y = "DNA Quantity (ng)", x = "Input concentration") +
  scale_y_continuous(
    limits = c(0, 5),
    breaks = seq(0, 5, by = 0.5),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

plot_Poetsch_repeat

# export data for supplementary information
write.csv(Poetsch_data_with_total, "./Results/Poetsch_repeat_data.csv", row.names = FALSE)

#### Final plot for article - Figure 1 ####
# Clean individual plots
plot_A <- plot_Poetsch_data + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

plot_B <- plot_Poetsch_repeat + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots with labels A and B
pCombined_Poetsch_pending <- ggarrange(plot_A, plot_B,
                                       labels = c("A", "B"),
                                       ncol = 2, nrow = 1,
                                       align = "hv",
                                       common.legend = FALSE,
                                       font.label = list(size = 12, color = "black"),
                                       hjust = -0.5, vjust = 1.2) +
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"))  # (Top, Right, Bottom, Left)

# Add shared y-axis label
pCombined_Poetsch <- annotate_figure(
  pCombined_Poetsch_pending,
  left = text_grob("DNA Quantity (ng)", rot = 90, vjust = 0.5, hjust = 0.5, size = 12),
  bottom = NULL  # You can add shared x-axis here if you want
)

# Show the combined plot
pCombined_Poetsch

# Save the figure
ggsave("./Results/Poetsch_combined_plot.png", pCombined_Poetsch,
       width = 10, height = 5, dpi = 600, units = "in")

# ------------------------------------------------------------------------
# Section 2: Fonnelop et al.
# ------------------------------------------------------------------------
#### Data from the article ####
fonnelop_data <- read_csv("./Data/Approximate_fonnelop_dna_quantities.csv") %>%
  mutate(Study = "Fonnelop et al.") %>%
  rename(Total_DNA = Quantity)

# Plot
plot_fonnelop_data <- ggplot(fonnelop_data, aes(x = "Fonnelop <i>et al.</i>", y = Total_DNA)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, colour = "black", fill = NA) +
  geom_jitter(aes(colour = Gender, shape = Gender), width = 0.15, size = 2, alpha = 0.7) +
  scale_colour_manual(values = c("F" = "grey40", "M" = "steelblue")) +
  scale_shape_manual(values = c("F" = 16, "M" = 4)) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 0.5),
    limits = c(0, 4)
  ) +
  labs(x = "\nStudy", y = "DNA Quantity (ng)") +
  theme_bw() +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", colour = "black", size = 0.7),
    legend.margin = margin(7,9,7,8),
    legend.key.size = unit(1.2, "lines"),
    axis.text.x = element_markdown(family = "Arial", size = 12),
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14)
  )

# Show plot
plot_fonnelop_data

#### Repeated data with Proxy ####
# Remove rows with unused data
Fonnelop_repeat <- Fonnelop_repeat %>%
  filter(`Sample Name` != "LH_3" & `Sample Name` != "RH_3" & `Sample Name` != "C_2" & `Sample Name` != "C")

# Select the columns of interest
Fonnelop_repeat <- Fonnelop_repeat %>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Fonnelop_data_intermediate <- Fonnelop_repeat %>%
  separate(`Sample Name`, into = c("Study", "Operator","Concentration", "Hand", "Repeat"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study)) %>%
  select(Study, Operator, Concentration, Hand,Repeat, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Fonnelop_data_averaged <- Fonnelop_data_intermediate %>%
  group_by(Study, Operator, Concentration, Hand, `Target Name`, Repeat) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Rename Target Name values
Fonnelop_data_renamed <- Fonnelop_data_averaged %>%
  mutate(`Target Name` = recode(`Target Name`,
                                "Mouse" = "Cell_DNA",
                                "TROUT 2" = "Cell_free_DNA"))

# Rename Concentration values
Fonnelop_data_renamed <- Fonnelop_data_renamed %>%
  mutate(Concentration = recode(Concentration,
                                "L" = "Low",
                                "A" = "Average",
                                "H" = "High",))

# Calculate total DNA quantity per row per hand
Fonnelop_data_with_total <- Fonnelop_data_renamed %>%
  pivot_wider(names_from = `Target Name`, values_from = Total_Quantity) %>%
  mutate(
    Cell_DNA = replace_na(Cell_DNA, 0),
    Cell_free_DNA = replace_na(Cell_free_DNA, 0),
    Total_DNA = Cell_DNA + Cell_free_DNA
  )

# force the x-axis (Concentration) to appear in the order Low, Average, High instead of alphabetical (A, H, L)
Fonnelop_data_with_total <- Fonnelop_data_with_total %>%
  mutate(Concentration = factor(Concentration, levels = c("Low", "Average", "High")))

# Plot repeated data
plot_fonnelop_repeat <- ggplot(Fonnelop_data_with_total, aes(x = Concentration, y = Total_DNA, fill = Concentration)) +
  geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8, shape = 21, stroke = 0.3, colour = "black") +
  scale_fill_manual(values = c("Low" = "#C6DBEF", "Average" = "#6BAED6", "High" = "#2171B5")) +
  scale_y_continuous(
    limits = c(0, 4),
    breaks = seq(0, 4, by = 0.5)
  ) +
  labs(x = "\nInput concentration", y = "DNA Quantity (ng)") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Show plot
plot_fonnelop_repeat

# export data for supplementary information
write.csv(Fonnelop_data_with_total, "./Results/Fonnelop_repeat_data.csv", row.names = FALSE)

#### Final plot for article - Figure 2 ####
# Clean individual plots
plot_A <- plot_fonnelop_data + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

plot_B <- plot_fonnelop_repeat + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots with labels A and B
pCombined_Fonnelop_pending <- ggarrange(plot_A, plot_B,
                                       labels = c("A", "B"),
                                       ncol = 2, nrow = 1,
                                       align = "hv",
                                       common.legend = FALSE,
                                       font.label = list(size = 12, color = "black"),
                                       hjust = -0.5, vjust = 1.2) +
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"))  # (Top, Right, Bottom, Left)

# Add shared y-axis label
pCombined_Fonnelop <- annotate_figure(
  pCombined_Fonnelop_pending,
  left = text_grob("DNA Quantity (ng)", rot = 90, vjust = 0.5, hjust = 0.5, size = 12),
  bottom = NULL  # You can add shared x-axis here if you want
)

# Show the combined plot
pCombined_Fonnelop

# Save the figure
ggsave("./Results/Fonnelop_combined_plot.png", pCombined_Fonnelop,
       width = 10, height = 5, dpi = 600, units = "in")

# ------------------------------------------------------------------------
# Section 3: Goray et al.
# ------------------------------------------------------------------------
#### Data from the article ####
Goray_data <- read_csv("./Data/Goray_tidy_data.csv") %>%
  mutate(Study = "Goray et al.") %>%
  rename(Total_DNA = Quantity)

# Convert Hand to full labels for clarity
Goray_data$Hand <- recode(Goray_data$Hand, "L" = "Left", "R" = "Right")

# Plot
plot_Goray_data <- ggplot(Goray_data, aes(x = Hand, y = Total_DNA)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, fill = "white", colour = "black") +  # now first
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6, colour = "grey30") +            # now second
  scale_y_continuous(
    breaks = seq(0, 10, by = 2),
    limits = c(0, 10)) +
  labs(x = "\nHand",y = "DNA Quantity (ng)") +
  theme_bw() +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", colour = "black", size = 0.7),
    legend.margin = margin(7,9,7,8),
    legend.key.size = unit(1.2, "lines"),
    axis.text.x = element_markdown(family = "Arial", size = 12),
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14)
  )

# Show plot
plot_Goray_data

#### Repeated data with Proxy ####
# Select the columns of interest
Goray_repeat <- Goray_repeat %>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Goray_data_intermediate <- Goray_repeat %>%
  separate(`Sample Name`, into = c("Study", "Operator","Concentration", "Hand", "Repeat"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study)) %>%
  select(Study, Operator, Concentration, Hand,Repeat, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Goray_data_averaged <- Goray_data_intermediate %>%
  group_by(Study, Operator, Concentration, Hand, `Target Name`, Repeat) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Convert Hand to full labels for clarity
Goray_data_averaged$Hand <- recode(Goray_data_averaged$Hand, "L" = "Left", "R" = "Right")

# Rename Target Name values
Goray_data_renamed <- Goray_data_averaged %>%
  mutate(`Target Name` = recode(`Target Name`,
                                "Mouse" = "Cell_DNA",
                                "TROUT 2" = "Cell_free_DNA"))

# Rename Concentration values
Goray_data_renamed <- Goray_data_renamed %>%
  mutate(Concentration = recode(Concentration,
                                "L" = "Low",
                                "A" = "Average",
                                "H" = "High",))

# Calculate total DNA quantity per row per hand
Goray_data_with_total <- Goray_data_renamed %>%
  pivot_wider(names_from = `Target Name`, values_from = Total_Quantity) %>%
  mutate(
    Cell_DNA = replace_na(Cell_DNA, 0),
    Cell_free_DNA = replace_na(Cell_free_DNA, 0),
    Total_DNA = Cell_DNA + Cell_free_DNA
  )

# force the x-axis (Concentration) to appear in the order Low, Average, High instead of alphabetical (A, H, L)
Goray_data_with_total <- Goray_data_with_total %>%
  mutate(Concentration = factor(Concentration, levels = c("Low", "Average", "High")))

# Plot repeated data
plot_Goray_repeat <- ggplot(Goray_data_with_total, aes(x = Hand, y = Total_DNA)) +
  geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6, colour = "grey30") +
  scale_y_continuous(
    breaks = seq(0, 10, by = 2),
    limits = c(0, 10)) +
  labs(x = "\nHand", y = "DNA Quantity (ng)") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Show plot
plot_Goray_repeat

# export data for supplementary information
write.csv(Goray_data_with_total, "./Results/Goray_repeat_data.csv", row.names = FALSE)

#### Statistic analysis ####
#Data from the article
goray_paired <- Goray_data %>%
  pivot_wider(names_from = Hand, values_from = Total_DNA) %>%
  drop_na()  # Remove incomplete pairs

#Check normality of differences
shapiro.test(goray_paired$Left - goray_paired$Right)
# p-value < 2.2e-16
# Shapiro-Wilk test indicates that the differences between left and right hand DNA quantities are not normally distributed.

#Perform Wilcoxon signed-rank test (non-parametric, for paired data)
wilcox.test(goray_paired$Left, goray_paired$Right, paired = TRUE)
# p-value = 0.01668
# There is a statistically significant difference in DNA quantities recovered from the left and right hands in the Goray dataset (p < 0.05).

# Repeated data with proxy
# Summarise per participant, concentration, and repeat
paired_data <- Goray_data_with_total %>%
  select(Operator, Concentration, Repeat, Hand, Total_DNA) %>%
  pivot_wider(names_from = Hand, values_from = Total_DNA) %>%
  drop_na()

# Check normality of differences
shapiro.test(paired_data$Left - paired_data$Right)
# p-value = 1.127e-08
# Shapiro-Wilk test indicates that the differences between left and right hand DNA quantities are not normally distributed.

# If not normal → use Wilcoxon
wilcox.test(paired_data$Left, paired_data$Right, paired = TRUE)
# p-value is greater than 0.05, there is no statistically significant difference in DNA quantity recovered between left and right hands.

#### Final plot for article - Figure 3 ####
# Clean individual plots
plot_A <- plot_Goray_data + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

plot_B <- plot_Goray_repeat + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots with labels A and B
pCombined_Goray_pending <- ggarrange(plot_A, plot_B,
                                        labels = c("A", "B"),
                                        ncol = 2, nrow = 1,
                                        align = "hv",
                                        common.legend = FALSE,
                                        font.label = list(size = 12, color = "black"),
                                        hjust = -0.5, vjust = 1.2) +
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"))  # (Top, Right, Bottom, Left)

# Add shared y-axis label
pCombined_Goray <- annotate_figure(
  pCombined_Goray_pending,
  left = text_grob("DNA Quantity (ng)", rot = 90, vjust = 0.5, hjust = 0.5, size = 12),
  bottom = NULL  # You can add shared x-axis here if you want
)

# Show the combined plot
pCombined_Goray

# Save the figure
ggsave("./Results/Goray_combined_plot.png", pCombined_Goray,
       width = 10, height = 5, dpi = 600, units = "in")

# ------------------------------------------------------------------------
# Section 4: Meakin et al.
# ------------------------------------------------------------------------
#### Data from the article ####
Meakin_data <- read_csv("./Data/Approximate_Meakin_dna_quantities.csv")

# Plot
plot_Meakin_data <- ggplot(Meakin_data, aes(x = Participants, y = Total_DNA)) +
  geom_bar(stat = "identity", fill = "grey90", colour = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Total_DNA, ymax = Total_DNA + SD),
                width = 0.2) +
  scale_y_continuous(
    limits = c(0, 15),
    breaks = seq(0, 15, by = 2),
    expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "\nParticipants",
    y = "DNA Quantity (ng)"
  ) +
  theme_bw()+
  theme(
    text = element_text(family = "Arial", size = 14),  # Set overall font
    axis.title = element_text(size = 14),              # Axis title font size
    axis.text = element_text(size = 12),               # Axis tick labels
  )

# Show plot
plot_Meakin_data

#### Repeated data with Proxy ####
# Select the columns of interest
Meakin_repeat <- Meakin_repeat%>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity", "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Meakin_data_intermediate <- Meakin_repeat %>%
  separate(`Sample Name`, into = c("Study", "Knife","Operator"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study)) %>%
  select(Study, Operator, Knife, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Meakin_data_averaged <- Meakin_data_intermediate %>%
  group_by(Study, Operator, Knife,`Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Rename Target Name values
Meakin_data_renamed <- Meakin_data_averaged %>%
  mutate(`Target Name` = recode(`Target Name`,
                                "Mouse" = "Cell_DNA",
                                "TROUT 2" = "Cell_free_DNA"))

# Calculate total DNA quantity 
Meakin_data_with_total <- Meakin_data_renamed %>%
  pivot_wider(names_from = `Target Name`, values_from = Total_Quantity) %>%
  mutate(
    Cell_DNA = replace_na(Cell_DNA, 0),
    Cell_free_DNA = replace_na(Cell_free_DNA, 0),
    Total_DNA = Cell_DNA + Cell_free_DNA
  )

# Plot repeated data as a boxplot
plot_Meakin_repeat <- Meakin_data_with_total %>%
  ggplot(aes(x = Operator, y = Total_DNA)) +
  geom_boxplot(fill = "grey90", colour = "black", width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  labs(y = "DNA Quantity (ng)", x = "\nParticipants") +
  scale_y_continuous(
    limits = c(0, 15),
    breaks = seq(0, 15, by = 2),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Show plot
plot_Meakin_repeat

# export data for supplementary information
write.csv(Meakin_data_with_total, "./Results/Meakin_repeat_data.csv", row.names = FALSE)

#### Final plot for article - Figure 4 ####
#average SD
meanSD_Meakin_data <- mean(Meakin_data$SD);meanSD_Meakin_data
#average SD
sd_per_operator <- Meakin_data_with_total %>%
  group_by(Operator) %>%
  summarise(SD_Total_DNA = sd(Total_DNA, na.rm = TRUE))
meanSD_Meakin_repeat <- mean(sd_per_operator$SD_Total_DNA, na.rm = TRUE);meanSD_Meakin_repeat

# Clean individual plots
plot_A <- plot_Meakin_data + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

plot_B <- plot_Meakin_repeat + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots with labels A and B
pCombined_Meakin_pending <- ggarrange(plot_A, plot_B,
                                     labels = c("A", "B"),
                                     ncol = 2, nrow = 1,
                                     align = "hv",
                                     common.legend = FALSE,
                                     font.label = list(size = 12, color = "black"),
                                     hjust = -0.5, vjust = 1.2) +
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"))  # (Top, Right, Bottom, Left)

# Add shared y-axis label
pCombined_Meakin <- annotate_figure(
  pCombined_Meakin_pending,
  left = text_grob("DNA Quantity (ng)", rot = 90, vjust = 0.5, hjust = 0.5, size = 12),
  bottom = NULL  # You can add shared x-axis here if you want
)

# Show the combined plot
pCombined_Meakin

# Save the figure
ggsave("./Results/Meakin_combined_plot.png", pCombined_Meakin,
       width = 10, height = 5, dpi = 600, units = "in")
# ------------------------------------------------------------------------
# Section 5: Daly et al.
# ------------------------------------------------------------------------
#### Data from the article ####
# Create the data frame directly
Daly_data <- tibble(
  Material = "Fabric",
  Mean_DNA_ng = 1.23,
  Range_ng = 14.8 - 1.23  # Using the range: max - mean
)

# Add max value column
Daly_data <- Daly_data %>%
  mutate(Max_DNA = Mean_DNA_ng + Range_ng)

# Plot
plot_Daly_data <- ggplot(Daly_data, aes(x = Material, y = Mean_DNA_ng)) +
  geom_bar(stat = "identity", fill = "grey80", colour = "black", width = 0.6) +
  
  # Add a range line beside the bar (e.g. slightly to the right)
  geom_segment(aes(x = 1.4, xend = 1.4, y = 0, yend = Max_DNA),
               colour = "black", linewidth = 1, linetype = "dotted") +
  
  # Optional: add a top cap at max
  geom_point(aes(x = 1.4, y = Max_DNA), shape = 95, size = 5) +
  
  scale_y_continuous(
    limits = c(0, 15),
    breaks = seq(0, 15, by = 0.5),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = expression(paste("Daly ", italic("et al."))),
    y = "DNA Quantity (ng)"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
plot_Daly_data

y <- plot_Daly_data+scale_y_break(c(2, 13))+
  theme(axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_blank())
y

#### Repeated data with Proxy ####
# Select the columns of interest
Daly_repeat <- Daly_repeat%>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity", "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Daly_data_intermediate <- Daly_repeat %>%
  separate(`Sample Name`, into = c("Study", "Operator","Textile", "Hand", "Repeat"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study)) %>%
  select(Study, Operator, Textile, Hand,Repeat, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Daly_data_averaged <- Daly_data_intermediate %>%
  group_by(Study, Operator, Textile, Hand,Repeat, `Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Rename Target Name values
Daly_data_renamed <- Daly_data_averaged %>%
  mutate(`Target Name` = recode(`Target Name`,
                                "Mouse" = "Cell_DNA",
                                "TROUT 2" = "Cell_free_DNA"))


# Calculate total DNA quantity 
Daly_data_with_total <- Daly_data_renamed %>%
  pivot_wider(names_from = `Target Name`, values_from = Total_Quantity) %>%
  mutate(
    Cell_DNA = replace_na(Cell_DNA, 0),
    Cell_free_DNA = replace_na(Cell_free_DNA, 0),
    Total_DNA = Cell_DNA + Cell_free_DNA
  )

# Clean textile labels
Daly_data_with_total <- Daly_data_with_total %>%
  mutate(Textile_Group = recode(Textile,
                                "Grey-blue" = "Grey/Blue, woven",
                                "Yellow" = "Yellow, knitted"))

# Get max values per group
range_data <- Daly_data_with_total %>%
  group_by(Textile_Group) %>%
  summarise(Max_DNA = max(Total_DNA, na.rm = TRUE)) %>%
  mutate(x_pos = as.numeric(factor(Textile_Group)) + 0.5)

# Plot repeated data as a boxplot
plot_Daly_repeat <- Daly_data_with_total %>%
  ggplot(aes(x = Textile_Group, y = Total_DNA, fill = Textile_Group)) +
  geom_boxplot(colour = "black", width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  
  # Add side range lines
  geom_segment(data = range_data,
               aes(x = x_pos, xend = x_pos, y = 0, yend = Max_DNA),
               inherit.aes = FALSE,
               linetype = "dotted", linewidth = 1) +
  
  geom_point(data = range_data,
             aes(x = x_pos, y = Max_DNA),
             inherit.aes = FALSE,
             shape = 95, size = 5) +
  
  scale_fill_manual(values = c("Grey/Blue, woven" = "#377eb8", "Yellow, knitted" = "#FFD700")) +
  labs(y = "DNA Quantity (ng)", x = expression(paste("Current study", italic("")))) +
  scale_y_continuous(
    limits = c(0, 15),
    breaks = seq(0, 15, by = 0.5),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

x <- plot_Daly_repeat+scale_y_break(c(2, 13))+
  theme(axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_blank())
x

# Show plot
plot_Daly_repeat

# export data for supplementary information
write.csv(Daly_data_with_total, "./Results/Daly_repeat_data.csv", row.names = FALSE)

#### Final plot for article - Figure 5 ####
# Clean individual plots
plot_A <- y + 
  theme(plot.title = element_text(hjust = 0.5))
plot_A

plot_B <- x + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))
plot_B

# Combine plots with labels A and B
pCombined_Daly_pending <- (plot_A | plot_B) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

pCombined_Daly_pending

# Save the figure
ggsave("./Results/Daly_combined_plot.png", pCombined_Daly_pending,
       width = 10, height = 5, dpi = 600, units = "in")

# ------------------------------------------------------------------------
# Section 6: Bowman et al.
# ------------------------------------------------------------------------
#### Data from the article ####
Bowman_data <- read_csv("./Data/Bowman_dna_quantities.csv")

# Plot
plot_Bowman_data <- ggplot(Bowman_data, aes(x = Scenario, y = DNA_Quantity_ng)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, fill = "white", colour = "black") +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6, colour = "grey30") +
  scale_y_continuous(
    breaks = seq(0, 6, by = 2),
    limits = c(0, 6)
  ) +
  labs(
    x = "\nScenario",
    y = "DNA Quantity (ng)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_markdown(family = "Arial", size = 12),
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )

# Show plot
plot_Bowman_data

#### Repeated data with Proxy ####
# Select the columns of interest
Bowman_repeat <- Bowman_repeat %>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Bowman_data_intermediate <- Bowman_repeat %>%
  separate(`Sample Name`, into = c("Study", "Operator","Hand", "Repeat", "Contact"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study)) %>%
  select(Study, Operator, Contact, Hand, Repeat, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Bowman_data_averaged <- Bowman_data_intermediate %>%
  group_by(Study, Operator, Contact, Hand, Repeat, `Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Rename Target Name values
Bowman_data_renamed <- Bowman_data_averaged %>%
  mutate(`Target Name` = recode(`Target Name`,
                                "Mouse" = "Cell_DNA",
                                "TROUT 2" = "Cell_free_DNA"))

# Calculate total DNA quantity per row per hand
Bowman_data_with_total <- Bowman_data_renamed %>%
  pivot_wider(names_from = `Target Name`, values_from = Total_Quantity) %>%
  mutate(
    Cell_DNA = replace_na(Cell_DNA, 0),
    Cell_free_DNA = replace_na(Cell_free_DNA, 0),
    Total_DNA = Cell_DNA + Cell_free_DNA
  )

Bowman_data_with_total <- Bowman_data_with_total %>%
  mutate(Contact = recode(Contact,
                          "3s" = "A",
                          "15s" = "B"))

Bowman_data_with_total <- Bowman_data_with_total %>%
  mutate(Contact = factor(Contact, levels = c("A", "B")))

# Plot repeated data
plot_Bowman_repeat <- ggplot(Bowman_data_with_total, aes(x = Contact, y = Total_DNA)) +
  geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6, colour = "grey30") +
  scale_y_continuous(
    breaks = seq(0, 6, by = 2),
    limits = c(0, 6)) +
  labs(x = "\nScenario", y = "DNA Quantity (ng)") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Show plot
plot_Bowman_repeat

# export data for supplementary information
write.csv(Bowman_data_with_total, "./Results/Bowman_repeat_data.csv.csv", row.names = FALSE)

#### Final plot for article - Figure 6 ####
# Clean individual plots
plot_A <- plot_Bowman_data + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

plot_B <- plot_Bowman_repeat + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots with labels A and B
pCombined_Bowman_pending <- ggarrange(plot_A, plot_B,
                                    labels = c("A", "B"),
                                    ncol = 2, nrow = 1,
                                    align = "hv",
                                    common.legend = FALSE,
                                    font.label = list(size = 12, color = "black"),
                                    hjust = -0.5, vjust = 1.2) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))  # (Top, Right, Bottom, Left)

# Add shared y-axis label
pCombined_Bowman <- annotate_figure(
  pCombined_Bowman_pending,
  left = text_grob("DNA Quantity (ng)", rot = 90, vjust = 0.5, hjust = 0.5, size = 12),
  bottom = NULL  # You can add shared x-axis here if you want
)

# Show the combined plot
pCombined_Bowman

# Save the figure
ggsave("./Results/Bowman_combined_plot.png", pCombined_Bowman,
       width = 10, height = 5, dpi = 600, units = "in")

# ------------------------------------------------------------------------
# Section 7: Thomasma and Foran
# ------------------------------------------------------------------------
#### Data from the article ####
Thomasma_data <- read_csv("./Data/Approximate_Thomasma_dna_quantities.csv")

# Define elution volume in μL from Thomasma and Foran
elution_volume <- 20

# Convert concentration (pg/μL) to total DNA (ng)
Thomasma_data <- Thomasma_data %>%
  mutate(
    Total_DNA_ng = (Concentration * elution_volume) / 1000,
    Error_ng = (Error * elution_volume) / 1000
  )

# View result
print(Thomasma_data)

# Plot
plot_Thomasma_data <- ggplot(Thomasma_data, aes(x = Participant, y = Total_DNA_ng)) +
  geom_bar(stat = "identity", fill = "grey90", colour = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Total_DNA_ng, ymax = Total_DNA_ng + Error_ng),
                width = 0.2) +
  scale_y_continuous(
    limits = c(0, 6),  # adjust depending on your data
    breaks = seq(0, 6, by = 2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "\nParticipants",
    y = "DNA Quantity (ng)"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Show plot
plot_Thomasma_data

#### Repeated data with Proxy ####
# Select the columns of interest
Thomasma_repeat <- Thomasma_repeat%>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity", "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Thomasma_data_intermediate <- Thomasma_repeat %>%
  separate(`Sample Name`, into = c("Study", "Operator", "Finger", "Repeat"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study)) %>%
  select(Study, Operator, Finger,Repeat, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Thomasma_data_averaged <- Thomasma_data_intermediate %>%
  group_by(Study, Operator, Finger,Repeat, `Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Rename Target Name values
Thomasma_data_renamed <- Thomasma_data_averaged %>%
  mutate(`Target Name` = recode(`Target Name`,
                                "Mouse" = "Cell_DNA",
                                "TROUT 2" = "Cell_free_DNA"))

# Calculate total DNA quantity 
Thomasma_data_with_total <- Thomasma_data_renamed %>%
  pivot_wider(names_from = `Target Name`, values_from = Total_Quantity) %>%
  mutate(
    Cell_DNA = replace_na(Cell_DNA, 0),
    Cell_free_DNA = replace_na(Cell_free_DNA, 0),
    Total_DNA = Cell_DNA + Cell_free_DNA
  )

# Rename Target Name values
Thomasma_data_with_total <- Thomasma_data_with_total %>%
  mutate(Operator = recode(Operator,
                           "Op1" = "A",
                           "Op2" = "B",
                           "Op3" = "C"))

# Plot repeated data as a boxplot
plot_Thomasma_repeat <- Thomasma_data_with_total %>%
  ggplot(aes(x = Operator, y = Total_DNA)) +
  geom_boxplot(fill = "grey90", colour = "black", width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  labs(y = "DNA Quantity (ng)", x = "\nParticipants") +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, by = 2),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Show plot
plot_Thomasma_repeat

# export data for supplementary information
write.csv(Thomasma_data_with_total, "./Results/Thomasma_repeat_data.csv.csv", row.names = FALSE)

#### Final plot for article - Figure 7 ####
# Clean individual plots
plot_A <- plot_Thomasma_data + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

plot_B <- plot_Thomasma_repeat + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots with labels A and B
pCombined_Thomasma_pending <- ggarrange(plot_A, plot_B,
                                      labels = c("A", "B"),
                                      ncol = 2, nrow = 1,
                                      align = "hv",
                                      common.legend = FALSE,
                                      font.label = list(size = 12, color = "black"),
                                      hjust = -0.5, vjust = 1.2) +
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"))  # (Top, Right, Bottom, Left)

# Add shared y-axis label
pCombined_Thomasma <- annotate_figure(
  pCombined_Thomasma_pending,
  left = text_grob("DNA Quantity (ng)", rot = 90, vjust = 0.5, hjust = 0.5, size = 12),
  bottom = NULL  # You can add shared x-axis here if you want
)

# Show the combined plot
pCombined_Thomasma

# Save the figure
ggsave("./Results/Thomasma_combined_plot.png", pCombined_Thomasma,
       width = 10, height = 5, dpi = 600, units = "in")

# ------------------------------------------------------------------------
# Section 8: Lim et al.
# ------------------------------------------------------------------------
#### Data from the article ####
# Load the data
Lim_data <- read_csv("./Data/Lim_dna_quantities.csv")

# Define elution volume in μL
elution_volume <- 20

# Convert concentration (ng/μL) to total DNA (ng)
Lim_data <- Lim_data %>%
  mutate(Total_DNA_ng = Concentration * elution_volume)

# View result
print(Lim_data)

# Plot
plot_Lim_data <- ggplot(Lim_data, aes(x = "", y = Total_DNA_ng)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, fill = "white", colour = "black") +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6, colour = "grey30") +
  scale_y_continuous(
    breaks = seq(0, 4, by = 0.5),
    limits = c(0, 4)
  ) +
  labs(
    x = "Fingermarks",
    y = "DNA Quantity (ng)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_markdown(family = "Arial", size = 12),
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )

# Show plot
plot_Lim_data

#### Repeated data with Proxy ####
# Select the columns of interest
Lim_repeat <- Lim_repeat %>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Finger
Lim_data_intermediate <- Lim_repeat %>%
  separate(`Sample Name`, into = c("Study", "Operator","Finger", "Repeat"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study)) %>%
  select(Study, Operator, Finger, Repeat, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Lim_data_averaged <- Lim_data_intermediate %>%
  group_by(Study, Operator, Finger, Repeat, `Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Rename Target Name values
Lim_data_renamed <- Lim_data_averaged %>%
  mutate(`Target Name` = recode(`Target Name`,
                                "Mouse" = "Cell_DNA",
                                "TROUT 2" = "Cell_free_DNA"))

# Calculate total DNA quantity per row per Finger
Lim_data_with_total <- Lim_data_renamed %>%
  pivot_wider(names_from = `Target Name`, values_from = Total_Quantity) %>%
  mutate(
    Cell_DNA = replace_na(Cell_DNA, 0),
    Cell_free_DNA = replace_na(Cell_free_DNA, 0),
    Total_DNA = Cell_DNA + Cell_free_DNA
  )


# Plot repeated data
plot_Lim_repeat <- ggplot(Lim_data_with_total, aes(x = "", y = Total_DNA)) +
  geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6, colour = "grey30") +
  scale_y_continuous(
    breaks = seq(0, 6, by = 2),
    limits = c(0, 6)) +
  labs(x = "Fingermarks", y = "DNA Quantity (ng)") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Show plot
plot_Lim_repeat

# export data for supplementary information
write.csv(Lim_data_with_total, "./Results/Lim_repeat_data.csv.csv", row.names = FALSE)

#### Final plot for article - Figure 8 ####
# Clean individual plots
plot_A <- plot_Lim_data + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

plot_B <- plot_Lim_repeat + 
  rremove("ylab") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots with labels A and B
pCombined_Lim_pending <- ggarrange(plot_A, plot_B,
                                      labels = c("A", "B"),
                                      ncol = 2, nrow = 1,
                                      align = "hv",
                                      common.legend = FALSE,
                                      font.label = list(size = 12, color = "black"),
                                      hjust = -0.5, vjust = 1.2) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))  # (Top, Right, Bottom, Left)

# Add shared y-axis label
pCombined_Lim <- annotate_figure(
  pCombined_Lim_pending,
  left = text_grob("DNA Quantity (ng)", rot = 90, vjust = 0.5, hjust = 0.5, size = 12),
  bottom = NULL  # You can add shared x-axis here if you want
)

# Show the combined plot
pCombined_Lim

# Save the figure
ggsave("./Results/Lim_combined_plot.png", pCombined_Lim,
       width = 10, height = 5, dpi = 600, units = "in")