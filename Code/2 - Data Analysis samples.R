##########################################################################
#####                      Data analysis washes                      #####
##########################################################################

# This R script is the first step to export all the data needed to generate the figure in the article.
# This includes the data from: 
# XXX

# ------------------------------------------------------------------------
# Section 1: XXX
# ------------------------------------------------------------------------
# Remove rows where 'Cт' column has "Undetermined"
Data <- Data %>%
  filter(Cт != "Undetermined")

# Create a list of dataframes split by Study name
data_by_study <- split(Data, Data$Study)

# Extract them individually as named dataframes
Poetsch_repeat <- data_by_study[["Poetsch"]]
Fonnelop_repeat <- data_by_study[["Fonnelop"]]
Goray_repeat <- data_by_study[["Goray"]]
Thomasma_repeat <- data_by_study[["Thomasma"]]
Thomasma_repeat$Quantity_Total_pgul <- ((Thomasma_repeat$Quantity_Total)/5)*1000 # 5ul for the PCR and 1000 to convert ng to pg
Meakin_repeat <- data_by_study[["Meakin"]]
Daly_repeat <- data_by_study[["Daly"]]
Lim_repeat <- data_by_study[["Lim"]]

# ------------------------------------------------------------------------
# Section 2: Poetsch et al.
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

# Plot repeated data
plot_Poetsch_repeat <- Poetsch_data_with_total %>%
  ggplot(aes(x = Concentration, y = Total_DNA)) +
  stat_summary(fun = mean, geom = "bar", fill = "#C6DBEF", colour = "black", width = 0.7) +
  stat_summary(fun.data = function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    data.frame(y = mean_val, ymin = mean_val, ymax = mean_val + sd_val)
    },
    geom = "errorbar", width = 0.2) +
  labs(y = "DNA Quantity (ng)", x = "Input concentration") +
  scale_y_continuous(
    limits = c(0, 5),
    breaks = seq(0, 5, by = 0.5),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw()+
  theme(
    text = element_text(family = "Arial", size = 14),  # Set overall font
    axis.title = element_text(size = 14),              # Axis title font size
    axis.text = element_text(size = 12),               # Axis tick labels
  )

plot_Poetsch_repeat

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
# Section 5: Meakin et al.
# ------------------------------------------------------------------------
#### Data from the article ####
Meakin_data <- read_csv("./Data/Approximate_Meakin_dna_quantities.csv")

# Plot
plot_Meakin_data <- ggplot(Meakin_data, aes(x = Knife, y = Total_DNA)) +
  geom_bar(stat = "identity", fill = "grey90", colour = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Total_DNA, ymax = Total_DNA + SD),
                width = 0.2) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Knife",
    y = "DNA Quantity (ng)"
  ) +
  theme_bw()

# Show plot
plot_Meakin_data

#### Repeated data with Proxy ####
# Select the columns of interest
Meakin_data <- Meakin_data %>%
  select(`Sample Name`, `Target Name`, "Quantity", "Quantity_Total","Study")

#Apply averaging explicitly outside of main pipeline
Meakin_data_averaged <- Meakin_data %>%
  group_by(`Sample Name`,`Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Pivot wider to obtain final format
Meakin_data <- Meakin_data_averaged %>%
  pivot_wider(names_from = c(`Target Name`),values_from = Total_Quantity,names_sep = "_") %>%
  rename("Cell_free" = `TROUT 2`,"Cell_DNA" = Mouse) %>%
  arrange(`Sample Name`)

# Calculate total DNA (Cell DNA + Cell-free DNA)
# treating NA values as 0
Meakin_data_total <- Meakin_data %>%
  mutate(Total = rowSums(cbind(Cell_free, Cell_DNA), na.rm = TRUE))

# Prepare your study data
Meakin_data_our_study <- Meakin_data_total %>%
  select(`Sample Name`, Total) %>%
  rename(Knife = `Sample Name`, Total_DNA = Total) %>%
  mutate(
    Knife = gsub("_", " ", Knife),
    SD = 0,
    Source = "Proxy study"
  )

# Combine real values from Meakin with your study values
Meakin_study_data <- tibble(
  Knife = c("Knife 1", "Knife 2", "Knife 3", "Knife 4"),
  Total_DNA = c(3.4, 0.9, 1.2, 10.4),   # Mean DNA concentration in ng
  SD = c(0.5, 0.8, 0.5, 3.7), # Standard deviation in ng
  Source = "Meakin et al."
)

# Combine both datasets
combined_data <- bind_rows(Meakin_study_data, Meakin_data_our_study)

# Plot
ggplot(combined_data, aes(x = Knife, y = Total_DNA, fill = Source)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Total_DNA - SD, ymax = Total_DNA + SD),
                width = 0.2) +
  facet_wrap(~Source, nrow = 1) +  # One row → side-by-side
  labs(
    title = "DNA Recovery per Knife",
    x = "Knife",
    y = "Total DNA (ng)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

# ------------------------------------------------------------------------
# Section 6: Daly et al.
# ------------------------------------------------------------------------
#### Data from the article ####
# Create the data frame directly
Daly_data <- tibble(
  Material = "Fabric",
  Mean_DNA_ng = 1.23,
  Error_ng = 14.8 - 1.23  # Using the range: max - mean
)

# Plot
plot_Daly_data <- ggplot(Daly_data, aes(x = Material, y = Mean_DNA_ng)) +
  geom_bar(stat = "identity", fill = "grey80", colour = "black", width = 0.6) +
  geom_errorbar(aes(ymin = 0, ymax = Mean_DNA_ng + Error_ng), width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Daly et al.",
    y = "DNA Quantity (ng)"
  ) +
  theme_bw()

# Show plot
plot_Daly_data

#### Repeated data with Proxy ####

# ------------------------------------------------------------------------
# Section 7: Bowman et al.
# ------------------------------------------------------------------------
#### Data from the article ####
Bowman_data <- read_csv("./Data/Bowman_dna_quantities.csv")

# Plot
plot_Bowman_data <- ggplot(Bowman_data, aes(x = Scenario, y = DNA_Quantity_ng)) +
  geom_boxplot(fill = c("white", "grey"), colour = "black", width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, shape = 21, fill = "black") +
  labs(
    x = "Scenario",
    y = "DNA Quantity (ng)"
  ) +
  theme_bw()

# Show plot
plot_Bowman_data

# ------------------------------------------------------------------------
# Section 4: Thomasma and Foran
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
  geom_bar(stat = "identity", fill = "grey70", colour = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Total_DNA_ng, ymax = Total_DNA_ng + Error_ng),
                width = 0.2) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Participant",
    y = "DNA Quantity (ng)"
  ) +
  theme_bw

# Show plot
plot_Thomasma_data

#### Repeated data with Proxy ####
# Remove rows with unused data
Thomasma_data <- Thomasma_data %>%
  filter(`Sample Name` != "99_L" & `Sample Name` != "C10_L" & `Sample Name` != "C_L" & `Sample Name` != "C_M" & `Sample Name` != "C_T")

# load correction list
CorrectionSample <- read.csv("Data/Studies/Thomasma_sampleNameCorrected.txt", sep = "\t", header = TRUE)
CorrectionSample <- as.data.frame(CorrectionSample)
Thomasma_data$`Sample Name Corrected` <- gsr(as.character(Thomasma_data$`Sample Name`),as.character(CorrectionSample$Sample),as.character(CorrectionSample$Sample_corrected))

# Select the columns of interest
Thomasma_data <- Thomasma_data %>%
  select(`Sample Name Corrected`, `Target Name`, "Quantity", "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Thomasma_data_intermediate <- Thomasma_data %>%
  separate(`Sample Name Corrected`, into = c("Operator", "Finger", "Solution"), sep = "_", remove = FALSE) %>%
  select(Operator, Solution, Finger,`Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Thomasma_data_averaged <- Thomasma_data_intermediate %>%
  group_by(Operator, Solution, Finger,`Target Name`) %>%
  summarise(Avg_Quantity = custom_average(Quantity_Total), .groups = "drop")

# Pivot wider to obtain final format
Thomasma_data <- Thomasma_data_averaged %>%
  pivot_wider(names_from = c(`Target Name`),values_from = Avg_Quantity,names_sep = "_") %>%
  rename("Cell_free" = `TROUT 2`,"Cell_DNA" = Mouse)

# Calculate total DNA (Cell DNA + Cell-free DNA)
# treating NA values as 0
Thomasma_data_total <- Thomasma_data %>%
  mutate(Total = rowSums(cbind(Cell_free, Cell_DNA), na.rm = TRUE))

# Compute descriptive statistics clearly
descriptive_stats <- Thomasma_data_total %>%
  group_by(Operator, Solution) %>%
  summarise(
    Mean = mean(Total, na.rm = TRUE),
    SD = sd(Total, na.rm = TRUE),
    Median = median(Total, na.rm = TRUE),
    Min = min(Total, na.rm = TRUE),
    Max = max(Total, na.rm = TRUE),
    .groups = 'drop'
  )

# View results
print(descriptive_stats)

# ------------------------------------------------------------------------
# Section 4: Lim et al.
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
  geom_boxplot(fill = "grey80", colour = "black", width = 0.3, outlier.shape = NA) +
  geom_jitter(aes(x = ""), width = 0.1, size = 2, shape = 21, fill = "black") +
  scale_y_continuous(
    breaks = seq(0, 4, by = 0.5),
    limits = c(0, 4),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "",
    y = "DNA Quantity (ng)"
  ) +
  theme_bw()

# Show plot
plot_Lim_data
