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
Poetsch_data <- data_by_study[["Poetsch"]]
Fonnelop_data <- data_by_study[["Fonnelop"]]
Goray_data <- data_by_study[["Goray"]]
Thomasma_data <- data_by_study[["Thomasma"]]
Thomasma_data$Quantity_Total_pgul <- ((Thomasma_data$Quantity_Total)/5)*1000 # 5ul for the PCR and 1000 to convert ng to pg
Meakin_data <- data_by_study[["Meakin"]]
Daly_data <- data_by_study[["Daly"]]
Lim_data <- data_by_study[["Lim"]]

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
  geom_bar(stat = "identity", fill = "#6BAED6", colour = "black", width = 0.7) +
  geom_errorbar(aes(ymin = DNA_qt, ymax = DNA_qt + Error), width = 0.3) +
  scale_y_continuous(
    breaks = seq(0, 6, by = 0.5),   # Set ticks from 1 to 6 every 0.5
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  labs(x = "Age Group",y = "DNA Quantity (ng)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

# Show plot
plot_Poetsch_data

#### Repeated data with Proxy ####
# Remove rows with unused data
Poetsch_data <- Poetsch_data %>%
  filter(`Sample Name` != "LH" & `Sample Name` != "LH_1" & `Sample Name` != "LH_2"
         & `Sample Name` != "RH" & `Sample Name` != "RH_1" & `Sample Name` != "RH_2"
         & `Sample Name` != "C")

# Select the columns of interest
Poetsch_data <- Poetsch_data %>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Poetsch_data_intermediate <- Poetsch_data %>%
  separate(`Sample Name`, into = c("Study", "Operator","Concentration", "Hand"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study)) %>%
  select(Study, Operator, Concentration, Hand, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Poetsch_data_averaged <- Poetsch_data_intermediate %>%
  group_by(Study, Operator, Concentration, Hand, `Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Pivot wider to obtain final format
Poetsch_data <- Poetsch_data_averaged %>%
  pivot_wider(names_from = c(`Target Name`, Hand),values_from = Total_Quantity,names_sep = "_") %>%
  rename("Cell_free_L" = `TROUT 2_L`,"Cell_free_R" = `TROUT 2_R`,"Cell_DNA_L" = Mouse_L,"Cell_DNA_R" = Mouse_R) %>%
  arrange(Study)

# Calculate total DNA (Cell DNA + Cell-free DNA)
# treating NA values as 0
Poetsch_data_total <- Poetsch_data %>%
  mutate(Total_L = rowSums(cbind(Cell_free_L, Cell_DNA_L), na.rm = TRUE),
         Total_R = rowSums(cbind(Cell_free_R, Cell_DNA_R), na.rm = TRUE))

# Compute descriptive statistics clearly
descriptive_stats <- Poetsch_data_total %>%
  pivot_longer(cols = c(Total_L, Total_R),names_to = "Hand",values_to = "Total_DNA") %>%
  group_by(Hand) %>%
  summarise(
    Mean = mean(Total_DNA, na.rm = TRUE),
    SD = sd(Total_DNA, na.rm = TRUE),
    Median = median(Total_DNA, na.rm = TRUE),
    Min = min(Total_DNA, na.rm = TRUE),
    Max = max(Total_DNA, na.rm = TRUE),
    .groups = 'drop'
  )
# View results
print(descriptive_stats)

# Reshape the data into long format (raw values)
Poetsch_data_our_study <- Poetsch_data_total %>%
  pivot_longer(cols = c(Total_L, Total_R),
               names_to = "Hand", values_to = "Total_DNA") %>%
  mutate(Study = case_when(Hand == "Total_L" ~ "Proxy DNA – Left hand",Hand == "Total_R" ~ "Proxy DNA – Right hand"))

# Create a dummy entry for the other study
Poetsch_study_data <- tibble(
  Study = "Poetsch et al.",
  Total_DNA = NA_real_,
  Min = 0,
  Max = 4.8
)

ggplot() +
  # Boxplots for your study
  geom_boxplot(data = Poetsch_data_our_study,
               aes(x = Study, y = Total_DNA, fill = Study),
               width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(data = Poetsch_data_our_study,
              aes(x = Study, y = Total_DNA),
              width = 0.15, height = 0, size = 2, color = "black", alpha = 0.6) +
  
  # Bar with errorbar for other study
  geom_bar(data = Poetsch_study_data,
           aes(x = Study, y = (Min + Max)/2),
           stat = "identity", fill = "NA", width = 0.5, alpha = 0.8) +
  geom_errorbar(data = Poetsch_study_data,
                aes(x = Study, ymin = Min, ymax = Max),
                width = 0.2, size = 1) +
  
  labs(title = "Comparison of DNA transfer",
       subtitle = "Poetsch et al. study vs. Repeat",
       y = "DNA quantity (ng)", x = "") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "none")

rm(Poetsch_data_averaged, Poetsch_data_intermediate,Poetsch_data_total,Poetsch_data)

# ------------------------------------------------------------------------
# Section 2: Fonnelop et al.
# ------------------------------------------------------------------------
#### Data from the article ####
fonnelop_data <- read_csv("./Data/Approximate_fonnelop_dna_quantities.csv") %>%
  mutate(Study = "Fonnelop et al.") %>%
  rename(Total_DNA = Quantity)

# Plot
plot_fonnelop_data <- ggplot(fonnelop_data, aes(x = Study, y = Total_DNA, colour = Gender)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  geom_boxplot(width = 0.3, outlier.shape = NA, colour = "black", fill = NA) +
  scale_colour_manual(values = c("F" = "grey40", "M" = "steelblue")) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 0.5),
    limits = c(0, 4)) +
  labs(x = NULL,y = "DNA Quantity (ng)") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(face = "italic"))

# Show plot
plot_fonnelop_data

#### Repeated data with Proxy ####
# Remove rows with unused data
Fonnelop_data <- Fonnelop_data %>%
  filter(`Sample Name` != "LH_3" & `Sample Name` != "RH_3" & `Sample Name` != "C_2" & `Sample Name` != "C")

# Select the columns of interest
Fonnelop_data <- Fonnelop_data %>%
  select(`Sample Name`=`Sample Name corrected`, `Target Name`, "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Fonnelop_data_intermediate <- Fonnelop_data %>%
  separate(`Sample Name`, into = c("Study", "Operator","Concentration", "Hand"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study)) %>%
  select(Study, Operator, Concentration, Hand, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Fonnelop_data_averaged <- Fonnelop_data_intermediate %>%
  group_by(Study, Operator, Concentration, Hand, `Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Pivot wider to obtain final format
Fonnelop_data <- Fonnelop_data_averaged %>%
  pivot_wider(names_from = c(`Target Name`, Hand),values_from = Total_Quantity,names_sep = "_") %>%
  rename("Cell_free_L" = `TROUT 2_L`,"Cell_free_R" = `TROUT 2_R`,"Cell_DNA_L" = Mouse_L,"Cell_DNA_R" = Mouse_R) %>%
  arrange(Study)

# Calculate total DNA (Cell DNA + Cell-free DNA)
# treating NA values as 0
Fonnelop_data_total <- Fonnelop_data %>%
  mutate(Total_L = rowSums(cbind(Cell_free_L, Cell_DNA_L), na.rm = TRUE),
         Total_R = rowSums(cbind(Cell_free_R, Cell_DNA_R), na.rm = TRUE))

# Compute descriptive statistics clearly
descriptive_stats <- Fonnelop_data_total %>%
  pivot_longer(cols = c(Total_L, Total_R),names_to = "Hand",values_to = "Total_DNA") %>%
  group_by(Hand) %>%
  summarise(
    Mean = mean(Total_DNA, na.rm = TRUE),
    SD = sd(Total_DNA, na.rm = TRUE),
    Median = median(Total_DNA, na.rm = TRUE),
    Min = min(Total_DNA, na.rm = TRUE),
    Max = max(Total_DNA, na.rm = TRUE),
    .groups = 'drop'
  )
# View results
print(descriptive_stats)

# Reshape the data into long format (raw values)
Fonnelop_data_our_study <- Fonnelop_data_total %>%
  pivot_longer(cols = c(Total_L, Total_R),
               names_to = "Hand", values_to = "Total_DNA") %>%
  mutate(Study = case_when(Hand == "Total_L" ~ "Proxy DNA – Left hand",Hand == "Total_R" ~ "Proxy DNA – Right hand"))

# Combine real values from Fonnelop with your study values
fonnelop_study_data <- read_csv("./Data/approximate_fonnelop_dna_quantities.csv") %>%
  mutate(Study = "Fonnelop et al.") %>%
  rename(Total_DNA = Quantity)

# Combine the two datasets
combined_data <- bind_rows(
  fonnelop_study_data %>% select(Study, Total_DNA, Gender, Participant),
  Fonnelop_data_our_study %>% select(Study, Total_DNA) %>% mutate(Gender = NA, Participant = NA)
)

# Calculate the average
# Mean value of 0.64 ng reported in Fonnelop et al. (actual value); 
# plotted points represent approximated extraction from figure
combined_data %>%
  group_by(Study) %>%
  summarise(
    Mean_DNA = mean(Total_DNA, na.rm = TRUE),
    SD_DNA = sd(Total_DNA, na.rm = TRUE),
    n = sum(!is.na(Total_DNA))
  )

ggplot(combined_data, aes(x = Study, y = Total_DNA)) +
  geom_boxplot(aes(fill = Study), outlier.shape = NA, alpha = 0.7, width = 0.5, show.legend = FALSE) +
  geom_jitter(aes(color = Gender), width = 0.15, height = 0, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("F" = "salmon", "M" = "turquoise"), na.value = "black") +
  labs(
    title = "Comparison of DNA transfer",
    subtitle = "Fonnelop et al. study vs. proxy DNA repeat",
    y = "DNA quantity (ng)", x = "",
    color = "Gender"  # ensures legend is labeled clearly
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "right"
  )

rm(Fonnelop_data_averaged, Fonnelop_data_intermediate,Fonnelop_data_total,Fonnelop_data)

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
  labs(x = "Hand",y = "DNA Quantity (ng)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),axis.title = element_text(size = 11))

# Show plot
plot_Goray_data

#### Repeated data with Proxy ####
#### Visualisation ####
# Select the columns of interest
Goray_data <- Goray_data %>%
  select(`Sample Name`, `Target Name`, "Quantity", "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Goray_data_intermediate <- Goray_data %>%
  separate(`Sample Name`, into = c("Study", "Repeat", "Hand"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study, Repeat)) %>%
  select(Study, Hand, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Goray_data_averaged <- Goray_data_intermediate %>%
  group_by(Study, Hand, `Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Pivot wider to obtain final format
Goray_data <- Goray_data_averaged %>%
  pivot_wider(names_from = c(`Target Name`, Hand),values_from = Total_Quantity,names_sep = "_") %>%
  rename("Cell_free_L" = `TROUT 2_L`,"Cell_free_R" = `TROUT 2_R`,"Cell_DNA_L" = Mouse_L,"Cell_DNA_R" = Mouse_R) %>%
  arrange(Study)

# Calculate total DNA (Cell DNA + Cell-free DNA)
# treating NA values as 0
Goray_data_total <- Goray_data %>%
  mutate(Total_L = rowSums(cbind(Cell_free_L, Cell_DNA_L), na.rm = TRUE),
         Total_R = rowSums(cbind(Cell_free_R, Cell_DNA_R), na.rm = TRUE))

# Reshape the data into long format (raw values)
Goray_data_our_study <- Goray_data_total %>%
  pivot_longer(cols = c(Total_L, Total_R),
               names_to = "Hand", values_to = "Total_DNA") %>%
  mutate(Study = case_when(Hand == "Total_L" ~ "Proxy DNA",Hand == "Total_R" ~ "Proxy DNA"))

# Compute descriptive statistics clearly
descriptive_stats <- Goray_data_total %>%
  pivot_longer(cols = c(Total_L, Total_R),names_to = "Hand",values_to = "Total_DNA") %>%
  group_by(Hand) %>%
  summarise(
    Mean = mean(Total_DNA, na.rm = TRUE),
    SD = sd(Total_DNA, na.rm = TRUE),
    Median = median(Total_DNA, na.rm = TRUE),
    Min = min(Total_DNA, na.rm = TRUE),
    Max = max(Total_DNA, na.rm = TRUE),
    .groups = 'drop'
  )
# View results
print(descriptive_stats)

# Combine real values from Goray with your study values
Goray_study_data <- read_csv("./Data/Goray_tidy_data.csv") %>%
  mutate(Study = "Goray et al.") %>%
  rename(Total_DNA = Quantity)

# Remove one point for visualisation purpose
Goray_study_data <- Goray_study_data %>%
  filter(!(Replicate == "D" & Hand == "R" & Day == "3" & Time == 1 & Total_DNA == 38.2))

ggplot(Goray_study_data, aes(x = Hand, y = Total_DNA, fill = Hand)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, aes(color = Hand)) +
  scale_fill_manual(values = c("L" = "salmon", "R" = "turquoise")) +
  scale_color_manual(values = c("L" = "salmon", "R" = "turquoise")) +
  labs(
    title = "DNA quantities from Goray et al.",
    subtitle = "Grouped by hand (Left/Right) across all days",
    x = "Hand",
    y = "DNA quantity (ng)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Replace names
Goray_data_our_study <- Goray_data_our_study %>%
  mutate(Hand = recode(Hand, "Total_L" = "L", "Total_R" = "R"))

# Combine the two datasets
combined_data <- bind_rows(
  Goray_study_data %>% select(Study, Hand, Total_DNA),
  Goray_data_our_study %>% select(Study, Hand,Total_DNA))

ggplot(combined_data, aes(x = Study, y = Total_DNA, fill = Hand)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = Hand),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("L" = "salmon", "R" = "turquoise")) +
  scale_color_manual(values = c("L" = "salmon", "R" = "turquoise")) +
  labs(
    title = "Comparison of DNA",
    subtitle = "Goray et al. vs Proxy DNA",
    x = "Study",
    y = "DNA quantity (ng)",
    fill = "Hand",
    color = "Hand"
  ) +
  theme_minimal(base_size = 14)

rm(Goray_data_averaged, Goray_data_intermediate,Goray_data_total,Goray_data)

#### Statistics ####
#  Wilcoxon rank-sum test: Left Hand Only
left_hand <- combined_data %>%
  filter(Hand == "L")

wilcox.test(Total_DNA ~ Study, data = combined_data %>% filter(Hand == "L"))

# Wilcoxon rank-sum test: right Hand Only
right_hand <- combined_data %>%
  filter(Hand == "R")

wilcox.test(Total_DNA ~ Study, data = combined_data %>% filter(Hand == "R"))

# Wilcoxon rank-sum test: Both Hands Combined
wilcox.test(Total_DNA ~ Study, data = combined_data)

# ------------------------------------------------------------------------
# Section 5: Meakin et al.
# ------------------------------------------------------------------------
#### Data from the article ####
Meakin_data <- read_csv("./Data/Approximate_Meakin_dna_quantities.csv")

# Plot
plot_Meakin_data <- ggplot(Meakin_data, aes(x = Knife, y = Total_DNA)) +
  geom_bar(stat = "identity", fill = "grey70", colour = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Total_DNA, ymax = Total_DNA + SD),
                width = 0.2) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Knife",
    y = "DNA Quantity (ng)"
  ) +
  theme_minimal()

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
