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
# Section 2: Poetsch_data
# ------------------------------------------------------------------------
# Remove rows with unused data
Poetsch_data <- Poetsch_data %>%
  filter(`Sample Name` != "LH" & `Sample Name` != "LH_1" & `Sample Name` != "LH_2"
         & `Sample Name` != "RH" & `Sample Name` != "RH_1" & `Sample Name` != "RH_2"
         & `Sample Name` != "C")

# Select the columns of interest
Poetsch_data <- Poetsch_data %>%
  select(`Sample Name`, `Target Name`, "Quantity", "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Poetsch_data_intermediate <- Poetsch_data %>%
  separate(`Sample Name`, into = c("Study", "Repeat", "Hand"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study, Repeat)) %>%
  select(Study, Hand, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Poetsch_data_averaged <- Poetsch_data_intermediate %>%
  group_by(Study, Hand, `Target Name`) %>%
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
# Section 2: Fonnelop_data
# ------------------------------------------------------------------------
# Remove rows with unused data
Fonnelop_data <- Fonnelop_data %>%
  filter(`Sample Name` != "LH_3" & `Sample Name` != "RH_3" & `Sample Name` != "C_2" & `Sample Name` != "C")

# Select the columns of interest
Fonnelop_data <- Fonnelop_data %>%
  select(`Sample Name`, `Target Name`, "Quantity", "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Fonnelop_data_intermediate <- Fonnelop_data %>%
  separate(`Sample Name`, into = c("Study", "Repeat", "Hand"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study, Repeat)) %>%
  select(Study, Hand, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Fonnelop_data_averaged <- Fonnelop_data_intermediate %>%
  group_by(Study, Hand, `Target Name`) %>%
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

ggplot(combined_data_clean, aes(x = Study, y = Total_DNA)) +
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
# Section 3: Goray_data
# ------------------------------------------------------------------------
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
# Section 2: Thomasma_data
# ------------------------------------------------------------------------
# Remove rows with unused data
Fonnelop_data <- Fonnelop_data %>%
  filter(`Sample Name` != "LH_3" & `Sample Name` != "RH_3" & `Sample Name` != "C_2" & `Sample Name` != "C")

# Select the columns of interest
Fonnelop_data <- Fonnelop_data %>%
  select(`Sample Name`, `Target Name`, "Quantity", "Quantity_Total","Study")

# First, split 'Sample Name' into Participant, Repeat, and Hand
Fonnelop_data_intermediate <- Fonnelop_data %>%
  separate(`Sample Name`, into = c("Study", "Repeat", "Hand"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study, Repeat)) %>%
  select(Study, Hand, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
Fonnelop_data_averaged <- Fonnelop_data_intermediate %>%
  group_by(Study, Hand, `Target Name`) %>%
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

ggplot(combined_data_clean, aes(x = Study, y = Total_DNA)) +
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
