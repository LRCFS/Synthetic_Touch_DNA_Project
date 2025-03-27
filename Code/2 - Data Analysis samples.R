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
Fonnelop_data <- data_by_study[["Fonneløp"]]
Goray_data <- data_by_study[["Goray"]]
Thomasma_data <- data_by_study[["Thomasma"]]
ThomasmaSDS_data <- data_by_study[["ThomasmaSDS"]]
Meakin_data <- data_by_study[["Meakin"]]
Daly_data <- data_by_study[["Daly"]]
Lim_data <- data_by_study[["Lim"]]

# ------------------------------------------------------------------------
# Section 1.1: Poetsch_data
# ------------------------------------------------------------------------
# Remove rows with unused data
Poetsch_data <- Poetsch_data %>%
  filter(`Sample Name` != "LH" & `Sample Name` != "LH_1" & `Sample Name` != "LH_2"
         & `Sample Name` != "RH" & `Sample Name` != "RH_1" & `Sample Name` != "RH_2"
         & `Sample Name` != "C")

# Select the columns of interest
Poetsch_data <- Poetsch_data %>%
  select(`Sample Name`, `Target Name`, "Quantity", "Quantity_Total","Study")

# Custom averaging function
custom_average <- function(values) {
  values <- na.omit(values)
  if (all(values == 0)) return(0)
  non_zero <- values[values != 0]
  mean(non_zero)
}

# First, split 'Sample Name' into Participant, Repeat, and Hand
data_intermediate <- Poetsch_data %>%
  separate(`Sample Name`, into = c("Study", "Repeat", "Hand"), sep = "_", remove = FALSE) %>%
  mutate(Study = paste0(Study, Repeat)) %>%
  select(Study, Hand, `Target Name`, Quantity_Total)

#Apply averaging explicitly outside of main pipeline
data_averaged <- data_intermediate %>%
  group_by(Study, Hand, `Target Name`) %>%
  summarise(Total_Quantity = custom_average(Quantity_Total), .groups = 'drop')

# Pivot wider to obtain final format
Poetsch_data <- data_averaged %>%
  pivot_wider(names_from = c(`Target Name`, Hand),
              values_from = Total_Quantity,
              names_sep = "_") %>%
  rename(
    "Cell_free_L" = `TROUT 2_L`,
    "Cell_free_R" = `TROUT 2_R`,
    "Cell_DNA_L" = Mouse_L,
    "Cell_DNA_R" = Mouse_R
  ) %>%
  arrange(Study)

# Calculate total DNA (Cell DNA + Cell-free DNA)
Poetsch_data_total <- Poetsch_data %>%
  mutate(
    Total_L = Cell_free_L + Cell_DNA_L,
    Total_R = Cell_free_R + Cell_DNA_R
  )

# Step 2: Compute descriptive statistics clearly
descriptive_stats <- Poetsch_data_total %>%
  pivot_longer(cols = c(Total_L, Total_R),
               names_to = "Hand",
               values_to = "Total_DNA") %>%
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
print(Poetsch_data_total)
print(descriptive_stats)

# Add a column for clearer labels
descriptive_stats <- descriptive_stats %>%
  mutate(Study = ifelse(Hand == "Total_L", "Our study - Left", "Our study - Right")) %>%
  mutate(Mean = (Min + Max) / 2)

# Create the other study data
other_study <- tibble(
  Study = "Other study",
  Min = 0,
  Max = 4.8,
)

# Combine
combined_data <- bind_rows(
  descriptive_stats %>% select(Study, Min, Max, Mean),
  other_study
)

# Plot
ggplot(combined_data, aes(x = Study, ymin = Min, ymax = Max)) +
  geom_errorbar(width = 0.3, size = 1.2, color = "black") +
  geom_point(aes(y = Mean), size = 3, color = "red") +
  labs(
    title = "Comparison of DNA Yield per Handprint",
    subtitle = "Our study (L/R) vs. published study (general range)",
    y = "DNA yield (ng)",
    x = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# #### Trout ####
# # Calculate the average of Quantity_Total for each duplicate in Sample Name
# data_hand_trout <- data_hand_trout %>%
#   group_by(`Sample Name`, `Target Name`) %>% # Group the data by Sample Name and Target Name
#   filter(n() == 2) %>%   # Ensure each group has exactly 2 rows
#   summarize( # Calculate the average based on the following conditions:
#     average_quantity_total = if (all(Quantity_Total == 0)) {
#       0  # If both values are 0, return 0 as the average
#     } else {
#       mean(Quantity_Total[Quantity_Total > 0], na.rm = TRUE)  # Calculate the average of non-zero values
#     },
#     .groups = "drop"
#   )
# 
# # Check the number of rows = 36
# nrow(data_hand_trout)
# 
# # split sample name based on "_"
# result_data_hand_trout_Extended <- data.frame(str_split(data_hand_trout$`Sample Name`, "_", simplify=TRUE))
# names(result_data_hand_trout_Extended) <- c("Sebum","Individuals","Repeat","Hand")
# result_data_hand_trout <- cbind(result_data_hand_trout_Extended,`Target Name`=data_hand_trout$`Target Name`, average_quantity_total=data_hand_trout$average_quantity_total)
# rm(result_data_hand_trout_Extended)
# 
# # Group and summarize data by Individuals and Hand
# result_data_hand_trout <- result_data_hand_trout %>%
#   group_by(Individuals, Hand, Sebum) %>%
#   summarise(average_quantity = mean(average_quantity_total, na.rm = TRUE))
# 
# # to separate the different concentrations
# result_data_hand_trout_200 <- result_data_hand_trout %>%
#   filter(grepl("200",Sebum ))
# result_data_hand_trout_300 <- result_data_hand_trout %>%
#   filter(grepl("300",Sebum ))
# result_data_hand_trout_400 <- result_data_hand_trout %>%
#   filter(grepl("400",Sebum ))
# 
# # Create the bar plot
# Plot_result_data_hand_trout_200 <- ggplot(result_data_hand_trout_200, aes(x = as.factor(Individuals), y = average_quantity, fill = Hand)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "Operators (Individuals)", y = "Average quantity of DNA (ng)", fill = "Hand") +
#   scale_fill_manual(values = c("gray70", "gray30")) +  # Set custom gray shades
#   ylim(0,2.5)+
#   theme_bw()
# Plot_result_data_hand_trout_200
# 
# Plot_result_data_hand_trout_300 <- ggplot(result_data_hand_trout_300, aes(x = as.factor(Individuals), y = average_quantity, fill = Hand)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "Operators (Individuals)", y = "Average quantity of DNA (ng)", fill = "Hand") +
#   scale_fill_manual(values = c("gray70", "gray30")) +  # Set custom gray shades
#   ylim(0,2.5)+
#   theme_bw()
# Plot_result_data_hand_trout_300
# 
# Plot_result_data_hand_trout_400 <- ggplot(result_data_hand_trout_400, aes(x = as.factor(Individuals), y = average_quantity, fill = Hand)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "Operators (Individuals)", y = "Average quantity of DNA (ng)", fill = "Hand") +
#   scale_fill_manual(values = c("gray70", "gray30")) +  # Set custom gray shades
#   ylim(0,2.5)+
#   theme_bw()
# Plot_result_data_hand_trout_400
# 
# #### Final graph  #### 
# pCombined_pending <- ggarrange(Plot_result_data_hand_trout_200+ rremove("ylab") + rremove("xlab"),
#                                Plot_result_data_hand_trout_300+ rremove("ylab") + rremove("xlab"),
#                                Plot_result_data_hand_trout_400+ rremove("ylab") + rremove("xlab"),
#                                labels = c("A   ","B   ","C   "),
#                                common.legend = T,legend = "right",
#                                ncol = 1, nrow = 3,
#                                font.label = list(size = 12, color = "black", family = NULL, position = "top"),
#                                hjust=0.5,vjust=1.5)+
#   theme(plot.margin = margin(0.5,1,0.5,1, "cm")) # in order (Top,right,bottom,left)
# pCombined_pending
# 
# pCombined_trout <- annotate_figure(pCombined_pending, left = textGrob("Average quantity of DNA (ng)", rot = 90, vjust = 0.5, hjust = 0.5, gp = gpar(cex =1)),
#                                bottom = textGrob("Operators (Individuals)", vjust = 0.5, hjust = 0.5,gp = gpar(cex = 1)))
# pCombined_trout
# 
# # to save the graph
# ggsave("Figure 7 - Shedding with manual pressure.png", pCombinedHP, width =8, height = 7, units = "in", dpi=600,path = "Results")
# #### Mouse ####
# # Calculate the average of Quantity_Total for each duplicate in Sample Name
# data_hand_mouse <- data_hand_mouse %>%
#   group_by(`Sample Name`, `Target Name`) %>% # Group the data by Sample Name and Target Name
#   filter(n() == 2) %>%   # Ensure each group has exactly 2 rows
#   summarize( # Calculate the average based on the following conditions:
#     average_quantity_total = if (all(Quantity_Total == 0)) {
#       0  # If both values are 0, return 0 as the average
#     } else {
#       mean(Quantity_Total[Quantity_Total > 0], na.rm = TRUE)  # Calculate the average of non-zero values
#     },
#     .groups = "drop"
#   )
# 
# # Check the number of rows = 36
# nrow(data_hand_mouse)
# 
# # split sample name based on "_"
# result_data_hand_mouse_Extended <- data.frame(str_split(data_hand_mouse$`Sample Name`, "_", simplify=TRUE))
# names(result_data_hand_mouse_Extended) <- c("Sebum","Individuals","Repeat","Hand")
# result_data_hand_mouse <- cbind(result_data_hand_mouse_Extended,`Target Name`=data_hand_mouse$`Target Name`, average_quantity_total=data_hand_mouse$average_quantity_total)
# rm(result_data_hand_mouse_Extended)
# 
# # Group and summarize data by Individuals and Hand
# result_data_hand_mouse <- result_data_hand_mouse %>%
#   group_by(Individuals, Hand) %>%
#   summarise(average_quantity = mean(average_quantity_total, na.rm = TRUE))
# 
# # Create the bar plot
# ggplot(result_data_hand_mouse, aes(x = as.factor(Individuals), y = average_quantity, fill = Hand)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "Operators (Individuals)", y = "Average quantity of DNA (ng)", fill = "Hand") +
#   scale_fill_manual(values = c("gray70", "gray30")) +  # Set custom gray shades
#   ylim(0,10)+
#   theme_bw()
