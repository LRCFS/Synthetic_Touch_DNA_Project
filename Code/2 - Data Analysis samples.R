##########################################################################
#####                      Data analysis washes                      #####
##########################################################################

# This R script is the first step to export all the data needed to generate the figure in the article.
# This includes the data from: 
# XXX

# ------------------------------------------------------------------------
# Section 1: XXX
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Section 1.1: XXX
# ------------------------------------------------------------------------
#### Trout ####
# Calculate the average of Quantity_Total for each duplicate in Sample Name
data_hand_trout <- data_hand_trout %>%
  group_by(`Sample Name`, `Target Name`) %>% # Group the data by Sample Name and Target Name
  filter(n() == 2) %>%   # Ensure each group has exactly 2 rows
  summarize( # Calculate the average based on the following conditions:
    average_quantity_total = if (all(Quantity_Total == 0)) {
      0  # If both values are 0, return 0 as the average
    } else {
      mean(Quantity_Total[Quantity_Total > 0], na.rm = TRUE)  # Calculate the average of non-zero values
    },
    .groups = "drop"
  )

# Check the number of rows = 36
nrow(data_hand_trout)

# split sample name based on "_"
result_data_hand_trout_Extended <- data.frame(str_split(data_hand_trout$`Sample Name`, "_", simplify=TRUE))
names(result_data_hand_trout_Extended) <- c("Sebum","Individuals","Repeat","Hand")
result_data_hand_trout <- cbind(result_data_hand_trout_Extended,`Target Name`=data_hand_trout$`Target Name`, average_quantity_total=data_hand_trout$average_quantity_total)
rm(result_data_hand_trout_Extended)

# Group and summarize data by Individuals and Hand
result_data_hand_trout <- result_data_hand_trout %>%
  group_by(Individuals, Hand) %>%
  summarise(average_quantity = mean(average_quantity_total, na.rm = TRUE))

# Create the bar plot
ggplot(result_data_hand_trout, aes(x = as.factor(Individuals), y = average_quantity, fill = Hand)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Operators (Individuals)", y = "Average quantity of DNA (ng)", fill = "Hand") +
  scale_fill_manual(values = c("gray70", "gray30")) +  # Set custom gray shades
  ylim(0,10)+
  theme_bw()

#### Mouse ####
# Calculate the average of Quantity_Total for each duplicate in Sample Name
data_hand_mouse <- data_hand_mouse %>%
  group_by(`Sample Name`, `Target Name`) %>% # Group the data by Sample Name and Target Name
  filter(n() == 2) %>%   # Ensure each group has exactly 2 rows
  summarize( # Calculate the average based on the following conditions:
    average_quantity_total = if (all(Quantity_Total == 0)) {
      0  # If both values are 0, return 0 as the average
    } else {
      mean(Quantity_Total[Quantity_Total > 0], na.rm = TRUE)  # Calculate the average of non-zero values
    },
    .groups = "drop"
  )

# Check the number of rows = 36
nrow(data_hand_mouse)

# split sample name based on "_"
result_data_hand_mouse_Extended <- data.frame(str_split(data_hand_mouse$`Sample Name`, "_", simplify=TRUE))
names(result_data_hand_mouse_Extended) <- c("Sebum","Individuals","Repeat","Hand")
result_data_hand_mouse <- cbind(result_data_hand_mouse_Extended,`Target Name`=data_hand_mouse$`Target Name`, average_quantity_total=data_hand_mouse$average_quantity_total)
rm(result_data_hand_mouse_Extended)

# Group and summarize data by Individuals and Hand
result_data_hand_mouse <- result_data_hand_mouse %>%
  group_by(Individuals, Hand) %>%
  summarise(average_quantity = mean(average_quantity_total, na.rm = TRUE))

# Create the bar plot
ggplot(result_data_hand_mouse, aes(x = as.factor(Individuals), y = average_quantity, fill = Hand)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Operators (Individuals)", y = "Average quantity of DNA (ng)", fill = "Hand") +
  scale_fill_manual(values = c("gray70", "gray30")) +  # Set custom gray shades
  ylim(0,10)+
  theme_bw()
