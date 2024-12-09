##########################################################################
#####                           Data upload                          #####
##########################################################################

# This R script is the first step to export all the data needed to generate the figure in the article.
# This includes the data from: 
# 1: Transfer from hands

# ------------------------------------------------------------------------
# Section 1: Transfer from hands
# ------------------------------------------------------------------------
# Specify the folder path
folder_path <- "Data/Hand transfer/"

# Get the list of .xls files in the folder
file_list <- list.files(path = folder_path, pattern = "*.xls", full.names = TRUE)

# Initialise an empty list to store dataframes
df_list <- list()

# Loop through each file, read the file, and remove the first 6 rows (header)
for (file in file_list) {
  # Read the file and skip the first 6 rows
  data <- read_xls(file, skip = 6)
  
  # Append the dataframe to the list
  df_list[[length(df_list) + 1]] <- data
}

# Combine all dataframes in the list into one dataframe using bind_rows
merged_data <- bind_rows(df_list)

# View the merged data (optional)
print(merged_data)

# Sort the data by the 'Sample Name' column
sorted_merged_data <- merged_data %>% arrange(`Sample Name`)

# Select the column of interest
sorted_data <- sorted_merged_data %>%
  select("Well",`Sample Name`,`Target Name`,"Quantity", "Task")

# Calculating the total amount of DNA in each sample
sorted_data$Quantity_Total <- sorted_data$Quantity * 50

# Create a new dataframe for rows with "STANDARD" in the Task column
Standard_hand <- sorted_data %>%
  filter(Task == "STANDARD")

# Remove rows with "STANDARD" in the Task column and where "Sample Name" contains "NTC"
data_filtered <- sorted_data %>%
  filter(Task != "STANDARD" & Task != "NTC")

# Replace all NA values with 0 in the data_filtered dataframe
data_filtered <- data_filtered %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

# Creating the different dataframes for analysis
data_hand_trout <- data_filtered %>%
  filter(grepl("^trout", `Target Name`))
data_hand_mouse <- data_filtered %>%
  filter(grepl("^mouse", `Target Name`))

# Removing the control from the above datasets
Controls_hand_trout <- data_hand_trout %>%
  filter(`Sample Name` %in% c("200_C", "300_C", "400_C"))
data_hand_trout <- data_hand_trout %>%
  filter(`Sample Name` != "200_C" & `Sample Name` != "300_C" & `Sample Name` != "400_C")

Controls_hand_mouse <- data_hand_mouse %>%
  filter(`Sample Name` %in% c("200_C", "300_C", "400_C"))
data_hand_mouse <- data_hand_mouse %>%
  filter(`Sample Name` != "200_C" & `Sample Name` != "300_C" & `Sample Name` != "400_C")

# Remove unnecessary dataframe
rm(data,data_filtered,df_list,merged_data, sorted_data,sorted_merged_data)
