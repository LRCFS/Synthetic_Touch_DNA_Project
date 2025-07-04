##########################################################################
#####                           Data upload                          #####
##########################################################################

# This R script is the first step to export all the data needed to generate the figure in the article.
# This includes the data from: 
# 1: Transfer from hands

# ------------------------------------------------------------------------
# Section 1: Data from studpy replication
# ------------------------------------------------------------------------
# Specify the folder path
folder_path <- "Data/Studies"

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

# Replace names of Samples
CorrectionSample <- read.csv("CorrectionLists/SampleNameCorrected.txt", sep = "\t", header = TRUE)
CorrectionSample <- as.data.frame(CorrectionSample)
merged_data$`Sample Name corrected` <- gsr(as.character(merged_data$`Sample Name`),as.character(CorrectionSample$Name),as.character(CorrectionSample$NameCorrected))

# Sort the data by the 'Sample Name' column
sorted_filtered_data <- merged_data %>% arrange(`Sample Name corrected`)

# Select the columns of interest
sorted_data <- sorted_filtered_data %>%
  select(`Sample Name`,`Sample Name corrected`,"Well", "Cт", `Target Name`, "Quantity", "Task", "Study")

# Create a new dataframe for rows with "STANDARD" in the Task column
Standard <- sorted_data %>%
  filter(Task == "STANDARD")

# Remove rows with "STANDARD" or "NTC" in the Task column
Data <- sorted_data %>%
  filter(Task != "STANDARD" & Task != "NTC")

# Updated calculation for Quantity_Total based on Study volumes
factor_lookup <- c("Poetsch" = 50, "Thomasma" = 50, "Fonnelop" = 50, "Goray" = 50, "Meakin" = 10, "Daly" = 10, "Lim" = 50, "Bowman" = 50)

Data <- Data %>%
  mutate(Factor = factor_lookup[Study]) %>%
  mutate(Quantity_Total = Quantity * Factor)

rm(data, df_list, merged_data, sorted_filtered_data, sorted_data)


# ------------------------------------------------------------------------
# Section 2: Create dataframes for analysis
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
Bowman_repeat <- data_by_study[["Bowman"]]
