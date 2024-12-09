##########################################################################
#####                     Data analysis controls                     #####
##########################################################################

# This R script is the first step to export all the data needed to generate the figure in the article.
# This includes the data from: 
# XXX

# ------------------------------------------------------------------------
# Section 1: Direct to buffer
# ------------------------------------------------------------------------
# Calculate the average Quantity_Total per Target Name
avg_quantity <- Control_dtb %>%
  group_by(`Target Name`) %>%
  summarize(avg_quantity = mean(Quantity_Total, na.rm = TRUE))

avg_quantity$avg_quantity_pg <- avg_quantity$avg_quantity *1000
