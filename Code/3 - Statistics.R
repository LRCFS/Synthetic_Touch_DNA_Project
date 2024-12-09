##########################################################################
#####                      Data analysis washes                      #####
##########################################################################

# This R script is the first step to export all the data needed to generate the figure in the article.
# This includes the data from: 
# XXX

# ------------------------------------------------------------------------
# Section 1: Comparison between donor garments
# ------------------------------------------------------------------------
## Difference between operators
# Gentle pressure
result_W_W1_donor %>%
  shapiro_test(Area.mm2) # p-value = 0.196, normal distribution

Handpressure_O4 %>%
  shapiro_test(Area.mm2) # p-value = 0.000109, non-parametric distribution

Handpressure_O5_gentle <- filter(Handpressure_O5, Pressure == "gentle")
Handpressure_O5_gentle %>%
  shapiro_test(Area.mm2) # p-value = 0.191, normal distribution

forstats_gentle_operators <- rbind(Handpressure_O1,Handpressure_O4,Handpressure_O5_gentle)
forstats_gentle_operators <- forstats_gentle_operators %>%
  select("Operator","Area.mm2", "Garment")

### sphericity test
#Levene’s Test
#H0: All sample variances are equal
#H1: At least one group has a variance that is not equal to the rest.
forstats_gentle_operators$Operator <- as.factor(forstats_gentle_operators$Operator)
leveneTest(Area.mm2 ~ Operator, forstats_gentle_operators) # p-value = 0.08089 , equal variances