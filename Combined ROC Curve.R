# Load required libraries
library(cutpointr)
library(foreign)
library(pROC)
library(Publish)

# Load your data and preprocess it as you've done before

# Create ROC curves for each variable
roc_MSI <- roc(DEATHltgt48hrs ~ MSI, data = Data)
roc_MSI
roc_SI <- roc(DEATHltgt48hrs ~ SI, data = Data)
roc_SI
roc_ASI <- roc(DEATHltgt48hrs ~ ASI, data = Data)
roc_ASI
roc_qSOFA <- roc(DEATHltgt48hrs ~ qSOFA, data = Data)
roc_qSOFA

# Calculate 95% confidence intervals for AUC
ci_MSI <- ci.auc(roc_MSI)
ci_MSI
ci_SI <- ci.auc(roc_SI)
ci_SI
ci_ASI <- ci.auc(roc_ASI)
ci_ASI
ci_qSOFA <- ci.auc(roc_qSOFA)
ci_qSOFA

# Create a new plot
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "ROC Curves for Mortality")

# Add ROC curves to the plot
lines(1 - roc_MSI$specificities, roc_MSI$sensitivities, col = "blue", lwd = 2, label = "MSI")
lines(1 - roc_SI$specificities, roc_SI$sensitivities, col = "red", lwd = 2, label = "SI")
lines(1 - roc_ASI$specificities, roc_ASI$sensitivities, col = "green", lwd = 2, label = "ASI")
lines(1 - roc_qSOFA$specificities, roc_qSOFA$sensitivities, col = "purple", lwd = 2, label = "qSOFA")

# Add legend
legend("bottomright", legend = c("MSI", "SI", "ASI", "qSOFA"),
       col = c("blue", "red", "green", "purple"), lwd = 2)

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)

# Adding a grid for clarity
grid()




##############################################################################
##############################################################################

# Create ROC curves for each variable
roc_BASEDEFICIT <- roc(DEATHltgt48hrs ~ BASEDEFICIT, data = Data)
roc_BASEDEFICIT
roc_ANIONGAP <- roc(DEATHltgt48hrs ~ ANIONGAP, data = Data)
roc_ANIONGAP
roc_Metabolic_acidosis <- roc(DEATHltgt48hrs ~ Metabolic_acidosis, data = Data)
roc_Metabolic_acidosis
roc_LACTATE <- roc(DEATHltgt48hrs ~LACTATE, data = Data)
roc_LACTATE


# Calculate 95% confidence intervals for AUC
roc_BASEDEFICIT <- ci.auc(roc_BASEDEFICIT)
roc_BASEDEFICIT
roc_ANIONGAP <- ci.auc(roc_ANIONGAP)
roc_ANIONGAP
roc_Metabolic_acidosis <- ci.auc(roc_Metabolic_acidosis)
roc_Metabolic_acidosis
roc_LACTATE <- ci.auc(roc_LACTATE)
roc_LACTATE



# Create a new plot
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "ROC Curves for Mortality")

# Add ROC curves to the plot
lines(1 - roc_BASEDEFICIT$specificities, roc_BASEDEFICIT$sensitivities, col = "blue", lwd = 2, label = "MSI")
lines(1 - roc_ANIONGAP$specificities, roc_ANIONGAP$sensitivities, col = "red", lwd = 2, label = "SI")
lines(1 - roc_Metabolic_acidosis$specificities, roc_Metabolic_acidosis$sensitivities, col = "green", lwd = 2, label = "ASI")
lines(1 - roc_LACTATE$specificities, roc_LACTATE$sensitivities, col = "purple", lwd = 2, label = "qSOFA")

# Add legend
legend("bottomright", legend = c("BASEDEFICIT", "ANIONGAP", "Metabolic_acidosis", "LACTATE"),
       col = c("blue", "red", "green", "purple"), lwd = 2)

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)

# Adding a grid for clarity
grid()








###############################################################################
## Morbidity
################################################################################

# Create ROC curves for each variable
roc_MSI <- roc(Morbidity ~ MSI, data = Data)
roc_MSI
roc_SI <- roc(Morbidity ~ SI, data = Data)
roc_SI
roc_ASI <- roc(Morbidity ~ ASI, data = Data)
roc_ASI
roc_qSOFA <- roc(Morbidity ~ qSOFA, data = Data)
roc_qSOFA


# Calculate 95% confidence intervals for AUC
ci_MSI <- ci.auc(roc_MSI)
ci_MSI
ci_SI <- ci.auc(roc_SI)
ci_SI
ci_ASI <- ci.auc(roc_ASI)
ci_ASI
ci_qSOFA <- ci.auc(roc_qSOFA)
ci_qSOFA


# Create a new plot
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "ROC Curves for Morbidity")

# Add ROC curves to the plot
lines(1 - roc_MSI$specificities, roc_MSI$sensitivities, col = "blue", lwd = 2, label = "MSI")
lines(1 - roc_SI$specificities, roc_SI$sensitivities, col = "red", lwd = 2, label = "SI")
lines(1 - roc_ASI$specificities, roc_ASI$sensitivities, col = "green", lwd = 2, label = "ASI")
lines(1 - roc_qSOFA$specificities, roc_qSOFA$sensitivities, col = "purple", lwd = 2, label = "qSOFA")

# Add legend
legend("bottomright", legend = c("MSI", "SI", "ASI", "qSOFA"),
       col = c("blue", "red", "green", "purple"), lwd = 2)

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)

# Adding a grid for clarity
grid()




##############################################################################
##############################################################################
# Create ROC curves for each variable
roc_BASEDEFICIT <- roc(Morbidity ~ BASEDEFICIT, data = Data)
roc_BASEDEFICIT
roc_ANIONGAP <- roc(Morbidity ~ ANIONGAP, data = Data)
roc_ANIONGAP
roc_Metabolic_acidosis <- roc(Morbidity ~ Metabolic_acidosis, data = Data)
roc_Metabolic_acidosis
roc_LACTATE <- roc(Morbidity ~LACTATE, data = Data)
roc_LACTATE

# Calculate 95% confidence intervals for AUC
roc_BASEDEFICIT <- ci.auc(roc_BASEDEFICIT)
roc_BASEDEFICIT
roc_ANIONGAP <- ci.auc(roc_ANIONGAP)
roc_ANIONGAP
roc_Metabolic_acidosis <- ci.auc(roc_Metabolic_acidosis)
roc_Metabolic_acidosis
roc_LACTATE <- ci.auc(roc_LACTATE)
roc_LACTATE


# Create a new plot
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "ROC Curves for Morbidity")

# Add ROC curves to the plot
lines(1 - roc_BASEDEFICIT$specificities, roc_BASEDEFICIT$sensitivities, col = "blue", lwd = 2, label = "MSI")
lines(1 - roc_ANIONGAP$specificities, roc_ANIONGAP$sensitivities, col = "red", lwd = 2, label = "SI")
lines(1 - roc_Metabolic_acidosis$specificities, roc_Metabolic_acidosis$sensitivities, col = "green", lwd = 2, label = "ASI")
lines(1 - roc_LACTATE$specificities, roc_LACTATE$sensitivities, col = "purple", lwd = 2, label = "qSOFA")

# Add legend
legend("bottomright", legend = c("BASEDEFICIT", "ANIONGAP", "Metabolic_acidosis", "LACTATE"),
       col = c("blue", "red", "green", "purple"), lwd = 2)

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)

# Adding a grid for clarity
grid()
