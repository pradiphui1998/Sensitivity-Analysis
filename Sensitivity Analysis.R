#install.packages("cutpointr")
library(cutpointr)
library(foreign)
#install.packages("pROC")
library(pROC)
Data <- read.spss("Data.sav",to.data.frame = TRUE)
Data=as.data.frame(Data)
attach(Data)
library(Publish)

Data$Metabolic_acidosis=ifelse(Data$Metabolic_acidosis=="Yes",1,0)
attach(Data)
#+-----------------------------------------------------------------------------+
#+ For MSI

cp<-cutpointr(Data,MSI,DEATHltgt48hrs, method = maximize_metric, metric=youden, 
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(DEATHltgt48hrs~MSI,data=Data)
auc(abc)
ci.auc(abc)

# Calculating ROC curve
abc <- roc(DEATHltgt48hrs ~ MSI, data = Data)

# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using MSI")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()

Data$MSI_group=ifelse(Data$MSI>=1.9,"Positive","Negative")
Data$MSI_group
attach(Data)

# Create a data frame from the main data set.
stu_data = data.frame(DEATHltgt48hrs,MSI_group)
# Create a contingency table with the needed variables.		
stu_data = table(DEATHltgt48hrs,MSI_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))



#+-----------------------------------------------------------------------------+
#+ For SI

cp<-cutpointr(Data,SI,DEATHltgt48hrs, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(DEATHltgt48hrs~SI,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using SI")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()

Data$SI_group=ifelse(Data$SI>=1.5,"Positive","Negative")
Data$SI_group
attach(Data)
sutable(DEATHltgt48hrs~SI_group,data=Data)

# Create a data frame from the main data set.
stu_data = data.frame(DEATHltgt48hrs,SI_group)
# Create a contingency table with the needed variables.		
stu_data = table(DEATHltgt48hrs,SI_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))

#+-----------------------------------------------------------------------------+
#+ For ASI

cp<-cutpointr(Data,ASI,DEATHltgt48hrs, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(DEATHltgt48hrs~ASI,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using ASI")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()




Data$ASI_group=ifelse(Data$ASI>=88.4,"Positive","Negative")
Data$ASI_group
attach(Data)
sutable(DEATHltgt48hrs~ASI_group,data=Data)

# Create a data frame from the main data set.
stu_data = data.frame(DEATHltgt48hrs,ASI_group)
# Create a contingency table with the needed variables.		
stu_data = table(DEATHltgt48hrs,ASI_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))

#+-----------------------------------------------------------------------------+
#+ For qSOFA

cp<-cutpointr(Data,qSOFA,DEATHltgt48hrs, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(DEATHltgt48hrs~qSOFA,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using qSOFA")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()



Data$qSOFA_group=ifelse(Data$qSOFA>=2,"Positive","Negative")
Data$qSOFA_group
attach(Data)
sutable(DEATHltgt48hrs~qSOFA_group,data=Data)

# Create a data frame from the main data set.
stu_data = data.frame(DEATHltgt48hrs,qSOFA_group)
# Create a contingency table with the needed variables.		
stu_data = table(DEATHltgt48hrs,qSOFA_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))




#+-----------------------------------------------------------------------------+
#+ For ANIONGAP

cp<-cutpointr(Data,ANIONGAP,DEATHltgt48hrs, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(DEATHltgt48hrs~ANIONGAP,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using ANIONGAP")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()



Data$ANIONGAP_group=ifelse(Data$ANIONGAP>=20.7,"Positive","Negative")
Data$ANIONGAP_group
attach(Data)
sutable(DEATHltgt48hrs~ANIONGAP_group,data=Data)


# Create a data frame from the main data set.
stu_data = data.frame(DEATHltgt48hrs,ANIONGAP_group)
# Create a contingency table with the needed variables.		
stu_data = table(DEATHltgt48hrs,ANIONGAP_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))


#+-----------------------------------------------------------------------------+
#+ Base deficit

cp<-cutpointr(Data,BASEDEFICIT,DEATHltgt48hrs, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(DEATHltgt48hrs~BASEDEFICIT,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using Base Deficit")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()



Data$BASEDEFICIT_group=ifelse(Data$BASEDEFICIT>=17,"Positive","Negative") 
attach(Data)
# Create a data frame from the main data set.
stu_data = data.frame(DEATHltgt48hrs,BASEDEFICIT_group)
# Create a contingency table with the needed variables.		
stu_data = table(DEATHltgt48hrs,BASEDEFICIT_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))





#+-----------------------------------------------------------------------------+
#+ Lactate

cp<-cutpointr(Data,LACTATE,DEATHltgt48hrs, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(DEATHltgt48hrs~LACTATE,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using LACTATE")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()



Data$LACTATE_group=ifelse(Data$LACTATE>=7.48,"Positive","Negative")
attach(Data)
sutable(DEATHltgt48hrs~LACTATE_group,data=Data)
# Create a data frame from the main data set.
stu_data = data.frame(DEATHltgt48hrs,LACTATE_group)
# Create a contingency table with the needed variables.		
stu_data = table(DEATHltgt48hrs,LACTATE_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))




#+-----------------------------------------------------------------------------+
#+ Metabolic Acidosis

cp<-cutpointr(Data,Metabolic_acidosis,DEATHltgt48hrs, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(DEATHltgt48hrs~Metabolic_acidosis,data=Data)
auc(abc)
ci.auc(abc)



# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "p", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using Metabolic Acidosis")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =",0.513),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()
















####################################################################################
#Morbidity
###################################################################################



#+-----------------------------------------------------------------------------+
#+ For MSI

cp<-cutpointr(Data,MSI,Morbidity, method = maximize_metric, metric=youden, 
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(Morbidity~MSI,data=Data)
auc(abc)
ci.auc(abc)

# Calculating ROC curve
abc <- roc(Morbidity ~ MSI, data = Data)

# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using MSI")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()

Data$MSI_group=ifelse(Data$MSI>=1.71,"Positive","Negative")
Data$MSI_group
attach(Data)
sutable(Morbidity~MSI_group,data=Data)

# Create a data frame from the main data set.
stu_data = data.frame(Morbidity,MSI_group)
# Create a contingency table with the needed variables.		
stu_data = table(Morbidity,MSI_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))






#+-----------------------------------------------------------------------------+
#+ For SI

cp<-cutpointr(Data,SI,Morbidity, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(Morbidity~SI,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using SI")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()

Data$SI_group=ifelse(Data$SI>=1.6,"Positive","Negative")
Data$SI_group
attach(Data)
sutable(Morbidity~SI_group,data=Data)


# Create a data frame from the main data set.
stu_data = data.frame(Morbidity,SI_group)
# Create a contingency table with the needed variables.		
stu_data = table(Morbidity,SI_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))


#+-----------------------------------------------------------------------------+
#+ For ASI

cp<-cutpointr(Data,ASI,Morbidity, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(Morbidity~ASI,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using ASI")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()




Data$ASI_group=ifelse(Data$ASI>=81.7,"Positive","Negative")
Data$ASI_group
attach(Data)
sutable(Morbidity~ASI_group,data=Data)

# Create a data frame from the main data set.
stu_data = data.frame(Morbidity,ASI_group)
# Create a contingency table with the needed variables.		
stu_data = table(Morbidity,ASI_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))

#+-----------------------------------------------------------------------------+
#+ For qSOFA

cp<-cutpointr(Data,qSOFA,Morbidity, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(Morbidity~qSOFA,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using qSOFA")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()



Data$qSOFA_group=ifelse(Data$qSOFA>=3,"Positive","Negative")
Data$qSOFA_group
attach(Data)
sutable(Morbidity~qSOFA_group,data=Data)


# Create a data frame from the main data set.
stu_data = data.frame(Morbidity,qSOFA_group)
# Create a contingency table with the needed variables.		
stu_data = table(Morbidity,qSOFA_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))




#+-----------------------------------------------------------------------------+
#+ For ANIONGAP

cp<-cutpointr(Data,ANIONGAP,Morbidity, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(Morbidity~ANIONGAP,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using ANIONGAP")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()



Data$ANIONGAP_group=ifelse(Data$ANIONGAP>=16.7,"Positive","Negative")
Data$ANIONGAP_group
attach(Data)
sutable(Morbidity~ANIONGAP_group,data=Data)

# Create a data frame from the main data set.
stu_data = data.frame(Morbidity,ANIONGAP_group)
# Create a contingency table with the needed variables.		
stu_data = table(Morbidity,ANIONGAP_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))




#+-----------------------------------------------------------------------------+
#+ Base deficit

cp<-cutpointr(Data,BASEDEFICIT,Morbidity, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(Morbidity~BASEDEFICIT,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using Base Deficit")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()



Data$BASEDEFICIT_group=ifelse(Data$BASEDEFICIT>=5.4,"Positive","Negative")
attach(Data)
sutable(Morbidity~BASEDEFICIT_group,data=Data)

# Create a data frame from the main data set.
stu_data = data.frame(Morbidity,BASEDEFICIT_group)
# Create a contingency table with the needed variables.		
stu_data = table(Morbidity,BASEDEFICIT_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))






#+-----------------------------------------------------------------------------+
#+ Lactate

cp<-cutpointr(Data,LACTATE,Morbidity, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(Morbidity~LACTATE,data=Data)
auc(abc)
ci.auc(abc)


# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using LACTATE")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =", round(auc(abc), 2)),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()



Data$LACTATE_group=ifelse(Data$LACTATE>=3.8,"Positive","Negative")
attach(Data)
sutable(Morbidity~LACTATE_group,data=Data)

# Create a data frame from the main data set.
stu_data = data.frame(Morbidity,LACTATE_group)
# Create a contingency table with the needed variables.		
stu_data = table(Morbidity,LACTATE_group)
print(stu_data)
# applying chisq.test() function
print(chisq.test(stu_data))







#+-----------------------------------------------------------------------------+
#+ Metabolic Acidosis
#
# Data$Metabolic_acidosis1=ifelse(Data$Metabolic_acidosis=="Yes",1,0)
# Data$Metabolic_acidosis1
# attach(Data)
cp<-cutpointr(Data,Metabolic_acidosis,Morbidity, method = maximize_metric, metric=youden,
              pos_class="Yes",na.rm=TRUE)
summary(cp)

abc<- roc(Morbidity~Metabolic_acidosis,data=Data)
auc(abc)
ci.auc(abc)



# Plotting ROC curve with reversed x-axis and modified labels
plot(1 - abc$specificities, abc$sensitivities,
     type = "l", lwd = 2, col = "blue",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ROC Curve using Metabolic Acidosis")

# Adding a diagonal reference line
abline(a = 0, b = 1, col = "black", lty = 2,lwd=2)

# Adding legend
legend("bottomright", legend = paste("AUC =",0.587),
       col = "black", lty = 1, cex = 1)

# Adding a grid for clarity
grid()
