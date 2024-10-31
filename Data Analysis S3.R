ecodata <- read.csv("C:\\Users\\sos523\\Desktop\\Data S2.csv")


### Preprocessing ###
str(ecodata)
ecodata$Region <- as.factor(ecodata$Region)

any(is.na(data))
# no missing data except for historical pathogen index 

### visualizations ###

## histograms 

variables <- c("GDP", "No_Epidemics", "No_Disasters", "Death_Epidemics", "Death_Disasters", "Mortality_Epidemics", "Mortality_Disasters", "GCI")
par(mfrow = c(3, 3))
for (var in variables) {
  hist(ecodata[[var]], main = paste("Histogram of", var), xlab = var)
}
par(mfrow = c(1, 1))
# all predictors are right-skewed, heavy-tailed distributions 


## boxplots
variables <- c("No_Epidemics", "No_Disasters", "Death_Epidemics", "Death_Disasters", "Mortality_Epidemics", "Mortality_Disasters", "GDP")
par(mfrow = c(3, 3))
for (var in variables) {
  boxplot(ecodata[[var]], main = paste("Boxplot of", var), ylab = "Value")
}
par(mfrow = c(1, 1))
# all predictors include extreme values 


### Correlations ###
library(Hmisc)

# Pearson 
rcorr(as.matrix(ecodata[c("GCI", "No_Epidemics", "No_Disasters", "Death_Epidemics", "Death_Disasters", "Mortality_Epidemics", "Mortality_Disasters", "GDP")]))
# as data is right-skewed and contains potential outliers --> Spearman is better

#Spearman 
rcorr(as.matrix(ecodata[c("GCI", "No_Epidemics", "No_Disasters","Death_Epidemics", "Death_Disasters", "Mortality_Epidemics", "Mortality_Disasters", "GDP")]),
      type = "spearman")


### mixed-effect modeling


## transforming predictors 
# predictors are in different units and ranges; 
# to make predictors comparable and for model convergence with lmer transformation is required  



## Rank Based Inverse Normal Transformation 
library(RNOmni)

# variables to transform  
variables <- c("GDP", "No_Epidemics", "No_Disasters", "Death_Epidemics", "Death_Disasters", "Mortality_Epidemics", "Mortality_Disasters")

for(var in variables) {
  ecodata[[var]] <- RankNorm(ecodata[[var]], ties.method = "first")
}

# verify transformation
hist(ecodata$Mortality_Epidemics)



### modeling number of epidemics and disasters ###  
library(lmerTest)
model_number <- lmer(GCI ~ GDP + No_Epidemics + No_Disasters + (1|Region), data = ecodata, REML = TRUE)
summary(model_number)
confint(model_number)
## model diagnostics

# checking residuals
plot(model_number)
 
hist(residuals(model_number))
# residuals normally distributed 

# multicollinearity 
library(car)
vif(model_number) 
# VIFs below 1.69,  so no multicollinearity issues 


# outliers
library(influence.ME)

cooksd <- cooks.distance(model_number)
plot(cooksd, type="h", main="Cook's Distance")

cutoff <- 4/(nrow(ecodata) - length(fixef(model_number)) - 1)
outliers <- which(cooksd > cutoff)

new_data <- ecodata[-outliers,]
clean_model <- lmer(GCI ~ GDP + No_Epidemics + No_Disasters + (1|Region), data = new_data, REML = TRUE)
summary(clean_model)
# removing 26 influential observations based on cook's distance didn't change 
# the significance of predictors


# model comparison
model_number2 <- lmer(GCI ~ GDP + No_Epidemics + No_Disasters + (1|Region), data = ecodata, REML = FALSE)
model_number2_reduced <- lmer(GCI ~ GDP + No_Disasters + (1|Region), data = ecodata, REML = FALSE)

anova(model_number2_reduced, model_number2, test = F)
# adding No_Epidemics didn't significantly improved the model


#########################################################

### modeling deaths caused by epidemics and disasters 

model_death <- lmer(GCI ~ GDP + Death_Epidemics + Death_Disasters + (1|Region), data = ecodata, REML = TRUE)
summary(model_death)

plot(model_death)

vif(model_death) 

cooksd <- cooks.distance(model_death)
plot(cooksd, type="h", main="Cook's Distance")

cutoff <- 4/(nrow(ecodata) - length(fixef(model_death)) - 1)
outliers <- which(cooksd > cutoff)

new_data <- ecodata[-outliers,]
clean_model <- lmer(GCI ~ GDP + Death_Epidemics + Death_Disasters + (1|Region), data = new_data, REML = TRUE)
summary(clean_model)
# removing outliers didn't change the significance of predictors 

confint(model_death)

#########################################################
## Mortality of epidemics and disasters
model_mortality <- lmer(GCI ~ GDP + Mortality_Epidemics + Mortality_Disasters + (1|Region), data = ecodata, REML = TRUE)
summary(model_mortality)

# model diagnostics
plot(model_mortality)
hist(residuals(model_mortality))
#residuals are normally distributed 

confint(model_mortality)

vif(model_mortality) 
# all VIFs below 1.45

# outliers
cooksd <- cooks.distance(model_mortality)
plot(cooksd, type="h", main="Cook's Distance")

cutoff <- 4/(nrow(ecodata) - length(fixef(model_mortality)) - 1)
outliers <- which(cooksd > cutoff)

new_data <- ecodata[-outliers,]
clean_model <- lmer(GCI ~ GDP + Mortality_Epidemics + Mortality_Disasters + (1|Region), data = new_data, REML = TRUE)
summary(clean_model)
# removing 27 influential observations didn't change significance of predictors



# extra analysis
model_quadratic <- lmer(GCI ~ GDP + No_Epidemics + No_Disasters + I(GDP^2) + I(No_Epidemics^2) + I(No_Disasters^2) + (1|Region), data = ecodata, REML = FALSE)  
summary(model_quadratic)

anova(model_number, model_quadratic)
# no sig. difference between linear and quadratic so 
# the relationship is not non-linear 



model_interaction <- lmer(GCI ~ GDP + No_Epidemics + No_Disasters + GDP:No_Epidemics + GDP:No_Disasters + No_Epidemics:No_Disasters + (1|Region), data = ecodata, REML = FALSE)
summary(model_interaction)

anova(model_number, model_interaction)
# no interaction either


######
# Historic Pathogen Index 
ecodata$Log_GDP = log(ecodata$GDP)
ecodata$HPI_7 <- as.numeric(ecodata$HPI_7)
ecodata <- ecodata[!is.na(ecodata$HPI_7), ]


model_HPI <- lmer(GCI ~ Log_GDP + HPI_7 + (1|Region), data = ecodata, REML = TRUE)
summary(model_HPI)
# HPI_7 was not significant 


## the end ## 