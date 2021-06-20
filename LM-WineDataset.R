rm(list=ls())
# Libraries
if(!require(corrplot)){
  install.packages("corrplot")
  library(corrplot)
}
if(!require(caret)){
  install.packages("caret")
  library(caret)
}
if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}
if(!require(MASS)){
  install.packages("MASS")
  library(MASS)
}
if(!require(nnet)){
  install.packages("nnet")
  library(nnet)
}
if(!require(car)){
  install.packages("car")
  library(car)
}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(knitr)){
  install.packages("knitr")
  library(knitr)
}
if(!require(lmtest)){
  install.packages("lmtest")
  library(lmtest)
}
if(!require(stargazer)){
  install.packages("stargazer")
  library(stargazer)
}
if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)
}
if(!require(broom)){
  install.packages("broom")
  library(broom)
}
if(!require(gtools)){
  install.packages("gtools")
  library(gtools)
}
set.seed(112358)

# Import data
data <- read.csv('winequality-red.csv', header=T)

head(data, n=4)
View(data)
# Descriptive Statistics 

# Summary Statistics of the data
# We remove the summary statistics of the response because it is ordinal
summary(data)

correlation <- cor(data)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlation, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         diag=FALSE 
)

# Scatter plots
featurePlot(x = data[, 1:11], 
            y = as.numeric(data$quality))

# Histograms of features
ggplot(gather(data), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

# Statistical Modeling 

# Full Linear Model
lm1 <- lm(quality ~ ., data = data)
summary(lm1)
summary(lm1)$r.squared
summary(lm1)$adj.r.squared
summary(lm1)$sigma

# Variable Selection 

## Backward step selection
step(lm1, direction="backward")

# Reduced OLS
lm2 <- lm(formula = quality ~ volatile.acidity + chlorides + free.sulfur.dioxide + 
            total.sulfur.dioxide + pH + sulphates + alcohol, data = data)
summary(lm2)  
summary(lm2)
summary(lm2)$adj.r.squared
summary(lm2)$sigma

### VIF ###
vif(lm1)
vif(lm2)

### Model Comparison ###

anova(lm2, lm1, test="F")

### Constant Variance Check ###

lm2.df <- augment(lm2)

# Increasing trend and clear pattern which makes sense because the response is ordinal
ggplot(lm2.df, aes(x=lm2.df$.fitted, y=lm2.df$.std.resid))+ xlab('Fitted')+ylab("Studentized Residuals")+
  geom_point(shape=1) + geom_hline(yintercept = 0, color="red")

### Component + Residual Plot ###

crPlots(lm2, main="")

### QQ Plot ###
qqPlot(lm2$residuals, "norm", ylab="Residuals")


## Outiers 

distances = data.frame(distance=cooks.distance(lm2))
ggplot(distances, aes(x=1:nrow(distances), y=distances$distance)) + geom_point(shape=1) +
  labs(x="Index", y="Cooks Distance")

# number of regressors in the model
k <- 7
n <- nrow(data)
cutoff.cook <- k/(n-k-1)
inflpoints <- cooks.distance(lm2)[cooks.distance(lm2) > cutoff.cook] 
length(inflpoints)

indexinfpoints <- as.numeric(names(inflpoints))

# Remove outliers
data.noOutliers <- data[-indexinfpoints,]

# Refit the model without outliers
lm3 <- lm(formula = quality ~ volatile.acidity + chlorides + free.sulfur.dioxide + 
            total.sulfur.dioxide + pH + sulphates + alcohol, data = data.noOutliers)

summary(lm3)  

### QQ plot for model without outliers ###
qqPlot(lm3$residuals, "norm", ylab="Residuals")

# Transform for logistic regression and Ordinal logistic
my_data <- data
my_data$quality <- factor(my_data$quality, ordered = T)
my_data$binqual <- recode(my_data$quality, "c(3,4,5) = 0")
my_data$binqual <- recode(my_data$binqual, "c(6,7,8) = 1")
my_data$binqual <- factor(my_data$binqual)
bin_data <- my_data[,-12]
my_data <- my_data[,-13]
my_data1 <- my_data[which(my_data$quality == 5),]
my_data2 <- my_data[which(my_data$quality == 6),]
my_data3 <- my_data[which(my_data$quality == 7),]
my_data <- rbind(my_data1,my_data2,my_data3)
my_data$quality <- ordered(my_data$quality)

# Training and test sets for both models

tr <- sample(1:nrow(bin_data),0.7*nrow(bin_data))
tr2 <- sample(1:nrow(my_data),0.7*nrow(my_data))
train1 <- my_data[tr2,]
test1 <- my_data[-tr2,]
train2 <- bin_data[tr,]
test2 <- bin_data[-tr,]
summary(bin_data)
summary(my_data)

#Logistic Reg (Using binary transformation data set)
logMod1 <- glm(binqual~., data = train2, family = "binomial")

#Backward elimination 
b_back <- stepAIC(logMod1) 
b_best <- formula(b_back)

#Keep model indicated by elimination
logMod2 <- glm(formula = b_best, data = train2, family = "binomial")
summary(logMod2)

#Omit free.sulphur.dioxide as indicated by summary:
logMod3 <- glm(binqual~fixed.acidity + volatile.acidity + citric.acid + residual.sugar + chlorides + 
                 total.sulfur.dioxide + density + sulphates + alcohol, 
         data = train2, family = "binomial")
summary(logMod3)

# As we can see AIC values contradict the summary's p-value for free.sulphur.dioxide 
#so use LRT
lrtest(logMod2,logMod3)

# p-value from LRT and Wald test from summary a bit different
#both indicate to omit free.sulphur.dioxide 

#investigate with C.I:
confint(logMod2)

# Calculate AIC/BIC values 
AIC1 <- summary(logMod2)$aic
AIC2 <- summary(logMod3)$aic
BIC1 <- BIC(logMod2)
BIC2 <- BIC(logMod3)
aicbic <- cbind(rbind(AIC1,AIC2),rbind(BIC1,BIC2))
rownames(aicbic) <- c("Model with free.sulfur.dioxide","Model without free.sulfur.dioxide")
colnames(aicbic) <- c("AIC","BIC")
aicbic

# Continue with most parsimonious model
final_tab <- cbind(summary(logMod3)$coefficient, confint(logMod3) )
final_tab

#Calculate Accuracy of final model
fitted.results <- predict(logMod3,newdata=test2, type = "response")
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test2$binqual)
print(paste("Accuracy is",1-misClasificError))

#Accuracy of training set to get a hint for consistency of the model
fitted.results1 <- predict(logMod3,newdata=train2, type = "response")
fitted.results1 <- ifelse(fitted.results1 > 0.5,1,0)
misClasificError1 <- mean(fitted.results1 != train2$binqual)
print(paste("Accuracy is",1-misClasificError1))

##Ordinal 
#PROPORTIONAL ODDS ASSUMPTION CHECK
aa1 <- logLik(multinom(quality ~ ., data = my_data))
aa2 <- logLik(polr(quality~ ., data = my_data))
G <- -2*(aa2[1] - aa1[1])
pchisq(G,22-12, lower.tail = F)

# Full Ordinal model
mod1 <- polr(quality ~ . , data = train1, Hess = T)

#Backward elimination
bbest <- step(mod1)
use1 <- formula(bbest)

#Keep model indicated by previous procedure
mod2 <- polr(use1, data = train1, Hess = T)
summary(mod2)

#Compute p-values and visualize alongside with the estimates

coefs <- coef(summary(mod2))
p_val <- pnorm(abs(coefs[,"t value"]), lower.tail = F)*2
coefs <- cbind(coefs,"p-value" = p_val)
coefs

#Omit chlorides as highest p-value and refit
mod3 <- polr(quality ~ fixed.acidity + volatile.acidity + citric.acid + residual.sugar + total.sulfur.dioxide 
             + density + pH + sulphates + alcohol, data = train1, Hess = T)
summary(mod3)

#Compute p-values as before and visualize
coefs1 <- coef(summary(mod3))
p_val1 <- pnorm(abs(coefs1[,"t value"]), lower.tail = F)*2
coefs1 <- cbind(coefs1,"p-value" = p_val1)
coefs1

#Perform LRT and compute AIC-BIC in order to reach conclusion on chlorides
lrtest(mod2,mod3)
aicmod2 <- AIC(mod2) 
aicmod3 <- AIC(mod3) 
modaic <- cbind(rbind(aicmod2,aicmod3),rbind(BIC(mod2),BIC(mod3)))
rownames(modaic) <- c("Model with Chlorides", "Model without Chlorides")
colnames(modaic) <- c("AIC","BIC")
modaic

#Use test set to compute predicted classes and predicted probabilities, confussion matrix
#and misclassification error
predictedScores2 <- predict(mod3, test1, type="probs")
predictedClass2 <- predict(mod3, test1) 

#conf-matrix
conf_mat2 <- table(test1$quality, predictedClass2)
conf_mat2
#miss error
ms_err2 <- 1 - sum(diag(conf_mat2))/sum(conf_mat2)
ms_err2

#Use pred. probabilities based on training set for the plot
predictedScores1 <- predict(mod3, train1, type="probs")

#Needed fictitious values used in plot 
fix.acid <- seq(from = min(my_data$fixed.acidity),to = max(my_data$fixed.acidity), length.out = nrow(predictedScores1))
vol.acid <- seq(from = min(my_data$volatile.acidity), to = max(my_data$volatile.acidity), length.out = nrow(predictedScores1))
cit.acid <- seq(from = min(my_data$citric.acid), to = max(my_data$citric.acid), length.out = nrow(predictedScores1))
res.sug <- seq(from = min(my_data$residual.sugar),to = max(my_data$residual.sugar), length.out = nrow(predictedScores1))
tot.sl <- seq(from = min(my_data$total.sulfur.dioxide),to = max(my_data$total.sulfur.dioxide), length.out = nrow(predictedScores1))
dens <- seq(from = min(my_data$density),to = max(my_data$density), length.out = nrow(predictedScores1))
ph <- seq(from = min(my_data$pH),to = max(my_data$pH), length.out = nrow(predictedScores1))
sulph <- seq(from = min(my_data$sulphates),to = max(my_data$sulphates), length.out = nrow(predictedScores1))
alc <- seq(from = min(my_data$alcohol),to = max(my_data$alcohol), length.out = nrow(predictedScores1))
b_x <- fix.acid*coefs1[1,1] + vol.acid*coefs1[2,1] + cit.acid*coefs1[3,1] + res.sug*coefs1[4,1] + 
  tot.sl*coefs1[5,1] + dens*coefs1[6,1] + ph*coefs1[7,1] + sulph*coefs1[8,1] + alc*coefs1[9,1]

#Calculating the probabilities for each category
logfunc <- function(x) {
  return( 1/(1+exp(-x) ) )
}
p1 <- logfunc(coefs1[10,1] - b_x)
p2 <- logfunc(coefs1[11,1] - b_x) - logfunc(coefs1[10,1] - b_x)
p3 <- 1 - logfunc(coefs1[11,1] - b_x)

#Plotting the results
plot(b_x, p3, type='l', ylab='Prob', xlab = 'sum(beta*x)')
lines(b_x, p2, col='red')
lines(b_x, p1, col='blue')
legend("topleft", lty=1, col=c("black", "red", "blue"), 
       legend=c("category 7", "category 6", "category 5"))

