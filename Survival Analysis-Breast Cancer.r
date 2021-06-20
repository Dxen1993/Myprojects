#Case Study Assignment 5
#QUESTION 1---Median follow-up
rm(list=ls())
install.packages("TH.data")
library(TH.data)
data(GBSG2)
head(GBSG2)
help(GBSG2)

library(survival)
install.packages("KMsurv")
library(KMsurv)

#Creating a variable status wich is the cens reversed because in follow-up time death times are treated as censored 
GBSG2$status = rep(0,nrow(GBSG2))
for(i in 1:nrow(GBSG2)){
  if(GBSG2$cens[i] == 1){
    GBSG2$status[i] <- 0
  }
  else{
    GBSG2$status[i] <- 1
  }
}

sf <- survfit(Surv(time,status) ~ 1 , GBSG2)
sf

sf2 <- survfit(Surv(time, cens) ~ 1, GBSG2)
par(mfrow=c(1,2))
plot(sf2, main = "Kaplan-Meier original", col = "red",xlab ="Time(Days)",ylab ="Survival Probability")
#Now plot the reverse Kaplan-Meier plot
plot(sf,main="Kaplan-Meier follow-up",col="blue", xlab = "Time(Days)" ,ylab ="Probability of not being censored")
summary(sf)
#interpretation...

#QUESTION2--Kaplan Meier for the 2 treatment groups + Log-rank test. Time dependent?? Truncated??

data_therapy <- GBSG2[which(GBSG2$horTh == "yes"), ]
data_notherapy <- GBSG2[which(GBSG2$horTh == "no"), ]
therapy <- survfit(Surv(time,cens) ~ 1 , data_therapy)
therapy
notherapy <- survfit(Surv(time,cens) ~ 1, data_notherapy)
notherapy

#Kaplan-Meier plot with 95% confidence interval
plot(therapy,main = "Kaplan-Meier for the two treatment groups",col="blue")
lines(notherapy,col="red")

#OR JUST
par(mfrow = c(1,1))
plot(survfit(Surv(time,cens) ~ horTh, data = GBSG2), col=1:2 , main="Kaplan-Meier plot of the two treatment groups")
legend("topright",c("Therapy","No therapy"),lwd=1,col=2:1)

#Conclusion: Hormonal treatment improves probability of survival

#Log-Rank test
survdiff(Surv(time,cens) ~ horTh, data = GBSG2) 
#There is a significant difference in survival between the 2 groups

#Cumulative hazards
fit1 <- summary(therapy)
fit2 <- summary(notherapy)
Hhat1 <- -log(fit1$surv)
Hhat2 <- -log(fit2$surv)
plot(Hhat1,main = "cummulative hazards",col="blue")
lines(Hhat2,col="red")

#Hazard Ratio
summary(coxph(Surv(time,cens) ~ horTh, GBSG2))
#The hazard ratio of a non-therapy patient relative to a therapy patient is exp(coef) = 0.6949 .
#A patient that doesn't receive therapy has 30% increase of risk of dying than a therapy-patient .

#QUESTION 3 --Cox Regression with forward selection
#TIED DATA ?

#step1
coxph(Surv(time,cens) ~ age, GBSG2) #p-value = 0.446 -> not significant
coxph(Surv(time,cens) ~ menostat, GBSG2) #pvalue = 0.595 -> not significant
coxph(Surv(time,cens) ~ tsize, GBSG2) #p-value = 2.3e-0.5 
coxph(Surv(time,cens) ~ factor(tgrade), GBSG2) #p-value = 1e-0.5, 0.047
coxph(Surv(time,cens) ~ pnodes, GBSG2) #p-value < 2e-16
coxph(Surv(time,cens) ~ progrec, GBSG2) #p-value = 1.5e-0.6
coxph(Surv(time,cens) ~ estrec, GBSG2) #p-value = 0.041

#the covariate pnodes is the most significant here. Hence, it will be our first term in our 
#model.

#step2
coxph(Surv(time,cens) ~ pnodes + tsize, GBSG2) #p-value = 0.067-> not significant
coxph(Surv(time,cens) ~ pnodes + factor(tgrade), GBSG2) #p-value = 6e-05, 0.099
coxph(Surv(time,cens) ~ pnodes + progrec, GBSG2) #p-value = 3.7e-0.6
coxph(Surv(time,cens) ~ pnodes + estrec, GBSG2) #p-value = 0.078 -> not significant

#next, we include the covariate progrec into our model

#step3
final.model1 <- coxph(Surv(time,cens) ~ pnodes + progrec + tsize, GBSG2) #p-value = 0.047
coxph(Surv(time,cens) ~ pnodes + progrec + factor(tgrade), GBSG2) #p-value = 0.0016, 0.1212

#step4
final.model2 <- coxph(Surv(time,cens) ~ pnodes + progrec + tsize + factor(tgrade), GBSG2)
summary(final.model2)
anova(final.model1, final.model2)

#step5
final.model3 <- coxph(Surv(time,cens) ~ pnodes + progrec + tsize + factor(tgrade)+ horTh, GBSG2)
anova(final.model3, final.model2)

summary(final.model3)

#QUESTION 4 -- Risk Scores Caluculation 
coeff <- final.model2$coefficients
nonames.coeff <-unname(coeff[c(3, 1, 2)])
as.vector(nonames.coeff)
new.dataset <- GBSG2[, c( 4, 6, 7)]


gradeInd <- rep(0, nrow(new.dataset))
for(i in 1:nrow(GBSG2)){
  if((GBSG2$tgrade[i] == "II") == TRUE)
    gradeInd[i] <- 0.583887363
  else if((GBSG2$tgrade[i] == "III") == TRUE)
    gradeInd[i] <-  -0.188032768 
  else gradeInd[i] <- 0
}
  
risk.score <- as.matrix(new.dataset) %*% as.matrix(nonames.coeff) + gradeInd
risk.score
hist(risk.score)

#5
quantile(risk.score, c(1/3, 2/3, 1))
summary(risk.score)
risk.group <- rep(0, nrow(GBSG2))
for(i in 1:nrow(GBSG2)){
  if(risk.score[i] < 0.3){ 
    risk.group[i] <- 1
    }
  else if(0.3 <= risk.score[i] && risk.score[i] < 0.8){
    risk.group[i] <- 2
  }
  else{ risk.group[i] <- 3}
}
table(risk.group)
GBSG2_new <- cbind(GBSG2, risk.group, risk.score)

#Kaplan-Meier plot of the risk groups
plot(survfit(Surv(time,cens) ~ risk.group, data = GBSG2_new), col=1:3 , xlab = "Time(Days)", ylab = "Survival Probability", main="Kaplan-Meier plot of the risk groups")
legend("bottomleft",c("low","middle", "high"),lwd=1,col=1:3)

#log rank test
survdiff(Surv(time,cens) ~ factor(risk.group), data = GBSG2_new) 

#hazard ratio
summary(coxph(Surv(time,cens) ~ factor(risk.group), GBSG2_new))


