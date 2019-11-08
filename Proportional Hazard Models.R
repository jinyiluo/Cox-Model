---
  title: "Proportional Hazard Models: Estimation, Prediction and Evaluation"
output: html_notebook
---
# Install Library

library(survival)
library(survC1)
library(readxl)
library(nricens)
library(cmprsk)
library(dplyr)
library(survIDINRI)###
library(mstate)###
library(PredictABEL)###
IMRAWandISTnov <- read_excel("IMRAWandISTnov.xls")
View(IMRAWandISTnov)


# 1 Computation Problems
"
IMRAW(International Myelodysplasia Risk Analysis Workshop)
Ai=AGE
Si=GENDER, 0=female or 1=male
Ni=NEUTRO(ANC(Absolute Neutrophil Count))
Pi=PLATE(platelets)
Ii=IST treatment(Immunosuppressive Therapy), NIH=1
"

# Extract Dataset, where Zi=(Ai,Si,Ni,Pi,Ii)^T

IMRAW_IST<-na.omit(IMRAWandISTnov[,c("survival","DIED","AGE","GENDER","NEUTRO","PLATE","NIH")])
# Change list to double: Set 0 for female and 1 formale
data1<-model.matrix(~ survival + DIED + AGE + GENDER + NEUTRO + PLATE + NIH,data=IMRAW_IST)[,-1]

#1.1 Whether IST treatment is predictive for overall survival after adjusting Zi^0=(Ai,Si,Ni,Pi)^T
# 1.1.1 & 1.1.2

# Inference of Delta C
model111<-Inf.Cval.Delta(mydata=data1[,c("survival","DIED")], # mydata: Input data. The 1st column should be time-to-event, and the 2nd column is event indicator (1=event, 0=censor).
                         covs0 = data1[,c("AGE","GENDERM","NEUTRO","PLATE")], # Model 1
                         covs1 = data1[,c("AGE","GENDERM","NEUTRO","PLATE","NIH")], # Model 2
                         tau = 3*365.25, # Truncation time (3 x 365 + 366)/4 = 365.25 days.
                         itr = 100) # Iteration of perturbation-resampling

round(model111,4)

# 1.1.3 & 1.1.4

# Inference of Delta C
model113<-Inf.Cval.Delta(mydata=data1[,c("survival","DIED")], # mydata: Input data. The 1st column should be time-to-event, and the 2nd column is event indicator (1=event, 0=censor).
                         covs0 = data1[,c("AGE","GENDERM","NEUTRO","PLATE")], # Standard prediction Model 0 (Model 1): "AGE","GENDERM","NEUTRO","PLATE"
                         covs1 = data1[,c("AGE","GENDERM","NEUTRO","PLATE","NIH")], # New prediction Model 1 (Model 2): "AGE","GENDERM","NEUTRO","PLATE","NIH"
                         tau = 5*365.25, # Truncation time (3 x 365 + 366)/4 = 365.25 days.
                         itr = 100) # Iteration of perturbation-resampling

round(model113,4)

#2 Whether IST treatment is useful for changing the risk categories of death after adjusting   Zi^0=(Ai,Si,Ni,Pi)^T

# 2.1.1

data2<-as.data.frame(data1)
time = data2$survival
event = data2$DIED
# Fit Cox Model
# Standard prediction Model 1: "AGE","GENDERM","NEUTRO","PLATE"
mstd<-coxph(Surv(time, event) ~ ., data=data.frame(time, event,subset(data2, select = c("AGE","GENDERM","NEUTRO","PLATE"))), x=TRUE)
# New prediction Model 2: "AGE","GENDERM","NEUTRO","PLATE","NIH"
mnew<-coxph(Surv(time, event) ~ ., data=data.frame(time, event, subset(data2, select = c("AGE","GENDERM","NEUTRO","PLATE","NIH"))), x=TRUE)
# Calculation of the risk category the Net Reclassification Improvement (NRI) for tau = 3 years
nri1 <- nricens(mdl.std = mstd, mdl.new = mnew, t0 = 3*365.25, cut = c(0.1, 0.3),
                niter = 0, alpha = 0.05, msg=TRUE)
SE=((nri1$nri$Estimate[4]+nri1$nri$Estimate[5])/sum(event)+(nri1$nri$Estimate[6]+nri1$nri$Estimate[7])/(length(event)-sum(event)))^(1/2)
lower=nri1$nri$Estimate[1]-1.96*SE
upper=nri1$nri$Estimate[1]+1.96*SE
round(cbind(SE,lower,upper),4)

# 2.2.2

nri2 <- nricens(mdl.std = mstd, mdl.new = mnew, t0 = 3*365.25, cut = c(0.1, 0.3, 0.5),
                niter = 0, alpha = 0.05, msg=TRUE)
SE=((nri2$nri$Estimate[4]+nri2$nri$Estimate[5])/sum(event)+(nri2$nri$Estimate[6]+nri2$nri$Estimate[7])/(length(event)-sum(event)))^(1/2)
lower=nri2$nri$Estimate[1]-1.96*SE
upper=nri2$nri$Estimate[1]+1.96*SE
round(cbind(SE,lower,upper),4)

# 2.2.3

nri3 <- nricens(mdl.std = mstd, mdl.new = mnew, t0 = 5*365.25, cut = c(0.2, 0.4),
                niter = 0, alpha = 0.05, msg=TRUE)
SE=((nri3$nri$Estimate[4]+nri3$nri$Estimate[5])/sum(event)+(nri3$nri$Estimate[6]+nri3$nri$Estimate[7])/(length(event)-sum(event)))^(1/2)
lower=nri3$nri$Estimate[1]-1.96*SE
upper=nri3$nri$Estimate[1]+1.96*SE
round(cbind(SE,lower,upper),4)

# 2.2.4

nri4 <- nricens(mdl.std = mstd, mdl.new = mnew, t0 = 5*365.25, cut = c(0.2, 0.4, 0.6),
                niter = 0, alpha = 0.05, msg=TRUE)
SE=((nri4$nri$Estimate[4]+nri4$nri$Estimate[5])/sum(event)+(nri4$nri$Estimate[6]+nri4$nri$Estimate[7])/(length(event)-sum(event)))^(1/2)
lower=nri4$nri$Estimate[1]-1.96*SE
upper=nri4$nri$Estimate[1]+1.96*SE
round(cbind(SE,lower,upper),4)

# 3 Consider two competing events: "time to AML death" (Event 1) and "time to Non-AML death" (Event 2).

# 3.3.1

# Add "AML" and"AMLTIME" to the data
data3<-IMRAWandISTnov[,c("survival","DIED","AGE","GENDER","NEUTRO","PLATE","NIH","AML","AMLTIME")]
# Add an "event" column, and set 0 when DIED==0, 1 when DIED==1&AML==1, 2 when DIED==1&AML==0, otherwise 3.
data331_1<-data3 %>% mutate(event = ifelse(DIED==0,0,
                                           ifelse(DIED==1&AML==1,1,
                                                  ifelse(DIED==1&AML==0,2,3
                                                  )
                                           )
)
)
# Add a "ftime" column to represent failure time variable, and use pmin to get min for each survival vs AMLTIME
data331_2<-data331_1 %>% mutate(ftime=pmin(survival,AMLTIME))
data331_3<-data331_2[,c("NIH","event","ftime")]

# For t = 3 years, what are the estimates and their corresponding SE and 95% CI for F_1(3|IST) and F_2(3|IMRAW)?
fit<-cuminc(ftime=data331_3$ftime, fstatus=data331_3$event, group=data331_3$NIH)
tfit <- timepoints(fit, times=c(3*365.25))
conf.level=0.95
z <- qnorm(1-(1-conf.level)/2)
lower <- tfit$est ^ exp(-z*sqrt(tfit$var)/(tfit$est*log(tfit$est)))
upper <- tfit$est ^ exp(z*sqrt(tfit$var)/(tfit$est*log(tfit$est)))
t1<-round(matrix(c(tfit$est[2],sqrt(tfit$var)[2],lower[2],upper[2],tfit$est[3],sqrt(tfit$var)[3],lower[3],upper[3]),ncol=4,byrow=TRUE),4)
colnames(t1)<-c("Estimate","SE","Lower","Upper")
rownames(t1)<-c("F1","F2")
t1


# 3.3.2

fit2<-cuminc(ftime=data331_3$ftime, fstatus=data331_3$event, group=data331_3$NIH)
tfit2 <- timepoints(fit2, times=c(5*365.25))
conf.level=0.95
z <- qnorm(1-(1-conf.level)/2)
lower <- tfit2$est ^ exp(-z*sqrt(tfit2$var)/(tfit2$est*log(tfit2$est)))
upper <- tfit2$est ^ exp(z*sqrt(tfit2$var)/(tfit2$est*log(tfit2$est)))
t2<-round(matrix(c(tfit2$est[2],sqrt(tfit2$var)[2],lower[2],upper[2],tfit2$est[3],sqrt(tfit2$var)[3],lower[3],upper[3]),ncol=4,byrow=TRUE),4)
colnames(t2)<-c("Estimate","SE","Lower","Upper")
rownames(t2)<-c("F1","F2")
t2

# 3.3.3

data333_1<-na.omit(data331_2[,c("AGE","GENDER","NEUTRO","PLATE","NIH","event","ftime")])
data333_2<-data333_1 %>% mutate(GENDERM= ifelse(GENDER=="F",0,1)) # Change F=0, M=1
data333_3<-data333_2[,c("AGE","GENDERM","NEUTRO","PLATE","NIH")]
# Competing Risks Regression
c1<-crr(ftime=data333_2$ftime, fstatus=data333_2$event, cov1=data333_3, failcode=1, cencode=0)
summary(c1)
round(summary(c1)$coef[,c(1,2,4,5)],4)

# 3.3.4

# Competing Risks Regression
c2<-crr(ftime=data333_2$ftime, fstatus=data333_2$event, cov1=data333_3, failcode=2, cencode=0)
summary(c2)
round(summary(c2)$coef[,c(1,2,4,5)],4)














