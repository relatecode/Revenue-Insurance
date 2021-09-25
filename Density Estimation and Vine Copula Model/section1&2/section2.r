# R code for 'Calculation of Crop Insurance Premium:Based on a Density Ratio Model and Vine Copula¡®
# section 2
# take jinan for example
Y=res[dat$City.code ==1]
CHS = dat$CHS[dat$City.code ==1]
library(readxl)
Price <- read_excel("D:/Y/´ó¶¹/Price.xlsx")
P = Price$Price

cor.test(Y,CHS,method='kendall')
cor.test(Y,P,method='kendall')
cor.test(P,CHS,method='kendall')

#import the real data and find the GoF of distribution P, Y, CHS
library(fitdistrplus)
fit.weibull<-fitdist(P,'weibull')
summary(fit.weibull)
fit.lnorm<-fitdist(P,'lnorm')
summary(fit.lnorm)
fit.gamma<-fitdist(P,"gamma")
summary(fit.gamma)
fit.weibull<-fitdist(CHS,'weibull')
summary(fit.weibull)
fit.lnorm<-fitdist(CHS,'lnorm')
summary(fit.lnorm)
fit.gamma<-fitdist(CHS,"gamma")
summary(fit.gamma)

ks.test(P,'pweibull',9.377398,4.180464)
ks.test(CHS,'pweibull',2.991485,545.184940)

PU<-c()
PU=pweibull(P,9.377398,4.180464)
CU<-c()
CU=pweibull(CHS,2.991485,545.184940)

library(copula)
library(MASS)
library(VineCopula)
library("fitdistrplus")

  

BiCopCompare(YU[,1],PU, familyset=1:6)

BiCopCompare(YU[,1],CU, familyset=1:6)


# calculate F_{1|2}(u1|u2) F_{3|2}(u3|u2) and Copula_13|2
c12<-BiCop(3,0.06)
c32<-BiCop(1,0.09)
F_12 = BiCopHinv(YU[,1],PU,c12)$hinv2
F_32 = BiCopHinv(YU[,1],CU ,c32)$hinv2
BiCopCompare(F_12,F_32 , familyset=1:6)




##simulation pairs of (U1,U2,U3)
family <- c(3,1,1)
par <- c(0.06,0.09, 0.99)          
par2<-c(0,0,0)
matrix<-c(1,2,3)
RVM <-C2RVine(matrix,family,par,par2)
set.seed(1)
U<- RVineSim(10000, RVM)

##simulation P,Y,CHS
P2 = qweibull(U[,1], 9.377398,4.180464)
CHS2 = qweibull(U[,3],2.991485,545.184940)
Y= dat$Yield[dat$City.code ==1]
f000=logspline(res,lbound=min(res),ubound=max(res))
res2=qlogspline(U[,2],f000)
Y2=(res2+1)*mean(Y)
x<-cbind(P2,Y2,CHS2)
##standard rainfall index CHA
alpha =0.2
CHA = qweibull(alpha,2.991485,545.184940)
CHA
##expected revenue ER
ER=mean(dat$Yield[dat$City.code ==1])*mean(P)
##revenue insurance
##level=1
ft= function(p,y){
  Ip=ER
  I = p*y
  (Ip-I)*(Ip-I>0)
}
M<-ft(P2,Y2)
mean(M)

##level=0.95
ft= function(p,y){
  Ip=ER*0.95
  I = p*y
  (Ip-I)*(Ip-I>0)
}
M<-ft(P2,Y2)
mean(M)

##level=0.9
ft= function(p,y){
  Ip=ER*0.9
  I = p*y
  (Ip-I)*(Ip-I>0)
}
M<-ft(P2,Y2)
mean(M)
