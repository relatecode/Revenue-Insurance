#循环（pr_mat模型四结果和不循环时不同(但是对照了一下p_T和y_T是完全一样的呀)，破案了！pr_seq的式子y_T没加i
#但是sn_mat结果却一致）
rm(list=ls())
library(ADGofTest)
library(readxl)
library(MASS)
library(splines)
require(copula)

# Cramer-von Mises Statistics
cal_sn<-function(copula,data,...){
  N_obs<-nrow(data)
  cfit<-fitCopula(copula,data,...)
  cn<-rep(0,N_obs)
  for (i in 1:N_obs){
    cn[i]<-sum((data[,1]<=data[i,1])&(data[,2]<=data[i,2]))
  }   #sum计算复合条件的数据个数
  cn<-cn/N_obs
  ctheta<-pCopula(data,cfit@copula)
  return(sum((cn-ctheta)^2))
}

pr_bootstrap<-function(copula,qy,qp,lambda,B=10000){
  r_cp<-rCopula(B,copula)   #random generation,生成B*2矩阵，Q2生成的本来就应该是分位数而不是概率呀，为什么要使用分位数函数呢
  y_r<-rep(0,B)
  p_r<-y_r
  for (i in 1:B){
    temp<<-r_cp[i,1]
    y_r[i]<-qy(r_cp[i,1])
    p_r[i]<-qp(r_cp[i,2])*mean(price_seq)
    if (floor(i/B*10)!=floor((i+1)/B*10))
      print(paste((floor(i/B*10)+1)*10,'%',sep=''))  
  }
  er<-mean((y_r)*(p_r))   
  L<-mean(sapply((lambda*er-(y_r)*(p_r)),max,0))
  return(L/er/lambda)
}

dat<-read_excel('山东省大豆数据.xlsx',sheet="Sheet2")
year<-dat$year-min(dat$year)+1
T<-length(year)
rf<-dat$rainfall
y<-dat$yield
p<-dat$price
im<-dat$import

knots_rf <- quantile(rf, c(1/4, 2/4, 3/4))
rf_knot <- bs(rf, knots = knots_rf, degree = 1)

knots_im <- quantile(im, c(1/4, 2/4, 3/4))
im_knot <- bs(im, knots = knots_im, degree = 1)

y_T<-matrix(NA,T,4)
p_T<-matrix(NA,T,4)
#模型一
y_fit1<-lm(y~year)
ycoef1<-y_fit1$coefficients
y_T[,1]<-y_fit1$residuals+ycoef1[1]+ycoef1[2]*year[T]

p_fit1<-lm(p~year)
pcoef1<-p_fit1$coefficients
p_T[,1]<-p_fit1$residuals+pcoef1[1]+pcoef1[2]*year[T]

#模型二
y_fit2<-lm(y~year+rf_knot)
ycoef2<-y_fit2$coefficients
y_T[,2]<-y_fit2$residuals+ycoef2[1]+ycoef2[2]*year[T]+ycoef2[3]*rf_knot[T,1]+ycoef2[4]*rf_knot[T,2]+ycoef2[5]*rf_knot[T,3]+ycoef2[6]*rf_knot[T,4]

p_fit2<-lm(p~year)
pcoef2<-p_fit2$coefficients
p_T[,2]<-p_fit2$residuals+pcoef2[1]+pcoef2[2]*year[T]

#模型三
y_fit3<-lm(y~year)
ycoef3<-y_fit3$coefficients
y_T[,3]<-y_fit3$residuals+ycoef3[1]+ycoef3[2]*year[T]

p_fit3<-lm(p~year+im_knot)
pcoef3<-p_fit3$coefficients
p_T[,3]<-p_fit3$residuals+pcoef3[1]+pcoef3[2]*year[T]+pcoef3[3]*im_knot[T,1]+pcoef3[4]*im_knot[T,2]+pcoef3[5]*im_knot[T,3]+pcoef3[6]*im_knot[T,4]

#模型四
y_fit4<-lm(y~year+rf_knot)
ycoef4<-y_fit4$coefficients
y_T[,4]<-y_fit4$residuals+ycoef4[1]+ycoef4[2]*year[T]+ycoef4[3]*rf_knot[T,1]+ycoef4[4]*rf_knot[T,2]+ycoef4[5]*rf_knot[T,3]+ycoef4[6]*rf_knot[T,4]

p_fit4<-lm(p~year+im_knot)
pcoef4<-p_fit4$coefficients
p_T[,4]<-p_fit4$residuals+pcoef4[1]+pcoef4[2]*year[T]+pcoef4[3]*im_knot[T,1]+pcoef4[4]*im_knot[T,2]+pcoef4[5]*im_knot[T,3]+pcoef4[6]*im_knot[T,4]
#
cpl_para<-rep(0,4)
sn_mat<-matrix(0,4,3)
pr_seq<-matrix(0,7,4)
for (i in 1:4){
  price_seq<-p_T[,i]
  std_p<-price_seq/mean(price_seq)   #计算Pt tuta
  p_pit<-plnorm(std_p,mean(log(std_p)),sd(log(std_p))) 
  y_pit_i<-pnorm(y_T[,i],mean(y_T[,i]),sd(y_T[,i]))
  temp<-cbind(y_pit_i,p_pit)
  sn_mat[i,1]<-cal_sn(normalCopula(),temp)  #cal_sn是pr_boot中定义的函数,输入copula函数和边缘分布的概率（累积概率，返回Cramer-von Mises Statistics Scores
  sn_mat[i,2]<-cal_sn(frankCopula(),temp)
  sn_mat[i,3]<-cal_sn(claytonCopula(),temp,'itau')
  res_cpl<-fitCopula(frankCopula(),temp)
  cpl_para[i]<-res_cpl@estimate
  for (j in 1:7)
    pr_seq[j,i]<-pr_bootstrap(res_cpl@copula,function(q_f1) qnorm(q_f1,mean(y_T[,i]),sd(y_T[,i])),function(p) qlnorm(p,mean(log(std_p)),sd(log(std_p))),1-(j-1)*0.05)
  cat(i,'  ')
}

par(oma=c(0.1,0.1,0.1,0.1)) 
plot(year+2007,y,type = "o",col=1,ylim=c(1000,3000),xlab='年',ylab='产量')
lines(year+2007,y_T[,1],type = "o",col='red')
lines(year+2007,y_T[,2],type = "o",col='blue')
legend('bottomright',legend = c('原产量','模型一/三','模型二/四'),col=c(1,'red','blue'),lty=1,pch=1)


plot(year+2007,p,type = "o",col=1,ylim=c(2000,5000),xlab='年',ylab='价格')
lines(year+2007,p_T[,1],type = "o",col='red')
lines(year+2007,p_T[,3],type = "o",col='blue')
legend('bottomright',legend = c('原价格','模型一/二','模型三/四'),col=c(1,'red','blue'),lty=1,pch=1)

# y_T[,1]<-y_fit4$residuals+ycoef4[1]+ycoef4[2]*year[T]+ycoef4[3]*rf_knot[which(rf==max(rf)),1]+ycoef4[4]*rf_knot[which(rf==max(rf)),2]+ycoef4[5]*rf_knot[which(rf==max(rf)),3]+ycoef4[6]*rf_knot[which(rf==max(rf)),4]
# y_T[,2]<-y_fit4$residuals+ycoef4[1]+ycoef4[2]*year[T]+ycoef4[3]*rf_knot[which(rf==min(rf)),1]+ycoef4[4]*rf_knot[which(rf==min(rf)),2]+ycoef4[5]*rf_knot[which(rf==min(rf)),3]+ycoef4[6]*rf_knot[which(rf==min(rf)),4]
# y_T[,3]<-y_fit4$residuals+ycoef4[1]+ycoef4[2]*year[T]+ycoef4[3]*rf_knot[T,1]+ycoef4[4]*rf_knot[T,2]+ycoef4[5]*rf_knot[T,3]+ycoef4[6]*rf_knot[T,4]
# y_T[,4]<-y_fit4$residuals+ycoef4[1]+ycoef4[2]*year[T]+ycoef4[3]*rf_knot[T,1]+ycoef4[4]*rf_knot[T,2]+ycoef4[5]*rf_knot[T,3]+ycoef4[6]*rf_knot[T,4]
# 
# p_T[,1]<-p_fit4$residuals+pcoef4[1]+pcoef4[2]*year[T]+pcoef4[3]*im_knot[T,1]+pcoef4[4]*im_knot[T,2]+pcoef4[5]*im_knot[T,3]+pcoef4[6]*im_knot[T,4]
# p_T[,2]<-p_fit4$residuals+pcoef4[1]+pcoef4[2]*year[T]+pcoef4[3]*im_knot[T,1]+pcoef4[4]*im_knot[T,2]+pcoef4[5]*im_knot[T,3]+pcoef4[6]*im_knot[T,4]
# p_T[,3]<-p_fit4$residuals+pcoef4[1]+pcoef4[2]*year[T]+pcoef4[3]*im_knot[which(im==max(im)),1]+pcoef4[4]*im_knot[which(im==max(im)),2]+pcoef4[5]*im_knot[which(im==max(im)),3]+pcoef4[6]*im_knot[which(im==max(im)),4]
# p_T[,4]<-p_fit4$residuals+pcoef4[1]+pcoef4[2]*year[T]+pcoef4[3]*im_knot[which(im==min(im)),1]+pcoef4[4]*im_knot[which(im==min(im)),2]+pcoef4[5]*im_knot[which(im==min(im)),3]+pcoef4[6]*im_knot[which(im==min(im)),4]


# rf[11]<-120
# rf[12]<-40
# im[11]<-12000
# im[12]<-2000
# knots_rf <- quantile(rf, c(1/4, 2/4, 3/4))
# rf_knot <- bs(rf, knots = knots_rf, degree = 1)
# knots_im <- quantile(im, c(1/4, 2/4, 3/4))
# im_knot <- bs(im, knots = knots_im, degree = 1)
# 
# y_T[,1]<-y_fit4$residuals+ycoef4[1]+ycoef4[2]*year[T]+ycoef4[3]*rf_knot[T+1,1]+ycoef4[4]*rf_knot[T+1,2]+ycoef4[5]*rf_knot[T+1,3]+ycoef4[6]*rf_knot[T+1,4]
# y_T[,2]<-y_fit4$residuals+ycoef4[1]+ycoef4[2]*year[T]+ycoef4[3]*rf_knot[T+2,1]+ycoef4[4]*rf_knot[T+2,2]+ycoef4[5]*rf_knot[T+2,3]+ycoef4[6]*rf_knot[T+2,4]
# y_T[,3]<-y_fit4$residuals+ycoef4[1]+ycoef4[2]*year[T]+ycoef4[3]*rf_knot[T,1]+ycoef4[4]*rf_knot[T,2]+ycoef4[5]*rf_knot[T,3]+ycoef4[6]*rf_knot[T,4]
# y_T[,4]<-y_fit4$residuals+ycoef4[1]+ycoef4[2]*year[T]+ycoef4[3]*rf_knot[T,1]+ycoef4[4]*rf_knot[T,2]+ycoef4[5]*rf_knot[T,3]+ycoef4[6]*rf_knot[T,4]
# 
# p_T[,1]<-p_fit4$residuals+pcoef4[1]+pcoef4[2]*year[T]+pcoef4[3]*im_knot[T,1]+pcoef4[4]*im_knot[T,2]+pcoef4[5]*im_knot[T,3]+pcoef4[6]*im_knot[T,4]
# p_T[,2]<-p_fit4$residuals+pcoef4[1]+pcoef4[2]*year[T]+pcoef4[3]*im_knot[T,1]+pcoef4[4]*im_knot[T,2]+pcoef4[5]*im_knot[T,3]+pcoef4[6]*im_knot[T,4]
# p_T[,3]<-p_fit4$residuals+pcoef4[1]+pcoef4[2]*year[T]+pcoef4[3]*im_knot[T+1,1]+pcoef4[4]*im_knot[T+1,2]+pcoef4[5]*im_knot[T+1,3]+pcoef4[6]*im_knot[T+1,4]
# p_T[,4]<-p_fit4$residuals+pcoef4[1]+pcoef4[2]*year[T]+pcoef4[3]*im_knot[T+2,1]+pcoef4[4]*im_knot[T+2,2]+pcoef4[5]*im_knot[T+2,3]+pcoef4[6]*im_knot[T+2,4]
