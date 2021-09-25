rm(list=ls())
require('logspline')
require('copula')
require('cubature')
source('nume_util.r')
source('dr_func.r')
source('pr_boot.r')
#山东产量数据
dat<-read.csv('yield_sd.csv')
N_year<-length(unique(dat$Year))
N_county<-length(unique(dat$County))
cty_seq<-unique(dat$County)
m_cty<-c()
sd_cty<-c()
std_val<-dat$Value
for (cty in cty_seq){
  temp<-which(dat$County==cty)
  m_tp<-mean(std_val[temp])
  m_cty<-c(m_cty,m_tp)
  std_val[temp]<-(std_val[temp]-m_tp)/m_tp
  sd_cty<-c(sd_cty,sd(std_val[temp]))
}
hist(std_val,freq=F,ylim=c(0,8),main='',xlab=expression(epsilon/bar(Y)))
LB<-min(std_val)*2
UB<-max(std_val)*2
fbase<-logspline(std_val,lbound=LB,ubound=UB)
pit_seq<-plogspline(std_val,fbase)
lines(sort(std_val),dlogspline(sort(std_val),fbase),col=1)

res<-c()
avgApr<-rep(0,12)
avgOct<-rep(0,12)
# 山东各市价格数据
hd<-names(dat_price<-read.csv('soybean_price/d1-2006.csv'))
for (i in 2006:2017){
  dat_price<-read.csv(paste('soybean_price/d1-',i,'.csv',sep=''))
  names(dat_price)<-hd
  id_sel1<-dat_price$合约==paste('a',substr(paste(i+1),3,4),'01',sep='')&((dat_price$日期>=i*10000+0401&dat_price$日期<=i*10000+0430))
  id_sel2<-dat_price$合约==paste('a',substr(paste(i+1),3,4),'01',sep='')&((dat_price$日期>=i*10000+1001&dat_price$日期<=i*10000+1031))
  avgApr[i-2005]<-mean(dat_price$结算价[id_sel1])
  avgOct[i-2005]<-mean(dat_price$结算价[id_sel2])
  res<-rbind(res,dat_price[id_sel1|id_sel2,])
  print(i)
}
price_seq<-apply(cbind(avgApr,avgOct),1,max)
std_p<-price_seq/mean(price_seq)
p_pit<-plnorm(std_p,mean(log(std_p)),sd(log(std_p)))
# Cramer-von Mises Statistics Scores
sn_mat<-matrix(0,N_county,3)
hist(std_val,freq=F,ylim=c(0,12.5),main='',xlab=expression(epsilon/bar(Y)))
lines(seq(LB,UB,length=300),dlogspline(seq(LB,UB,length=300),fbase),type='o')
r_mat<-c()
for (i in 1:N_county){
  d_tempp<-function(x){dlogspline(x,fbase)*dnorm(x,0,sd_cty[i])}
  c_inte<-integrate(d_tempp,LB,UB)$value
  d_temp<-function(x){d_tempp(x)/c_inte}
  lines(seq(LB,UB,length=300),d_temp(seq(LB,UB,length=300)),col=i+1)
  p_temp<-function(x) sapply(x,function(x) integrate(d_temp,LB,x)$value)
  y_pit_i<-p_temp(std_val[seq(i,length(std_val),N_county)])
  temp<-cbind(y_pit_i,p_pit)
  sn_mat[i,1]<-cal_sn(normalCopula(),temp)
  sn_mat[i,2]<-cal_sn(frankCopula(),temp)
  sn_mat[i,3]<-cal_sn(claytonCopula(),temp,'itau')
  r_mat<-rbind(r_mat,rank(sn_mat[i,]))
  cat(i,'  ')
}
legend('topright',legend='logspline',col=1,lty=1,pch=1)
# number of times when each copula has max statistics
table(apply(sn_mat,1,function(x) which(x==max(x))))
# number of times when each copula has min statistics
table(apply(sn_mat,1,function(x) which(x==min(x))))

# pr_seq，pr_seq2？
pr_seq<-matrix(0,7,N_county)
pr_seq2<-pr_seq
for (i in 1:N_county){
  d_tempp<-function(x){dlogspline(x,fbase)*dnorm(x,0,sd_cty[i])}
  c_inte<-integrate(d_tempp,LB,UB)$value
  d_temp<-function(x){d_tempp(x)/c_inte}
  p_temp<-function(x) sapply(x,function(x) integrate(d_temp,LB,x)$value)
  y_pit_i<-p_temp(std_val[seq(i,length(std_val),N_county)])
  res_cpl<-fitCopula(claytonCopula(),cbind(y_pit_i,p_pit),'itau')
  for (j in 1:7){
    pr_seq2[j,i]<-premium_rate_yp(res_cpl@copula,p_temp,d_temp,function(x) plnorm(x,mean(log(std_p)),sd(log(std_p))),function(x) dlnorm(x,mean(log(std_p)),sd(log(std_p))),LB,UB,0,max(std_p)*2,1-(j-1)*0.05)
    pr_seq[j,i]<-pr_bootstrap(res_cpl@copula,q_f(p_temp,LB,UB),function(p) qlnorm(p,mean(log(std_p)),sd(log(std_p))),1-(j-1)*0.05)
  }
  # pr_seq2[i]<-premium_rate_y(d_temp,LB,UB,.9)
  cat(i,'  ')
}
print(pr_seq)
print(pr_seq2)



# N_cell<-round(N_year/2)
# freq_mat<-val2freq(pit_seq,dat$County,N_cell)
# bsfunc_list<-gene_bsfunc_list(6,legendre_basis,0,1)
# preg_res<-possion_reg(freq_mat,N_cell,bsfunc_list)
# pr_seq<-rep(0,N_county)
# # Cramer-von Mises Statistics Scores
# sn_mat<-matrix(0,N_county,3)
# for (i in 1:N_county){
#   d_temp<-d_indi(preg_res$b[i,],function(x)dlogspline(x,fbase),function(x)plogspline(x,fbase),bsfunc_list,-1,max(std_val)*3)
#   p_temp<-function(x) sapply(x,function(x) integrate(d_temp,-1,x)$value)
#   print(ks.test(std_val[seq(i,204,by=17)],pnorm,sd=sd_cty[i]))
#   print(ks.test(std_val[seq(i,204,by=17)],p_temp))
#   y_pit_i<-p_temp(std_val[seq(i,length(std_val),N_county)])
#   sn_mat[i,1]<-cal_sn(normalCopula(),cbind(y_pit_i,p_pit))
#   sn_mat[i,2]<-cal_sn(amhCopula(),cbind(y_pit_i,p_pit))
#   sn_mat[i,3]<-cal_sn(claytonCopula(),cbind(y_pit_i,p_pit),'itau')
#   cat(i,'  ')
# }



# hist(std_p,freq=F,main='',xlab=expression(P/bar(P)))
# lines(seq(0.6,1.2,by=0.01),dlnorm(seq(0.6,1.2,by=0.01),mean(log(std_p)),sd(log(std_p))),col=2)
# legend('topleft','log-normal',col=2,lty=1)