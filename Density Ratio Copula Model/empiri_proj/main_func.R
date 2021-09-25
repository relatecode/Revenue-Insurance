rm(list=ls())
require('logspline')
require('copula')
require('pracma')
source('nume_util.r')
source('dr_func.r')
source('price_process.r')
dat<-read.csv('yield_sd.csv')
N_year<-length(unique(dat$Year))
N_county<-length(unique(dat$County))
std_val<-yield_preprocess(dat$Value,dat$Year,dat$County)
# base density
fbase<-logspline(std_val,lbound=min(std_val)*1.1,ubound=max(std_val)*1.1)
pit_seq<-plogspline(std_val,fbase)
N_cell<-round(N_year/2)
freq_mat<-val2freq(pit_seq,dat$County,N_cell)
bsfunc_list<-gene_bsfunc_list(6,legendre_basis,0,1)
preg_res<-possion_reg(freq_mat,N_cell,bsfunc_list)
pr_seq<-rep(0,N_county)
# Cramer-von Mises Statistics Scores
sn_mat<-matrix(0,N_county,3)
for (i in 1:N_county){
  d_temp<-d_indi(preg_res$b[i,],function(x)dlogspline(x,fbase),function(x)plogspline(x,fbase),bsfunc_list,-1,max(std_val)*3)
  p_temp<-function(x) sapply(x,function(x) integrate(d_temp,-1,x)$value)
  y_pit_i<-p_temp(std_val[seq(i,length(std_val),N_county)])
  sn_mat[i,1]<-cal_sn(normalCopula(),cbind(y_pit_i,p_pit))
  sn_mat[i,2]<-cal_sn(amhCopula(),cbind(y_pit_i,p_pit))
  sn_mat[i,3]<-cal_sn(claytonCopula(),cbind(y_pit_i,p_pit),'itau')
  cat(i,'  ')
}
# number of times when each copula has max statistics
table(apply(sn_mat,1,function(x) which(x==max(x))))
# number of times when each copula has min statistics
table(apply(sn_mat,1,function(x) which(x==min(x))))
# Calculate Premium Rate
for (i in 1:N_county){
  d_temp<-d_indi(preg_res$b[i,],function(x)dlogspline(x,fbase),function(x)plogspline(x,fbase),bsfunc_list,-1,max(std_val)*3)
  p_temp<-function(x) sapply(x,function(x) integrate(d_temp,-1,x)$value)
  y_pit_i<-p_temp(std_val[seq(i,length(std_val),N_county)])
  res_cpl<-fitCopula(normalCopula(),cbind(y_pit_i,p_pit))
  pr_seq[i]<-premium_rate_yp(res_cpl@copula,p_temp,d_temp,function(x) pmixlnorm(x,mix_res),function(x) dmixlnorm(x,mix_res),4,.9)
  cat(i,'  ')
}
