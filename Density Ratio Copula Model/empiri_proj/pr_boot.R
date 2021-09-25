require(fitdistrplus)
require(copula)
##### MAIN FUNCTION #####
pr_cal<-function(yield_pit,qyield,price,lambda,verbose=T){
  # price data must be preprocessed
  # yield data must be on [0,1]
  # qyield is a quantile function of yield model
  nfseq<-c('lnorm','weibull')
  pfseq<-c(plnorm,pweibull)
  N_pmod<-length(pfseq)
  qfseq<-c(qlnorm,qweibull)
  res<-rep(0,N_pmod)
  # AIC for price model
  for (i in 1:N_pmod)
    res[i]<-fitdist(price,nfseq[i])$aic
  id_p<-which(res==min(res))[1]
  if (verbose)
    cat(nfseq[id_p],'for price model','\n')
  res<-fitdist(price,nfseq[id_p])$estimate
  N_pm<-length(res)
  str_pm<-''
  for (i in 1:N_pm)
    str_pm<-paste(str_pm,res[i],sep=',')
  # print(paste('price_pit<-pfseq[[id_p]](price',str_pm,')',sep=''))
  eval(parse(text=paste('price_pit<-pfseq[[id_p]](price',str_pm,')',sep='')))
  temp<-cbind(yield_pit,price_pit)
  # CM statistics for copula
  sn_seq<-rep(0,3)
  ncplseq<-c('norm','frank','clayton')
  cplseq<-c(normalCopula(),frankCopula(),claytonCopula())
  mlist<-c('mpl','mpl','itau')
  for (i in 1:3)
    sn_seq[i]<-cal_sn(cplseq[[i]],temp,mlist[i])
  id_c<-which(sn_seq==min(sn_seq))[1]
  if (verbose)
    cat(ncplseq[id_c],'for copula model','\n')
  res_cpl<-fitCopula(cplseq[[id_c]],cbind(y_pit_i,p_pit),mlist[id_c])
  # print(paste('qprice<-function(p){qfseq[[id_p]](p',str_pm,')}',sep=''))
  eval(parse(text=paste('qprice<-function(p){qfseq[[id_p]](p',str_pm,')}',sep='')))
  return(pr_bootstrap(res_cpl@copula,qyield,qprice,lambda))
}
pr_bootstrap<-function(copula,qy,qp,lambda,B=10000){
  r_cp<-rCopula(B,copula)
  y_r<-rep(0,B)
  p_r<-y_r
  for (i in 1:B){
    temp<<-r_cp[i,1]
    y_r[i]<-qy(r_cp[i,1])
    p_r[i]<-qp(r_cp[i,2]) 
    if (floor(i/B*10)!=floor((i+1)/B*10))
      print(paste((floor(i/B*10)+1)*10,'%',sep=''))
  }
  er<-mean((1+y_r)*(p_r))
  L<-mean(sapply((lambda*er-(1+y_r)*(p_r)),max,0))
  return(L/er/lambda)
}
q_f<-function(f,LB,UB,eps=1e-8,max.iter=100){
  l<-LB
  u<-UB
  count<-0
  f_res<-function(p){
    repeat{
      va<-f(l)-p
      vb<-f(u)-p
      vc<-f((l+u)/2)-p
      if(va*vc<=0){
        u<-(l+u)/2
      }else{
        l<-(l+u)/2
      }
      count<-count+1
      if(abs(va-vb)<eps|(count>=max.iter)){
        return((l+u)/2)
      }
    }
  }
  return(f_res)
}
# Cramer-von Mises Statistics
cal_sn<-function(copula,data,...){
  N_obs<-nrow(data)
  cfit<-fitCopula(copula,data,...)
  cn<-rep(0,N_obs)
  for (i in 1:N_obs){
    cn[i]<-sum((data[,1]<=data[i,1])&(data[,2]<=data[i,2]))
  }
  cn<-cn/N_obs
  ctheta<-pCopula(data,cfit@copula)
  return(sum((cn-ctheta)^2))
}
# pr_bootstrap(res_cpl@copula,q_f(p_temp,LB,UB),function(p) qlnorm(p,mean(log(std_p)),sd(log(std_p))),.9)