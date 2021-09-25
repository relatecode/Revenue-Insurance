require(fitdistrplus)
require(copula)
pr_cal<-function(yield_pit,price,lambda){
  pfseq<-c(plnorm,pweibull)
  N_pmod<-length(pfseq)
  qfseq<-c(qlnorm,qweibull)
  res<-rep(0,N_pmod)
  for (i in 1:N_pmod)
    res[i]<-fitdist(price,pfseq[i])$aic
  id_p<-which(res==min(res))[1]
  res<-fitdist(price,pfseq[i])
  price_pit<-pfseq[[id_p]](price,res$estimate)
  temp<-c(yield_pit,price_pit)
  sn_seq[1]<-cal_sn(normalCopula(),temp)
  sn_seq[2]<-cal_sn(frankCopula(),temp)
  sn_seq[3]<-cal_sn(claytonCopula(),temp,'itau')
  return(pr_bootstrap(copula,qy,qp,lambda))
}