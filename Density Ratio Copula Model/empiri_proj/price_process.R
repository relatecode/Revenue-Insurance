# dtst<-function(x){dnorm(x)}
# ptst<-function(q){pnorm(q)}
# qtst<-function(p){qnorm(p)}
# aa<-mvrnorm(10,c(0,0),diag(1,2))
# tstMvd <- mvdc(copula = normalCopula(),margins = c("tst", "tst"), paramMargins = list(list(),list()))
# fitmvdc<-fitMvdc(aa,tstMvd, start=c(sin(cor(mat_adj[,3],mat_adj[,4], method = "kendall")*pi/2)),optim.control = list(trace = F, maxit = 2000))
require('mixR')
require('Peacock.test')
p_tb<-read.csv('yp_data.csv')
plot(p_tb$Year,p_tb$Value)
cpi_tb<-read.csv('CPIAUCSL.csv')
id_s2<-which(cpi_tb$DATE=='1950-01-01')
id_e2<-which(cpi_tb$DATE=='2010-12-01')
p_raw<-p_tb$Value
p_adj<-p_raw
p_adj<-p_raw/cpi_tb[id_s2:id_e2,2]*cpi_tb[id_e2,2]

res<-quad_spline_reg(p_tb$Year,p_adj,p_tb$Period)
p_nor<-resid(res)/fitted(res)
plot(p_nor)
mix_res<-mixfit(p_nor+1,3,'lnorm')
pmixlnorm<-function(x,mix_res){
  res<-0
  for (i in 1:length(mix_res$pi))
    res<-res+plnorm(x,meanlog=mix_res$mulog[i],sdlog=mix_res$sdlog[i])*mix_res$pi[i]
  return(res)
}
dmixlnorm<-function(x,mix_res){
  res<-0
  for (i in 1:length(mix_res$pi))
    res<-res+dlnorm(x,meanlog=mix_res$mulog[i],sdlog=mix_res$sdlog[i])*mix_res$pi[i]
  return(res)
}
p_pit<-pmixlnorm(p_nor[p_tb$Period=='FEB']+1,mix_res)

