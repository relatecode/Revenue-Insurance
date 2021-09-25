N_lambda<-2
drcopula<-function(std_val,std_p,N_county,sd_cty){
  LB<-min(std_val)*2
  UB<-max(std_val)*2
  fbase<-logspline(std_val,lbound=LB,ubound=UB)
  pit_seq<-plogspline(std_val,fbase)
  pr_seq<-matrix(0,N_lambda,N_county)
  for (i in 1:N_county){
    p_pit<-plnorm(std_p[,i],mean(log(std_p[,i])),sd(log(std_p[,i])))
    d_tempp<-function(x){dlogspline(x,fbase)*dnorm(x,0,sd_cty[i])}
    c_inte<-integrate(d_tempp,LB,UB)$value
    d_temp<-function(x){d_tempp(x)/c_inte}
    p_temp<-function(x) sapply(x,function(x) integrate(d_temp,LB,x)$value)
    y_pit_i<-p_temp(std_val[seq(i,length(std_val),N_county)])
    res_cpl<-fitCopula(claytonCopula(),cbind(y_pit_i,p_pit),'itau')
    print(y_pit_i)
    print(p_pit)
    for (j in 1:N_lambda){
      pr_seq[j,i]<-premium_rate_yp(res_cpl@copula,p_temp,d_temp,function(x) plnorm(x,mean(log(std_p)),sd(log(std_p))),function(x) dlnorm(x,mean(log(std_p)),sd(log(std_p))),LB,UB,0,max(std_p)*2,1-(j-1)*0.05)
    }
    cat(i,'(dr)  ')
  }
  return(pr_seq)
}
empcopula<-function(std_val,std_p,N_county,sd_cty,B=10000){
  LB<-min(std_val)*2
  UB<-max(std_val)*2
  fbase<-logspline(std_val,lbound=LB,ubound=UB)
  pit_seq<-plogspline(std_val,fbase)
  pr_seq<-matrix(0,N_lambda,N_county)
  for (i in 1:N_county){
    std_yi<-std_val[seq(i,length(std_val),N_county)]
    lengthy<-length(std_yi)
    p_pit<-(rank(std_p[,i])-0.5)/lengthy
    y_pit_i<-(rank(std_yi)-0.5)/lengthy
    res_cpl<-fitCopula(claytonCopula(),cbind(y_pit_i,p_pit),'itau')
    r_cp<-rCopula(B,res_cpl@copula)
    for (j in 1:N_lambda){
      y_r<-rep(0,B)
      p_r<-y_r
      for (k in 1:B){
        y_r[k]<-quantile(std_yi,r_cp[k,1])
        p_r[k]<-quantile(std_p[,i],r_cp[k,2])
      }
      er<-mean((1+y_r)*(p_r))
      lbd<-1-(j-1)*0.05
      L<-mean(sapply((lbd*er-(1+y_r)*(p_r)),max,0))
      pr_seq[j,i]<-L/lbd/er
    }
    cat(i,'(dr)  ')
  }
  return(pr_seq)
}
simu_gene<-function(copula,qy,qp,N_size=15){
  r_cp<-rCopula(N_size,copula)
  y_r<-rep(0,N_size)
  p_r<-y_r
  for (i in 1:N_size){
    temp<<-r_cp[i,1]
    y_r[i]<-qy(r_cp[i,1])
    p_r[i]<-qp(r_cp[i,2])
  }
  return(cbind(y_r,p_r))
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
N_rep<-10
drcpl<-rep(0,N_rep)
norplug<-rep(0,N_rep)
empcpl<-rep(0,N_rep)
for (k in 1){
  cat(k,'\n')
  simu_y<-c()
  simu_p<-c()
  simu_sd_cty<-c()
  pr_seq<-matrix(0,N_lambda,N_county)
  pr_plugin<-pr_seq
  for (i in 1:N_county){
    d_tempp<-function(x){dlogspline(x,fbase)*dnorm(x,0,sd_cty[i])}
    c_inte<-integrate(d_tempp,LB,UB)$value
    d_temp<-function(x){d_tempp(x)/c_inte}
    p_temp<-function(x) sapply(x,function(x) integrate(d_temp,LB,x)$value)
    y_pit_i<-p_temp(std_val[seq(i,length(std_val),N_county)])
    res_cpl<-fitCopula(claytonCopula(),cbind(y_pit_i,p_pit),'itau')
    data_temp<-simu_gene(res_cpl@copula,q_f(p_temp,LB,UB),function(p) qlnorm(p,mean(log(std_p)),sd(log(std_p))))
    simu_y<-cbind(simu_y,data_temp[,1])
    simu_sd_cty<-c(sd_cty,sd(data_temp[,1]))
    simu_p<-cbind(simu_p,data_temp[,2])
    for (j in 1:N_lambda){
      y_r<-data_temp[,1]
      p_r<-data_temp[,2]
      er<-mean((1+y_r)*(p_r))
      lbd<-1-(j-1)*0.05
      L<-mean(sapply((lbd*er-(1+y_r)*(p_r)),max,0))
      pr_plugin[j,i]<-L/er/lbd
      pr_seq[j,i]<-premium_rate_yp(res_cpl@copula,p_temp,d_temp,function(x) plnorm(x,mean(log(std_p)),sd(log(std_p))),function(x) dlnorm(x,mean(log(std_p)),sd(log(std_p))),LB,UB,0,max(std_p)*2,1-(j-1)*0.05)
    }
    cat(i,'  ')
  }
  cat('plugin\n')
  print(pr_plugin)
  simu_y<-c(t(simu_y))
  pr_temp<-drcopula(simu_y,simu_p,N_county,simu_sd_cty)
  pr_empc<-empcopula(simu_y,simu_p,N_county,simu_sd_cty)
  cat('drcpl\n')
  print(pr_temp)
  cat('real\n')
  print(pr_seq)
  drcpl[k]<-mean((pr_seq-pr_temp)^2)
  empcpl[k]<-mean((pr_seq-pr_empc)^2)
  norplug[k]<-mean((pr_seq-pr_plugin)^2)
  cat('\n')
}