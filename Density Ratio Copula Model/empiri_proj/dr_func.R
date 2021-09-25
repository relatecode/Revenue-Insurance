# detrend and standardize
yield_preprocess<-function(yield,year,county,knot_seq=NULL){
  fit_qsp<-quad_spline_reg(year,yield,county,knot_seq)
  res<-resid(fit_qsp)/fitted(fit_qsp)
  return(res)
}
# obtain frequency for possion regression
val2freq<-function(pit_seq,county,N_cell=30){
  cell_seq<-cut(pit_seq,breaks=seq(0,1,length=N_cell+1))
  freq_seq<-matrix(unlist(xtabs(~cell_seq+county)),ncol=1)
  cy=unique(county)
  cy=sort(cy)
  cy=rep(cy,each=N_cell)
  freq_mat<-data.frame(freq=freq_seq,county=factor(cy))
}
possion_reg<-function(freq_mat,N_cell,bsfunc_list,N_county=length(unique(freq_mat$county))){
  temp<-seq(0,1,length=N_cell+1)
  temp<-rep(temp[-1]/2+temp[-(N_cell+1)]/2,N_county)
  l_bs_mat<-c()
  max_deg<-length(bsfunc_list)
  for (i in 1:max_deg){
    l_bs_mat<-cbind(l_bs_mat,bsfunc_list[[i]](temp))
  }
  base_str<-'freq~'
  aic_seq<-c()
  fit_sel<-NULL
  min_aic<-Inf
  for (b in 1:max_deg){
    if (b!=1)
      base_str<-paste(base_str,'+',sep='')
    base_str<-paste(base_str,'l_bs_mat[,',b,']:factor(county)',sep='')
    fit_temp<-glm(eval(parse(text=base_str)),data=freq_mat,family='poisson')
    aic_seq<-c(aic_seq,AIC(fit_temp))
    if(AIC(fit_temp)<min_aic){
      fit_sel<-fit_temp
      min_aic<-AIC(fit_temp)
    }
  }
  id<-which(aic_seq==min(aic_seq))
  b<-matrix(coef(fit_sel)[-1],ncol=id)
  return(list(K=id,fit_sel=fit_sel,aic_seq=aic_seq,b=b))
}
# individual's density
d_indi<-function(b,dbase,pbase,bsfunc_list,lower,upper){
  d_temp<-function(x){
    tilt=0
    for (i in 1:length(b))
    {
      temp<-bsfunc_list[[i]](pbase(x))
      tilt<-tilt+b[i]*temp
    }  
    d<-dbase(x)*exp(tilt)
    return(d)
  }
  c_inte<-integrate(d_temp,lower,upper)$value
  return(function(x){d_temp(x)/c_inte})
}
# premium rate with yield and price
premium_rate_yp<-function(copula,py,dy,pp,dp,lowery,uppery,lowerp,upperp,lambda){
  # djoint<-function(y,p) {
  #   N<-ncol(y)
  #   y<-c(y)
  #   p<-c(p)
  #   matrix(dCopula(cbind(py(y),pp(p)),copula)*dy(y)*dp(p),ncol=N)
  # }
  plot(seq(lowery,uppery,length=1000),dy(seq(lowery,uppery,length=1000)))
  plot(seq(lowerp,upperp,length=1000),dp(seq(lowerp,upperp,length=1000)))
  ejoint<-function(arg){
    y<-arg[1]
    p<-arg[2]
    # RECONSIDER!!!
    (1+y)*(p)*dCopula(cbind(py(y),pp(p)),copula)*dy(y)*dp(p)
  }
  intres<-cubintegrate(ejoint,c(lowery,lowerp),c(uppery,upperp),method='cuhre')
  er<-intres$integral
  # print(intres)
  prjoint<-function(arg){
    y<-arg[1]
    p<-arg[2]
    # RECONSIDER!!!
    max((lambda*er-(1+y)*(p)),0)*dCopula(cbind(py(y),pp(p)),copula)*dy(y)*dp(p)
  }
  xxx<-seq(lowery,uppery,length=100)
  yyy<-seq(lowerp,upperp,length=100)
  mattt<-matrix(0,100,100)
  for (i in 1:100)
    for(j in 1:100)
      mattt[i,j]<-prjoint(c(xxx[i],yyy[j]))
  contour(mattt)
  intres<-cubintegrate(prjoint,c(lowery,lowerp),c(uppery,upperp),method='cuhre')
  L<-intres$integral
  cat('L:',L,'\n')
  cat('er:',er,'\n')
  return(L/er/lambda)
}
# premium rate with only yield
premium_rate_y<-function(dy,lower,upper,lambda){
  # er<-integrate(function(x)(1+x)*dy(x),-1,upper)$value
  L<-integrate(function(x) sapply(lambda-x-1,max,0)*dy(x),lower,upper)$value
  return(L/lambda)
}
