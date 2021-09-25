legendre_basis<-function(deg,lower,upper,is.std=TRUE){
  l_func<-function(x){
    x<-(2*x-lower-upper)/(upper-lower)
    res<-0
    for (k in 0:deg)
      res<-res+x^k*choose(deg,k)*choose((deg+k-1)/2,deg)
    return(2^deg*res)
  }
  const_adj<-1
  if (is.std)
    const_adj<-integrate(function(x){l_func(x)^2},lower,upper)$value
  return(function(x){l_func(x)/sqrt(const_adj)})
}
# generate basis function sequence
gene_bsfunc_list<-function(max_deg,bsfunc,...){
  bsfunc_list<-c()
  for (d in 1:max_deg)
    bsfunc_list<-c(bsfunc_list,bsfunc(d,...))
  return(bsfunc_list)
}
# quadratic spline regression
quad_spline_reg<-function(x,y,fct,knot_seq=NULL,is.plot=F){
  if (is.null(knot_seq))
    knot_seq<-quantile(x,c(1/3,2/3))
  D_mat<-cbind(x,x^2)
  for (k in knot_seq)
    D_mat<-cbind(D_mat,sapply(x-k,max,0)^2)
  res<-glm(y~D_mat+fct)
  if (is.plot)
    plot(x,D_mat%*%coef(res)[2:(3+length(knot_seq))])
  return(res)
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
