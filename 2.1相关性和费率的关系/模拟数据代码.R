rm(list=ls())

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
  er<-mean((y_r)*(p_r))   
  L<-mean(sapply((lambda*er-(y_r)*(p_r)),max,0))
  return(L/er/lambda)
}
N<-30
pr_seq<-matrix(0,7,N)
for (i in 1:N){
  p<-2
  pho<--0.9
  diagn<-1
  n<-30
  M=matrix(rep(1:p,p),ncol=p,byrow=F)
  corM=pho^(abs(M-t(M)))
  diag(corM)<-diagn
  Z<-mvrnorm(n,rep(6,p),corM)
  yield<-Z[,1]
  price<-Z[,2]
  #plot(1:n,yield,ylim=c(0,9),type = "l")
  #lines(1:n,price)
  
  
  p_pit<-pnorm(price,mean(price),sd(price)) 
  y_pit_i<-pnorm(yield,mean(yield),sd(yield))
  temp<-cbind(y_pit_i,p_pit)
  
  res_cpl<-fitCopula(normalCopula(),temp)
  for (j in 1:7)
    pr_seq[j,i]<-pr_bootstrap(res_cpl@copula,function(q_f1) qnorm(q_f1,mean(yield),sd(yield)),function(p) qnorm(p,mean(price),sd(price)),1-(j-1)*0.05)
  cat(i,'  ')
}


rm(corM,M,res_cpl,Z,diagn,i,j,n,N,p,p_pit,pho,price,temp,y_pit_i,yield)

N<-30
pr_seq_m0.7<-matrix(0,7,N)
for (i in 1:N){
  p<-2
  pho<--0.7
  diagn<-1
  n<-30
  M=matrix(rep(1:p,p),ncol=p,byrow=F)
  corM=pho^(abs(M-t(M)))
  diag(corM)<-diagn
  Z<-mvrnorm(n,rep(6,p),corM)
  yield<-Z[,1]
  price<-Z[,2]
  #plot(1:n,yield,ylim=c(0,9),type = "l")
  #lines(1:n,price)
  
  
  p_pit<-pnorm(price,mean(price),sd(price)) 
  y_pit_i<-pnorm(yield,mean(yield),sd(yield))
  temp<-cbind(y_pit_i,p_pit)
  
  res_cpl<-fitCopula(normalCopula(),temp)
  for (j in 1:7)
    pr_seq_m0.7[j,i]<-pr_bootstrap(res_cpl@copula,function(q_f1) qnorm(q_f1,mean(yield),sd(yield)),function(p) qnorm(p,mean(price),sd(price)),1-(j-1)*0.05)
  cat(i,'  ')
}
rm(corM,M,res_cpl,Z,diagn,i,j,n,N,p,p_pit,pho,price,temp,y_pit_i,yield)

N<-30
pr_seq_0.7<-matrix(0,7,N)
for (i in 1:N){
  p<-2
  pho<-0.7
  diagn<-1
  n<-30
  M=matrix(rep(1:p,p),ncol=p,byrow=F)
  corM=pho^(abs(M-t(M)))
  diag(corM)<-diagn
  Z<-mvrnorm(n,rep(6,p),corM)
  yield<-Z[,1]
  price<-Z[,2]
  #plot(1:n,yield,ylim=c(0,9),type = "l")
  #lines(1:n,price)
  
  
  p_pit<-pnorm(price,mean(price),sd(price)) 
  y_pit_i<-pnorm(yield,mean(yield),sd(yield))
  temp<-cbind(y_pit_i,p_pit)
  
  res_cpl<-fitCopula(normalCopula(),temp)
  for (j in 1:7)
    pr_seq_0.7[j,i]<-pr_bootstrap(res_cpl@copula,function(q_f1) qnorm(q_f1,mean(yield),sd(yield)),function(p) qnorm(p,mean(price),sd(price)),1-(j-1)*0.05)
  cat(i,'  ')
}
rm(corM,M,res_cpl,Z,diagn,i,j,n,N,p,p_pit,pho,price,temp,y_pit_i,yield)

N<-30
pr_seq_0.9<-matrix(0,7,N)
for (i in 1:N){
  p<-2
  pho<-0.9
  diagn<-1
  n<-30
  M=matrix(rep(1:p,p),ncol=p,byrow=F)
  corM=pho^(abs(M-t(M)))
  diag(corM)<-diagn
  Z<-mvrnorm(n,rep(6,p),corM)
  yield<-Z[,1]
  price<-Z[,2]
  #plot(1:n,yield,ylim=c(0,9),type = "l")
  #lines(1:n,price)
  
  
  p_pit<-pnorm(price,mean(price),sd(price)) 
  y_pit_i<-pnorm(yield,mean(yield),sd(yield))
  temp<-cbind(y_pit_i,p_pit)
  
  res_cpl<-fitCopula(normalCopula(),temp)
  for (j in 1:7)
    pr_seq_0.9[j,i]<-pr_bootstrap(res_cpl@copula,function(q_f1) qnorm(q_f1,mean(yield),sd(yield)),function(p) qnorm(p,mean(price),sd(price)),1-(j-1)*0.05)
  cat(i,'  ')
}
#rm(corM,M,res_cpl,Z,diagn,i,j,n,N,p,p_pit,pho,price,temp,y_pit_i,yield)

boxplot(pr_seq[1,],pr_seq_m0.7[1,],pr_seq_0.7[1,],pr_seq_0.9[1,],names = c("pho=-0.9","pho=-0.7","pho=0.7","pho=0.9"),main="lambda=100%")
boxplot(pr_seq[3,],pr_seq_m0.7[3,],pr_seq_0.7[3,],pr_seq_0.9[3,],names = c("pho=-0.9","pho=-0.7","pho=0.7","pho=0.9"),main="lambda=90%")
boxplot(pr_seq[5,],pr_seq_m0.7[5,],pr_seq_0.7[5,],pr_seq_0.9[5,],names = c("pho=-0.9","pho=-0.7","pho=0.7","pho=0.9"),main="lambda=80%")



#一维下费率的计算
a<-yield*price 
#a<-a[(length(a)-24):length(a)]

con_exp<-function(y){
  df<-approxfun(density(a))
  y*df(y)
}

densit<-function(y){
  df<-approxfun(density(a))
  df(y)
}

m<-10
num<-9
pr<-c( )
expect<-quantile(a,5/m)
for (i in 1:7){
  lambda<-1-(i-1)*0.05
  Prob<-integrate(densit, lower =min(a)-0.4, upper =lambda*expect)$value
  Con_exp<-integrate(con_exp, lower =min(a)-0.4, upper =lambda*expect)$value
  pr[i]<-Prob*(lambda*expect-Con_exp/Prob)/(lambda*expect)
}



