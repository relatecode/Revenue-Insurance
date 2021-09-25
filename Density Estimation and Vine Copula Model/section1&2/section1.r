# R code for 'Calculation of Crop Insurance Premium:Based on a Density Ratio Model and Vine Copula¡®
# section 1

dat=read.delim("D:/Y/´ó¶¹/shandong.txt",header= TRUE,sep="\t", stringsAsFactors = FALSE)

# load the following libraries and functions
library('logspline')
library('lme4')
library('fBasics')



# fit a simple model with spline time trend and county dummies
# use data for year 2006 and onward
dat=dat[dat$Year>2005,]


# two point quadratic spline for time trend
t2=quantile(dat$Year,c(1/3,2/3))
year=cbind(dat$Year,dat$Year^2,(dat$Year-t2[1])^2*(dat$Year>t2[1]),(dat$Year-t2[2])^2*(dat$Year>t2[2]))
City=dat$City

# first estimate at the district level
fit=glm(Yield~year+factor(City),data=dat)

# check the time trend
t0=unique(dat$Year)
t0=sort(t0)
tt0=cbind(t0,t0^2,(t0-t2[1])^2*(t0>t2[1]),(t0-t2[2])^2*((t0>t2[2])))
trend=tt0%*%coef(fit)[2:5]
plot(t0,trend,type='b')


summary(resid(fit))
plot(logspline(fit$res))


# studentize residuals by predicted value
summary(fitted(fit))
res=resid(fit)/fitted(fit)
summary(res)
plot(density(res))


# use logspline for the baseline 
# as it is easier to calculate the normalizing factor here
f0=logspline(res,lbound=min(res),ubound=max(res))
y=plogspline(res,f0)


# divide 12 years of data into 6 cells
y2=cut(y,breaks=seq(0,1,length=7))
w=xtabs(~y2+dat$City.code)
w=unlist(w)
w=matrix(w,ncol=1)
cy=unique(dat$City.code)
cy=sort(cy)
cy=rep(cy,each=6)
temp=seq(0,1,length=7)
z=temp[-1]/2+temp[-7]/2


dat2=data.frame(cbind(w,cy,rep(z,17)))
names(dat2)=c('w','cy','z')
dat2$cy=factor(dat2$cy)


# generate Legendre basis 
zz=cbind(legendre(dat2$z,1,2),legendre(dat2$z,2,2),legendre(dat2$z,3,2),legendre(dat2$z,4,2),legendre(dat2$z,5,2),legendre(dat2$z,6,2),legendre(dat2$z,7,2),legendre(dat2$z,8,2))


# model specfication via Poisson regression
# Table 1 in the text is based on these models
# Model 1 (fit1) is preferred accroding to AIC.But!K = 1 is apparently too restrictive ,so we still choose K=2!
fit1=glm(w~zz[,1]:factor(cy),data=dat2,family='poisson')
fit2=glm(w~zz[,1]:factor(cy)+zz[,2]:factor(cy),data=dat2,family='poisson')
fit3=glm(w~zz[,1]:factor(cy)+zz[,2]:factor(cy)+zz[,3]:factor(cy),data=dat2,family='poisson')
fit4=glm(w~zz[,1]:factor(cy)+zz[,2]:factor(cy)+zz[,3]:factor(cy)+zz[,4]:factor(cy),data=dat2,family='poisson')
fit5=glm(w~zz[,1]:factor(cy)+zz[,2]:factor(cy)+zz[,3]:factor(cy)+zz[,4]:factor(cy)+zz[,5]:factor(cy),data=dat2,family='poisson')
fit6=glm(w~zz[,1]:factor(cy)+zz[,2]:factor(cy)+zz[,3]:factor(cy)+zz[,4]:factor(cy)+zz[,5]:factor(cy)+zz[,6]:factor(cy),data=dat2,family='poisson')
c(fit1$aic,fit2$aic,fit3$aic,fit4$aic,fit5$aic,fit6$aic)

# calculate f0
x0=seq(min(res),max(res),length=300)
out0=dlogspline(x0,f0)
plot(x0,out0,type='l')

# function to calculate densities
ff=function(x,b,f.logspline,id)
{
  d0=dlogspline(x,f.logspline)
  p0=plogspline(x,f.logspline)
  tilt=0
  for (i in 1:length(id))
  {
    y=legendre(p0,id[i],2)
    tilt=tilt+b[i]*y
  }  
  d=d0*exp(tilt)
  return(d)
}


# repeat this for each city
a1=matrix(0,ncol=1,nrow=17)
out1=matrix(0,ncol=17,nrow=300)
id=c(1,2)
b1=coef(fit2)[-1]
b1=matrix(b1,ncol=2)
for (i in 1:17)
{
  a1[i]=integrate(ff,lower=min(res),upper=max(res),b=b1[i,],f.logspline=f0,id=id)$value
  temp=ff(x0,b1[i,],f0,id)
  out1[,i]=temp/a1[i]
}
matplot(x0,out1,type='l')
lines(x0,out0,lwd=3)



# for each city, plot the district common density and individually estimated density 
out.city=matrix(0,ncol=17,nrow=length(x0))
for (i in seq(1,17))
{
  temp=res[dat$City.code==i]
  f00=logspline(temp)
  out.city[,i]=dlogspline(x0,f00)
}


matplot(x0,out.city,type='l')
lines(x0,out0,lwd=3)


## semi-parametric method to obtain the PIT of yield
YU=matrix(0,ncol=17,nrow=12)
for (i in 1:17)
  for(j in 1:12){
    YU[j,i]=integrate(ff,lower = min(res),upper = res[ j +12*(i-1)],b=b1[i,],f.logspline=f0,id=id)$value/a1[i]
  }

