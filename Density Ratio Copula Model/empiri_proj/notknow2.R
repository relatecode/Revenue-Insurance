res<-c()
avgApr<-rep(0,12)
avgOct<-rep(0,12)
hd<-names(dat_price<-read.csv('soybean_price/d1-2006.csv'))
for (i in 2006:2017){
  dat_price<-read.csv(paste('soybean_price/d1-',i,'.csv',sep=''))
  names(dat_price)<-hd
  id_sel1<-dat_price$合约==paste('a',substr(paste(i+1),3,4),'01',sep='')&((dat_price$日期>=i*10000+0401&dat_price$日期<=i*10000+0430))
  id_sel2<-dat_price$合约==paste('a',substr(paste(i+1),3,4),'01',sep='')&((dat_price$日期>=i*10000+1001&dat_price$日期<=i*10000+1031))
  avgApr[i-2005]<-mean(dat_price$收盘价[id_sel1])
  avgOct[i-2005]<-mean(dat_price$收盘价[id_sel2])
  res<-rbind(res,dat_price[id_sel1|id_sel2,])
  print(i)
}
price_seq<-apply(cbind(avgApr,avgOct),1,max)
temp<-price_seq/mean(price_seq)
p_pit<-plnorm(temp,mean(log(temp)),sd(log(temp)))
