function(y,flist){
  N_obs<-nrow(y)
  N_county<-ncol(y)
  N_flist<-length(flist)
  for (i in 1:N_county){
    temp<-y[,i]
    fn<-rank(temp)/N_obs
    Xmat<-matrix(0,N_obs,N_flist)
    for (i in 1:N_flist){
      Xmat[,i]<-flist[[i]](temp)
    }
    A<-t(rbind(diag(1,N_flist),diag(-1,N_flist)))
    CA<-cbind(A,matrix(c(rep(1,N_flist),rep(-1,N_flist)),nrow=N_flist))
    B<-c(rep(0,N_flist),rep(-1,N_flist),1,-1)
    solve.QP(t(Xmat)%*%Xmat,t(Xmat)%*%fn,CA,B)$solution
  }
}