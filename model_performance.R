getStat=function(m){
  vd_df=m$model$stat #validation statistics 
  t_df=m$prediction #prediction result
  t_df[t_df==1]='Responder'
  t_df[t_df==2]='Non_Responder'
  obs=t_df[1,]
  metrics=function(pred,obs){
    cm=confusionMatrix(pred,obs,positive="Responder")
    ms=c(cm$overall[1:2],cm$byClass[1:2])
    ms
  }
  test_stat=apply(t_df[-1,],1,metrics,t_df[1,])
  rbind(t(vd_df),test_stat)
}
all_stat=lapply(cf_models_TP53m,getStat)
names_stat=rownames(all_stat[[1]])
plot_stats=list()
ind=1
for(name in names_stat){
  plot_stats[[name]]=do.call(rbind,lapply(all_stat,function(x,ind) x[ind,],ind))
  ms=apply(plot_stats[[name]],2,mean)
  print(ms)
  if(ind==1){
    par(mar=c(7.1,5.1,3.1,3.1))
    plot(ms,pch=8,ylab="Average Score",xlab="",col=ind,ylim=c(0,1),main="Average score~Method",xaxt='n',cex.lab=1.5,cex.main=1.5)
    lines(ms,col=ind,lty=2)
    axis(1,at=1:length(ms),labels=names(ms),las=2,font=2)
  }else{
    points(ms,col=ind,pch=8)
    lines(ms,col=ind,lty=2)
  }
  ind=ind+1
  if(ind>3) break
}
legend(locator(1),names_stat[1:3],pch=8,col=1:3,lty=2,bty='n')
ind=4
for(name in names_stat[4:5]){
  plot_stats[[name]]=do.call(rbind,lapply(all_stat,function(x,ind) x[ind,],ind))
  ms=apply(plot_stats[[name]],2,mean)
  print(ms)
  if(ind==4){
    par(mar=c(7.1,5.1,3.1,3.1))
    plot(ms,pch=8,ylab="Average Score",xlab="",col=ind-3,ylim=c(0,1),main="Average score~Method",xaxt='n',cex.lab=1.5,cex.main=1.5)
    lines(ms,col=ind-3,lty=2)
    axis(1,at=1:length(ms),labels=names(ms),las=2,font=2)
  }else{
    points(ms,col=ind-3,pch=8)
    lines(ms,col=ind-3,lty=2)
  }
  ind=ind+1
}
legend(locator(1),names_stat[4:5],pch=8,col=1:2,lty=2,bty='n')
ind=6
for(name in names_stat[6:7]){
  plot_stats[[name]]=do.call(rbind,lapply(all_stat,function(x,ind) x[ind,],ind))
  ms=apply(plot_stats[[name]],2,mean)
  print(ms)
  if(ind==6){
    par(mar=c(7.1,5.1,3.1,3.1))
    plot(ms,pch=8,ylab="Average Score",xlab="",col=ind-5,ylim=c(0,1),main="Average score~Method",xaxt='n',cex.lab=1.5,cex.main=1.5)
    lines(ms,col=ind-5,lty=2)
    axis(1,at=1:length(ms),labels=names(ms),las=2,font=2)
  }else{
    points(ms,col=ind-5,pch=8)
    lines(ms,col=ind-5,lty=2)
  }
  ind=ind+1
  #if(ind>3) break
}
legend(locator(1),names_stat[6:7],pch=8,col=1:2,lty=2,bty='n')
metric=do.call(rbind,lapply(plot_stats,colMeans))
