KL=read.csv('KL_distance/TAK-931-KI.csv',h=T)
rownames(KL)=KL$X
plot(KL$kl,response[rownames(KL),'RTV_tc'])
cor(KL$kl,response[rownames(KL),'AUC'],method='spearman')
