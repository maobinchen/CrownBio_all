#!/public/software/bin/Rscript
load('lmm.Rdata')
dat1$CT=substr(dat1$PDX,1,2)
dat1=subset(dat1,PDX %in% colnames(log_filteredExpr3))
library(nlme)

run=function(ichunk,dat,ex,quadratic=T){
  require(nlme)
  m=length(ichunk)
  fit_table=data.frame(Intercept=numeric(m),Day=numeric(m),
                       Day_Treatment=numeric(m),Day_Gene=numeric(m),Day_Treatment_Gene=numeric(m))
  if(quadratic){
    fit_table=data.frame(Intercept=numeric(m),Day=numeric(m),DaySq=numeric(m),
                         Day_Treatment=numeric(m),Day_Gene=numeric(m),DaySq_Treatment=numeric(m),DaySq_Gene=numeric(m),
                         Day_Treatment_Gene=numeric(m),DaySq_Treatment_Gene=numeric(m))
  }		
  if('TV0' %in% names(dat)){
    if(quadratic){
      fit_table=data.frame(Intercept=numeric(m),tv0=numeric(m),Day=numeric(m),DaySq=numeric(m),
                           Day_Treatment=numeric(m),Day_Gene=numeric(m),DaySq_Treatment=numeric(m),DaySq_Gene=numeric(m),
                           Day_Treatment_Gene=numeric(m),DaySq_Treatment_Gene=numeric(m))
    }else{
      fit_table=data.frame(Intercept=numeric(m),tv0=numeric(m),Day=numeric(m),
                           Day_Treatment=numeric(m),Day_Gene=numeric(m),Day_Treatment_Gene=numeric(m))
    }
  }	
  ex_df=data.frame(sigma=numeric(m),MAE=numeric(m),MSE=numeric(m),AIC=numeric(m),BIC=numeric(m),logLik=numeric(m),p_val_tr=numeric(m),p_val_gtr=numeric(m))
  fit_table = cbind(fit_table,ex_df)
  j=0
  for(gene in ichunk){
    j = j+1
    gene_ex=data.frame(PDX=colnames(ex),gene=unlist(ex[gene,]))
    merged=merge(dat,gene_ex,by='PDX')
    f=paste("logTV ~ Day + Day:Treatment + Day:gene + Day:Treatment:gene",sep = "")		
    random_f=list(PDX = ~Day+1, Mouse = ~Day+1)
    if('TV0' %in% names(dat)){
      f="logTV ~ TV0+Day + Day:Treatment + Day:gene + Day:Treatment:gene"
    }
    if(quadratic){
      quadratic_term = "+ Day_sq + Day_sq:Treatment + Day_sq:gene + Day_sq:Treatment:gene"
      f=paste(f,quadratic_term,sep="")
      random_f = list(PDX = ~Day_sq+Day-1, Mouse = ~Day_sq+Day-1)	#to speed up
    }
    model3a.fit=NULL
    tryCatch(model3a.fit<-lme(as.formula(f),random = random_f,data=merged, method='REML',control = lmeControl(opt = "optim")),error=function(e) {print(paste("Error:",e))} )
    if(is.null(model3a.fit)) next
    sum_m=summary(model3a.fit)	
    p_values=sum_m$tTable[,5]
    dif=2
    if(quadratic) dif=3
    p_value = p_values[c(length(p_values)-dif,length(p_values))] #only save the p-value of three way interaction term
    residuals = model3a.fit$residuals[,3]
    sum_stat = c(model3a.fit$sigma,mean(abs(residuals)),sqrt(sum(residuals^2)/length(residuals)),sum_m$AIC,sum_m$BIC,sum_m$logLik,p_value)
    len_terms = length(t(sum_m$tTable[,1]))
    fit_table[j,1:len_terms]=t(sum_m$tTable[,1])
    fit_table[j,(len_terms+1):(len_terms+8)]=sum_stat	
  }
  row.names(fit_table)=ichunk
  fit_table
}
paraRun = function(cls,fun,b_vec,dat,ex,quadratic=T){
  n=length(b_vec)
  nc=length(cls)
  options(warn=-1)
  ichunks=split(b_vec,1:nc)
  options(warn=0)
  result=clusterApply(cls,ichunks,fun,dat,ex,quadratic)
}

ss_fit = function(dat,ex,nrep=10,pct=0.5){
  require(snow)
  nc=8
  b_vec = rownames(ex)
  pdxs = colnames(ex)
  fit_tables = list()
  for(i in 1:nrep){
    sample_pdxs = sample(pdxs,length(pdxs)*pct)
    dat_s = subset(dat,PDX %in% sample_pdxs)
    ex_s =ex[,sample_pdxs]
    cl=makeCluster(type="SOCK",rep("localhost",nc))
    pp=paraRun(cl,run,b_vec,dat_s,ex_s,quadratic=F)
    fit_tables[[i]]=do.call("rbind",pp)
    stopCluster(cl)
  }
  fit_tables
}

fit_tables=ss_fit(dat1,log_filteredExpr3[1:10,],nrep=1,pct=1)
save(fit_tables,file=paste(Sys.Date(),"fit_tables.RData",sep="_")) #save the LMM fit output 

post_process_fit_table = function(fit_table,ex_df,kw){
  fit_table=fit_table[rownames(ex_df),]
  mean_ex=apply(ex_df,1,mean)
  min_ex=apply(ex_df,1,min)
  max_ex=apply(ex_df,1,max)
  sd_ex=apply(ex_df,1,sd)
  fit_table = cbind(fit_table,'mean_ex'=mean_ex,'sd_ex'=sd_ex,'min_ex'=min_ex,'max_ex'=max_ex)
  write.csv(fit_table,paste(kw,"_fit.csv",sep=""))
  fit_table
}

#fit the lmm multiple times, see how robust those top 100 genes are
#lmm_genes_lst=lapply(fit_tables_train,function(x) rownames(x[order(x$AIC),][1:100,]))
#table(unlist(lmm_genes_lst))
load('2017-12-04_fit_tables.RData') 
fit_table=fit_tables[[1]]
fit_table$pgg_d_pd=fit_table$p_val_g/fit_table$p_val_day #gene effect on tumor growth
fit_table$pgt_d_pt=fit_table$p_val_gtr/fit_table$p_val_tr #treatment effect on gene 
fit_table=post_process_fit_table(fit_table,log_filteredExpr3,'lmm_1st') #top 100 genes enriched in Cell Cycle, even more enriched for Coef(gene*Day)>0
sel_genes=rownames(subset(fit_table,pgt_d_pt<0.01 & Day_Gene>0))
sel_genes=rownames(subset(fit_table,pgt_d_pt<0.01 & pgg_d_pd<0.1))

sel_genes_100=rownames(fit_table[order(fit_table$pgt_d_pt),])[1:100]
plot(log10(fit_table$pgt_d_pt),fit_table$mean_ex)
boxplot(fit_table$mean_ex~cut_number(fit_table$pgt_d_pt,10))
fit_table_cor=cbind(fit_table,t(cor_all)[rownames(fit_table),])
fit_table_cor$class=round(log10(fit_table$pgt_d_pt))
plot(log10(fit_table$pgt_d_pt),fit_table_cor$AUC)
boxplot(fit_table_cor$TGI~cut_number(fit_table_cor$pgt_d_pt,20))
boxplot(fit_table_cor$AUC~fit_table_cor$class)
#write.csv(fit_table_cor,paste('lmm_1st',"_fit_cor.csv",sep=""))
gs2Reactome=function(gs,kw='unknown',pval=0.01){
  require(ReactomePA)
  require("biomaRt")
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  gsInfo=getBM(attributes=c('hgnc_symbol','entrezgene','gene_biotype'),filters = 'hgnc_symbol', values = gs,mart=ensembl)
  eids=gsInfo[[2]]
  names(eids)=gsInfo[[1]]
  eids=eids[gs]
  require(ReactomePA)
  pwys <- enrichPathway(gene=eids,pvalueCutoff=pval, readable=T)
  print(barplot(pwys,showCategory = 15))
  dev.copy(png,paste0(kw,'_reactome_barplot.png'),width=1600,height=800, res=150)
  dev.off()
  print(dotplot(pwys))
  dev.copy(png,paste0(kw,'_reactome_dotplot.png'),width=1200,height=800, res=150)
  dev.off()
  print(enrichMap(pwys, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1))
  dev.copy(png,paste0(kw,'_enrichmap.png'),width=1600,height=1600, res=150)
  dev.off()
  write.csv(pwys@result,paste0(kw,"_reactome_pathways.csv"),row.names = F)
}  
gs2Reactome(sel_genes,'Pval_ratio_sel')
require(fgsea)
cp_pathways = gmtPathways(gmt.file = "c2.cp.v6.1.symbols.gmt")
pwy='REACTOME_CELL_CYCLE'
pwy_genes=intersect(cp_pathways[[pwy]],sel_genes)
apply(log_ex3[,pwy_genes],1,sum)
res_df=cbind('Responder'=response[rownames(log_ex3),'Responder'],log_ex3[,sel_genes])
res_df=subset(res_df,Responder!=2)
res_df$Responder=factor(ifelse(res_df$Responder==1,'Responder','Non_Responder'),levels=c('Responder','Non_Responder'))
cf_genes=feature_sel(res_df,min_f=10,iter=200,subsample = 1,nrep = 1)
cf_sel=res_df[,c('Responder',cf_genes$sel_genes)]
featurePlot(x =cf_sel[,-1],y =cf_sel[,1], between = list(x = 1, y = 1),type = c("g", "p"))

##model building
cf_genes_lmm=cbind('Responder'=response[rownames(log_ex3),'Responder'],log_ex3[,sel_genes_100])
cf_genes_lmm=subset(cf_genes_lmm,Responder!=2)
cf_genes_lmm$Responder=factor(ifelse(cf_genes_lmm$Responder==1,'Responder','Non_Responder'),levels=c('Responder','Non_Responder'))
cf_models_genes_lmm=trainModel(cf_genes_lmm,nrep=1,fs=TRUE,ROC=TRUE)


spl=createDataPartition(res_df$Responder,times=10,p=0.7)
save(spl,file='spl.Rdata')
train_pdxs_lst=lapply(spl,function(x) rownames(res_df[x,]))
save(train_pdxs_lst,file='train_pdxs_lst.Rdata')

load('2017-12-05_fit_tables.RData')
nrep=10
result=list()
for (i in 1:nrep){
  sel_pdxs=train_pdxs_lst[[i]]
  test_pdxs=rownames(res_df[-spl[[i]],])
  fit_table=fit_tables[[i]]
  fit_table$pgt_d_pt=fit_table$p_val_gtr/fit_table$p_val_tr
  sel_genes=rownames(fit_table[order(fit_table$pgt_d_pt),])[1:100]
  #sel_genes=rownames(subset(fit_table,pgt_d_pt<0.01))
  cf_genes_lmm=cbind('Responder'=response[rownames(log_ex3),'Responder'],log_ex3[,sel_genes])
  cf_genes_lmm=subset(cf_genes_lmm,Responder!=2)
  cf_genes_lmm$Responder=factor(ifelse(cf_genes_lmm$Responder==1,'Responder','Non_Responder'),levels=c('Responder','Non_Responder'))
  train_df=cf_genes_lmm[sel_pdxs,]
  test_df=cf_genes_lmm[test_pdxs,]
  filtered_genes=feature_sel(train_df,min_f=5,iter=200,nrep=10,m='spearman',subsample = 0.8,consistency = 0.5)
  train_sel=train_df[,c(names(train_df)[1],filtered_genes$sel_genes)]
  test_sel=test_df[,c(names(test_df)[1],filtered_genes$sel_genes)]
  print(featurePlot(x = train_sel[,-1],y = train_sel[,1], between = list(x = 1, y = 1),type = c("g", "p")))
  models=buildModels(train_sel,ROC=TRUE,all=FALSE)
  predicted=rbind(observed=test_sel[,1],do.call(rbind,lapply(lapply(models$models,predict,test_sel[,-1]),as.vector)))
  result[[i]]=list('model'=models,'prediction'=predicted,'test'=test_sel)
  print(featurePlot(x = test_sel[,-1],y = test_sel[,1], between = list(x = 1, y = 1),type = c("g", "p")))
  result
}  
lapply(result,getStat)
#cf_models_lmm_genes=result
#save(cf_models_lmm_genes,file='cf_models_lmm_genes.Rdata')


###load lmm pathway data
load('lmm_pwys_fit_tables.RData')
fit_table=fit_tables[[1]]
fit_table$pgg_d_pd=fit_table$p_val_g/fit_table$p_val_day #gene effect on tumor growth
fit_table$pgt_d_pt=fit_table$p_val_gtr/fit_table$p_val_tr #treatment effect on gene 
fit_table=post_process_fit_table(fit_table,pwys_df,'lmm_1st_pwys')
colnames(cor_all_pwy)=make.names(colnames(cor_all_pwy))
fit_table_cor=cbind(fit_table,t(cor_all_pwy)[rownames(fit_table),])

write.csv(fit_table_cor,'lmm_1st_pwys_fit_cor.csv')
sel_pwys=rownames(fit_table[order(fit_table$p_val_gtr),][1:20,])

sel_pwys_df=res_pwys_df[,c('Responder',sel_pwys)]
require(corrplot)
corrplot(cor(sel_pwys_df[,-1],use='pairwise.complete.obs',method='spearman'), 
         order="hclust",hclust.method="ward.D2",tl.col="black", tl.srt=45,addrect = 4, rect.col = "red")
res_ex_sel=cbind(response[rownames(log_ex3),],log_ex3[,tgi_genes$sel_genes])
ggplot(sel_pwys_df)+aes(factor(Responder),REACTOME_CHROMOSOME_MAINTENANCE)+geom_boxplot(notch=T,varwidth = T)+
  geom_jitter(width=0.1)+facet_wrap(~substr(rownames(sel_pwys_df),1,2))+mytheme
ggplot(sel_pwys_df)+aes(factor(Responder),REACTOME_TELOMERE_MAINTENANCE)+geom_boxplot(notch=T,varwidth = T)+
  geom_jitter(width=0.1)+facet_wrap(~substr(rownames(sel_pwys_df),1,2))+mytheme
ggplot(sel_pwys_df)+aes(factor(Responder),REACTOME_CELL_CYCLE)+geom_boxplot(notch=T,varwidth = T)+
  geom_jitter(width=0.1)+facet_wrap(~substr(rownames(sel_pwys_df),1,2))+mytheme
ggplot(sel_pwys_df)+aes(factor(Responder),BIOCARTA_P27_PATHWAY)+geom_boxplot(notch=T,varwidth = T)+
  geom_jitter(width=0.1)+facet_wrap(~substr(rownames(sel_pwys_df),1,2))+mytheme
ggplot(sel_pwys_df)+aes(factor(Responder),REACTOME_DEPOSITION_OF_NEW_CENPA_CONTAINING_NUCLEOSOMES_AT_THE_CENTROMERE)+geom_boxplot(notch=T,varwidth = T)+
  geom_jitter(width=0.1)+facet_wrap(~substr(rownames(sel_pwys_df),1,2))+mytheme+theme(axis.title.y=element_text(face="bold.italic",size=10, color="brown"))
ggplot(sel_pwys_df)+aes(factor(Responder),PID...fibrinolysis.pathway)+geom_boxplot(notch=T,varwidth = T)+
  geom_jitter(width=0.1)+facet_wrap(~substr(rownames(sel_pwys_df),1,2))+mytheme
cf_models_lmm_pwys_fsr=trainModel(sel_pwys_df,nrep=1,fs=TRUE,ROC=TRUE)


#######################################
load('lmm_cnv_fit_tables.RData')
fit_table=fit_tables[[1]]
fit_table$pgg_d_pd=fit_table$p_val_g/fit_table$p_val_day #gene effect on tumor growth
fit_table$pgt_d_pt=fit_table$p_val_gtr/fit_table$p_val_tr #treatment effect on gene 
fit_table=post_process_fit_table(fit_table,cnv_genes_info,'lmm_1st_cnv')
fit_table_cor=cbind(fit_table,t(cor_all_cnv)[rownames(fit_table),])

write.csv(fit_table_cor,'lmm_1st_cnv_fit_cor.csv')

###all model, 70% sampling rate, LMM for 10 times
load('2017-12-01_fit_tables.RData')
xx=do.call(rbind,lapply(fit_tables,function(x) x$p_val_gtr/x$p_val_tr))
colnames(xx)=row.names(fit_tables[[1]])
sigc=apply(xx,2,function(v) sum(v<0.05))
print(names(sigc[sigc>=8]),quote=F)
xx_r=apply(xx,1,rank)
top100c=apply(xx_r,1,function(v) sum(v<=100))

##

lmm_genes_lst=lapply(fit_tables,function(x) rownames(x[order(x$AIC),][1:200,]))
AIC_genes=rev(sort(table(unlist(lmm_genes_lst))))
require(fgsea)
t_val = fit_tables[[1]]$AIC
names(t_val)=rownames(fit_tables[[1]])
t_val=t_val[order(t_val)] 
cp_pathways = gmtPathways(gmt.file = "c2.cp.v6.1.symbols.gmt")
pathways=cp_pathways
fgseaRes <- fgsea(pathways, t_val, minSize=5, maxSize=500, nperm=1000)
#fgseaRes <- fgsea(pathways, t_val, minSize=10, maxSize=30, nperm=1000)
head(fgseaRes[order(pval), ],n=20)
sum(fgseaRes[, padj < 0.05])
fgseaRes=subset(fgseaRes,padj<0.05)
topPathwaysUp <- fgseaRes[ES > 0][head(rev(order(NES)), n=50), pathway]
