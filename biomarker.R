library(mct)
library(tidyverse)
source('./get_AUC_v2.R')
model_to_remove=c('OV0243','OV0273','CR0047')
dat=read_trt_vehicle("TAK931.csv","Vehicle.csv") #LU2512 (two passages used in this study)
mouse_nd=table(dat$Mouse)
mouse1=names(mouse_nd[mouse_nd<4])#at least four observations
dat1=subset(dat,!Mouse %in% mouse1 & !PDX %in% model_to_remove)
AUC=get_AUC(dat1,CI=F)
auc=mct::get_AUC(dat1,CI=F) #old AUC scores
RTV_tc=AUC$AUC_r#relative tumor volume at day 21 (treatment/control)
RTV_tr=AUC$RTV #relative tumor volume t21/t0 for treatment 
wtmc=AUC$wtmc #weighted kt minus kc(weight=2)
TGI=get_TGI(dat1) 
#cor(TGI,RTV_tc)
PDXs=unique(dat1$PDX)
cor(TGI[PDXs],AUC$AUC_r[PDXs],method='spearman')
#code: 1->CR/PR, 2->stable disease(SD) 3,4-> Non-responder
pdx_eff=data.frame('PDX'=PDXs,'TGI'=TGI[PDXs],'RTV_tc'=RTV_tc[PDXs],'RTV_tr'=RTV_tr[PDXs],'wtmc'=wtmc[PDXs],'AUC'=auc$AUC_r[PDXs])
pdx_eff$response='PD'
pdx_eff[pdx_eff$RTV_tr<=1.73,'response']='SD' #stable disease(SD)
pdx_eff[pdx_eff$RTV_tr<=1,'response']='PR' #tumor declining, partial response (CR/PR)
pdx_eff[pdx_eff$RTV_tr<=0.3,'response']='CR' #tumor declining, partial response (CR/PR)
pdx_eff[pdx_eff$RTV_tr>3 & pdx_eff$TGI < 0.5,'response']='CPD' #complete progress disease(PD)
write.csv(pdx_eff,"response_TAK931_auto.csv",row.names=F) #manually check responder/non-responder(1,0), code 2 indicates ambiguous case

sel_PDXs=PDXs
for(pdx in sel_PDXs){
  dat_p=subset(dat,PDX==pdx)
  title=paste0(pdx,':','RTV_tc=',round(RTV_tc[pdx],2),',RTV_tr=',round(RTV_tr[pdx],2),',TGI=',round(TGI[pdx],2))
  with(data=dat_p,print(xyplot(logTV~Day|Treatment,groups=Mouse,type=c('p','l'),main=title)))
  invisible(readline(prompt="Press [enter] to continue"))
}

####gene expression 
response=read.csv('response_TAK931.csv',h=T)
response=subset(response,!PDX %in% model_to_remove)
response$CT=substr(response$PDX,1,2)
response$response=factor(response$response,levels=c('CR','PR','SD','PD','CPD'))
rownames(response)=response$PDX
table(response$Responder,response$response)
table(response$CT,response$response)
ex <- read.csv("C:\\Users\\binchen.mao\\Documents\\crownbio\\HuSignature\\RNAseq_expr_RSEM_estimated_count.csv", header=TRUE)
RNASeq_pdxs=sub("P[0-9]+$","",colnames(ex))
length(unique(RNASeq_pdxs)) #1540 PDX models with RNASeq data
hu_prime_models=read.csv("C:\\Users\\binchen.mao\\Documents\\crownbio\\HuSignature\\HuPrimeModel.csv",h=T) #1691 PDX models
samples=unique(response$PDX)
samples[!samples %in% RNASeq_pdxs] 
samples[!samples %in% hu_prime_models$SAMPLE_NAME] 
samples=samples[samples %in% hu_prime_models$SAMPLE_NAME] 
ex_studied=ex[,sub("P[0-9]+$","",colnames(ex))%in% samples]
ex_studied_pdxs=sub("P[0-9]+$","",colnames(ex_studied))
duplicated_pdxs=ex_studied_pdxs[duplicated(ex_studied_pdxs)]
sort(names(ex_studied)[ex_studied_pdxs %in% duplicated_pdxs])
pdxs_to_remove=c('CR2506P3','CR2528P5','LU1656P9','OV0243P4','PA1194P5')
ex_studied=ex_studied[!names(ex_studied) %in% pdxs_to_remove]

colnames(ex_studied)=sub("P[0-9]+$","",colnames(ex_studied))
samples=samples[samples %in% colnames(ex_studied)]
ex_studied=ex_studied[,!duplicated(colnames(ex_studied))]
all_genes=gsub('\\|.*','',rownames(ex_studied))
ex_studied=ex_studied[!duplicated(all_genes),]
rownames(ex_studied)=gsub('\\|.*','',rownames(ex_studied))

rows <- apply(ex_studied , 1 , function(x) quantile(x,0.8)>=1 )
filteredExpr1 <- ex_studied[ rows , ]
filteredExpr3=filteredExpr1[!(grepl("-",row.names(filteredExpr1))),] 
prior=1
log_filteredExpr3=log(filteredExpr3+prior,2)
boxplot(log_filteredExpr3)

require(ggplot2)
mytheme <- theme(plot.title=element_text(face="bold.italic",size=16, color="brown"),
                 axis.title=element_text(face="bold.italic",size=14, color="brown"),
                 axis.text=element_text(face="bold", size=14,color="darkblue"),
                 axis.text.x=element_text(angle=0,face="bold", size=14,color="darkblue"),
                 panel.background=element_rect(fill="white",color="darkblue"),
                 panel.grid.major.y=element_line(color="grey",linetype=1),
                 panel.grid.minor.y=element_line(color="grey",linetype=2),
                 panel.grid.minor.x=element_blank(),
                 strip.text=element_text(face="bold.italic",size=16, color="brown"),
                 legend.position="top")
for (CT in levels(dat$CancerType)){
  dat_ct=subset(dat,CancerType==CT)
  n_pdx=length(unique(dat_ct$PDX))
  ggplot(dat_ct)+aes(Day,TV,color=Treatment,Group=as.factor(Mouse))+geom_point()+geom_line()+facet_wrap(~PDX,scales='free',ncol=4)+mytheme
  ggsave(paste0(CT,'_Growth_curve.png'),height=ceiling(n_pdx/4)*5,width=20,units='in',dpi=150,limitsize = F)
}  
dat %>% filter(Day==0 & Treatment=='Vehicle') %>% ggplot(aes(CancerType,TV))+geom_boxplot()

boxplot(response$TGI~response$CT)
boxplot(response$RTV_tc~response$CT)
boxplot(response$RTV_tr~response$CT)


#linear mixed model analysis
require(nlme)
lmm1=lme(log(TV+1)~Day+log(TV0)+Day:(Treatment*CancerType),random = list(PDX = ~Day+1,Mouse = ~Day+1),
         data=dat1, method='REML',control = lmeControl(opt = "optim"))
summary(lmm1) #Growth rate (OV>CR>PA>LU) under no treatment #drug inhibition LU>OV~=PA>CR
lmm2=lme(log(TV+1)~Day+log(TV0)+Day:Treatment,random = list(PDX = ~Day+1,Mouse = ~Day+1),
         data=dat1, method='REML',control = lmeControl(opt = "optim"))
randeff=random.effects(lmm2)[[1]]
cor(randeff[,2],response[rownames(randeff),'TGI'],method='spearman')
plot(randeff[,2],TGI[rownames(randeff)],xlab='randeff',ylab='TGI')
lmm3=lme(log(TV+1)~Day+log(TV0),random = list(PDX = ~Day+1,Mouse = ~Day+1),
         data=dat1[dat1$Treatment=='Treatment',], method='REML',control = lmeControl(opt = "optim"))
randeff=random.effects(lmm3)[[1]]
randeff_df=cbind(randeff,response[rownames(randeff),c('TGI','AUC','RTV_tr')])
randeff_df$PDX=rownames(randeff_df)
randeff_df=randeff_df[order(randeff_df$Day),]
ggplot(randeff_df)+aes(reorder(PDX,Day),Day)+geom_bar(stat='identity')+mytheme+
  theme(axis.text.x=element_text(angle=90,face="bold", size=8,color="darkblue"))+xlab('PDX')+ylab('Random Effects(PDX)')
ggplot(randeff_df)+aes(reorder(PDX,Day),TGI)+geom_bar(stat='identity')+mytheme+
  theme(axis.text.x=element_text(angle=90,face="bold", size=8,color="darkblue"))+xlab('PDX')
ggplot(randeff_df)+aes(reorder(PDX,Day),RTV_tr)+geom_bar(stat='identity')+mytheme+
  theme(axis.text.x=element_text(angle=90,face="bold", size=8,color="darkblue"))+xlab('PDX')
ggplot(randeff_df)+aes(reorder(PDX,Day),AUC)+geom_bar(stat='identity')+mytheme+
  theme(axis.text.x=element_text(angle=90,face="bold", size=8,color="darkblue"))+xlab('PDX')

save(list=c('dat1','log_filteredExpr3'),file='lmm.Rdata')
#fit linear mixed model using gene expression

###regression modeling & classification modeling
log_ex3=as.data.frame(t(log_filteredExpr3))
tgi_df=cbind('TGI'=response[rownames(log_ex3),'TGI'],log_ex3)
res_ex_df=cbind(response[rownames(log_ex3),],log_ex3)
res_df=cbind('Responder'=response[rownames(log_ex3),'Responder'],log_ex3)
res_df=subset(res_df,Responder!=2)
res_df$Responder=factor(ifelse(res_df$Responder==1,'Responder','Non_Responder'),levels=c('Responder','Non_Responder'))
ggplot(res_ex_df)+aes(BTG2,TGI)+geom_point()+geom_smooth(method = 'lm',se=F)+facet_wrap(~factor(substr(rownames(res_ex_df),1,2)))
ggplot(res_ex_df)+aes(response,BTG2)+geom_point()
tgi_genes=feature_sel(tgi_df,iter=200,subsample = 0.7,nrep = 10)
tgi_full_sel=tgi_df[,c('TGI',tgi_genes$sel_genes)]
featurePlot(x =tgi_full_sel[,-1],y =tgi_full_sel[,1], between = list(x = 1, y = 1),type = c("g", "p","smooth")) 
#HIST***, part of nucleosome, higher expression, better efficacy
cor(tgi_full_sel)
require(corrplot)
corrplot(cor(tgi_full_sel,use='pairwise.complete.obs',method='spearman'), 
         order="hclust",hclust.method="ward.D2",tl.col="black", tl.srt=45,addrect = 4, rect.col = "red")
res_ex_sel=cbind(response[rownames(log_ex3),],log_ex3[,tgi_genes$sel_genes])
cor(res_ex_sel[,2:5],res_ex_sel[,-(1:8)])
ggplot(res_ex_df)+aes(response,HIST2H4A)+geom_point()+geom_boxplot()
ggplot(res_ex_df)+aes(factor(Responder),NPDC1)+geom_boxplot(notch=T,varwidth = T)+geom_jitter(aes(color=CT),width=0.1)+mytheme
ggplot(res_ex_df)+aes(NPDC1,TGI,color=CT)+geom_point()+geom_smooth(method='rlm',se=F)+mytheme
#ggplot(res_ex_df)+aes(response,NPDC1)+geom_boxplot()+geom_point(position = 'dodge')
###ingore control, RTV_tr enrichment analysis (ribosome, peptide elongation->RPS** poorer efficacy)
ggplot(res_ex_df)+aes(response,BCR)+geom_point()+geom_boxplot()
ggplot(res_ex_df)+aes(EGFR,RTV_tr)+geom_point()+geom_smooth()

###
cor_all=cor(response[rownames(log_ex3),2:6],log_ex3,method='spearman')
write.csv(t(cor_all),'genes_eff_cor_all.csv')

###model training
trainModel=function(df,nrep=5,fs=FALSE,ROC=TRUE){
  result=list()
  for(i in 1:nrep){
    train_idx <- createDataPartition(df[,1], p = 0.7, list = FALSE,  times = 1)
    train_df=df[train_idx,]
    test_df=df[-train_idx,]
    train_sel=train_df
    test_sel=test_df
    if(fs){
      filtered_genes=feature_sel(train_df,min_f=5,iter=200,nrep=10,m='spearman',subsample = 0.7,consistency = 0.7)
      train_sel=train_df[,c(names(df)[1],filtered_genes$sel_genes)]
      test_sel=test_df[,c(names(df)[1],filtered_genes$sel_genes)]
      #filtered_genes=feature_sel(train_df,min_f=10,sig=0.01,iter=1000,nrep=1,m='spearman',subsample = 1)
      #train_sel=train_df[,c(names(df)[1],filtered_genes$high_cor_genes)]
      #test_sel=test_df[,c(names(df)[1],filtered_genes$high_cor_genes)]
    }	
    print(featurePlot(x = train_sel[,-1],y = train_sel[,1], between = list(x = 1, y = 1),type = c("g", "p","smooth")))
    #featurePlot(x = test_sel[,-1],y = test_sel[,1], between = list(x = 1, y = 1),type = c("g", "p"))
    models=buildModels(train_sel,ROC=ROC,all=FALSE)
    predicted=rbind(observed=test_sel[,1],do.call(rbind,lapply(lapply(models$models,predict,test_sel[,-1]),as.vector)))
    result[[i]]=list('model'=models,'prediction'=predicted,'test'=test_sel)
    print(featurePlot(x = test_sel[,-1],y = test_sel[,1], between = list(x = 1, y = 1),type = c("g", "p","smooth")))
  }	
  result
}
reg_models_genes_fsr=trainModel(tgi_df,nrep=10,fs=TRUE,ROC=FALSE)
save(reg_models_genes_fsr,file='tgi_model_genes.Rdata')
cf_models_genes_fsr=trainModel(res_df,nrep=2,fs=TRUE,ROC=TRUE)

require(GSVA)
library(GSEABase)
library(GSVAdata)
C2_cp=getGmt('c2.cp.v6.1.symbols.gmt',collectionType=BroadCollection(category="c2"),
             geneIdType=SymbolIdentifier())
ss2=gsva(as.matrix(ex_studied),C2_cp,mx.diff=TRUE,method="gsva",parallel.sz=4,kernel=TRUE)
pw_mat=ss2$es.obs
effm=TGI #use TGI as efficay measurement
#effm=wtmc
pw_cor=rev(sort(apply(pw_mat,1,function(x) cor(x,effm[colnames(pw_mat)],method='spearman'))))
top10=names(rev(sort(abs(pw_cor)))[1:30])
pw_cor[top10]
rand_cor=rev(sort(apply(pw_mat,1,function(x) cor(x,sample(effm),method='spearman'))))
summary(rand_cor)
boxplot(pw_cor,rand_cor) #significantly better corrleation than permutated correlations


#topology based pathway score
library(PerPAS)
require("biomaRt")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
gsInfo=getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),filters = 'hgnc_symbol', values = rownames(ex_studied),mart=ensembl)
eids=gsInfo[[2]]
names(eids)=gsInfo[[1]]
ex_studied_ens=ex_studied
rownames(ex_studied_ens)=make.names(eids[rownames(ex_studied)],unique = T)
norm.ex_studied=lcy.tORz.score(data=log2(ex_studied_ens+1),group=NULL,byrow=T,scale = T,method='median',type='zscore')
pas.score=lcy.pathway.scoring(data=norm.ex_studied,score.type = 'PDT')
pw_cor=rev(sort(apply(pas.score,1,function(x) cor(x,effm[colnames(pas.score)],method='spearman'))))
top10=names(rev(sort(abs(pw_cor)))[1:10])
pw_cor[top10]
rand_cor=rev(sort(apply(pas.score,1,function(x) cor(x,sample(effm),method='spearman'))))
boxplot(pw_cor,rand_cor)

#pathway modeling
tgi_pwys_df=as.data.frame(cbind('TGI'=response[colnames(pas.score),'TGI'],t(pas.score),t(pw_mat)))
#save(tgi_pwys_df,file='tgi_pwys_df.Rdata')
colnames(tgi_pwys_df)=make.names(colnames(tgi_pwys_df))
tgi_pwys=feature_sel(tgi_pwys_df,sig=0.01,min_f=5,iter=200,nrep=10,subsample = 0.7)
tgi_sel=tgi_pwys_df[,c('TGI',tgi_pwys$sel_genes)]
cor(tgi_sel,method='spearman')[,1]
ct=substr(rownames(tgi_sel),1,2)
x=by(tgi_sel,ct,cor,method='spearman')
featurePlot(x =tgi_sel[,-1],y =tgi_sel[,1], between = list(x = 1, y = 1),type = c("g", "p",'smooth'))
tgi_models_pwy_fsr=trainModel(tgi_pwys_df,nrep=10,fs=TRUE,ROC=FALSE)
save(tgi_models_pwy_fsr,file='tgi_model_pwys.Rdata')
tgi_models_pwy=trainModel(tgi_sel,nrep=1,fs=FALSE,ROC=FALSE)

#####cancer specific correlation
CR_tgi_pwys_df=tgi_pwys_df[substr(rownames(tgi_pwys_df),1,2)=='CR',]
tgi_pwys=feature_sel(CR_tgi_pwys_df,sig=0.01,min_f=5,iter=200,nrep=1,subsample = 1,m='spearman')
tgi_sel=CR_tgi_pwys_df[,c('TGI',tgi_pwys$sel_genes)]
cor(tgi_sel,method='spearman')[,1]
LU_tgi_pwys_df=tgi_pwys_df[substr(rownames(tgi_pwys_df),1,2)=='LU',]
tgi_pwys=feature_sel(LU_tgi_pwys_df,sig=0.01,min_f=5,iter=200,nrep=1,subsample = 1,m='spearman')
tgi_sel=LU_tgi_pwys_df[,c('TGI',tgi_pwys$sel_genes)]
cor(tgi_sel,method='spearman')[,1]
tgi_CR_models_pwy_fsr=trainModel(CR_tgi_pwys_df,nrep=10,fs=TRUE,ROC=FALSE)
tgi_LU_models_pwy_fsr=trainModel(LU_tgi_pwys_df,nrep=10,fs=TRUE,ROC=FALSE)

###DEG analysis
deg=function(x,f){
  t=t.test(x~f)
  w=wilcox.test(x~f)
  c(t$estimate,t$statistic,'pval'=t$p.value,'pval_mwu'=w$p.value)
}
res_pdxs=as.character(response[response$Responder != 2,'PDX'])
class=factor(response[res_pdxs,'Responder'])
#res_pdxs=response[response$response %in% c('CR','CPD'),'PDX']
#class=factor(response[res_pdxs,'Responder'])
Cohort_degTable=apply(log_filteredExpr3[,res_pdxs],1,function(x) deg(x,class))
Cohort_degTable=as.data.frame(t(Cohort_degTable))
Cohort_degTable$diff=Cohort_degTable[,2]-Cohort_degTable[,1] #difference of log2(RSEM+1)
Cohort_degTable$gs=rownames(Cohort_degTable)
Cohort_degTable$adj_p = p.adjust(Cohort_degTable$pval_mwu,method='fdr')
Cohort_Pos_high=Cohort_degTable %>% filter(diff>0 & adj_p<0.05)
Cohort_Pos_low=Cohort_degTable %>% filter(diff<0 & adj_p<0.05)
plot(Cohort_degTable$diff,-log10(Cohort_degTable$pval_mwu),xlab='log_RSEM_difference',ylab='-log10(pvalue)')
dev.copy(png,paste0('DEG',"_valcano_plot.png"))
dev.off()
require(ReactomePA)
gs2Reactome=function(gs,kw='unknown',pval=0.01){
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
Cohort_pos_up_gs=Cohort_Pos_high$gs
gs2Reactome(Cohort_pos_up_gs,paste0('DEG','_POS_UP'))
Cohort_pos_down_gs=Cohort_Pos_low$gs
gs2Reactome(Cohort_pos_down_gs,paste0('DEG','_POS_DOWN'))
require(fgsea)
t_val = Cohort_degTable$t
names(t_val)=rownames(Cohort_degTable)
#t_val = cor_all['TGI',]
t_val=t_val[order(t_val)] 
cp_pathways = gmtPathways(gmt.file = "c2.cp.v6.1.symbols.gmt")
pathways=cp_pathways
fgseaRes <- fgsea(pathways, t_val, minSize=5, maxSize=500, nperm=1000)
#fgseaRes <- fgsea(pathways, t_val, minSize=10, maxSize=30, nperm=1000)
head(fgseaRes[order(pval), ],n=20)
fgseaRes=subset(fgseaRes,padj<0.05)
pwinfo=as.data.frame(subset(fgseaRes,padj < 0.05))
write.csv(pwinfo[,1:7],'tgi_cor_sig_pathways_gsea.csv',row.names=F)
topPathwaysUp <- fgseaRes[ES > 0][head(rev(order(NES)), n=20), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(NES), n=20), pathway]
topPathways <- c(topPathwaysUp[1:10], topPathwaysDown[1:10])
topPathways=topPathways[!is.na(topPathways)]
#Ribosome(translation up)->non-responder; Integrein,centromere->responder
plotGseaTable(pathways[topPathways], t_val, fgseaRes,gseaParam = 0.5)
plotGseaTable(pathways[topPathwaysUp], t_val, fgseaRes,gseaParam = 0.5)
plotGseaTable(pathways[topPathwaysDown], t_val, fgseaRes,gseaParam = 0.5)
#use consensus cluster to separate positive & negative genes in a pathway
log_ex3_rank=apply(log_ex3,1,function(x) rank(-x))
pwy='PID_INTEGRIN1_PATHWAY'
pwy='REACTOME_G_ALPHA_Z_SIGNALLING_EVENTS'
pwy='REACTOME_DEPOSITION_OF_NEW_CENPA_CONTAINING_NUCLEOSOMES_AT_THE_CENTROMERE'
pwy='KEGG_RIBOSOME'
pwy_genes=intersect(cp_pathways[[pwy]],rownames(log_ex3_rank))
pwy_ex=log_ex3_rank[pwy_genes,]
#cor_pwy_ex=cor(pwy_ex,method='spearman')
#corrplot(cor_pwy_ex,order="hclust",hclust.method="ward.D2",tl.col="black", tl.srt=45,addrect = 2, rect.col = "red")
#require(ConsensusClusterPlus)
#results = ConsensusClusterPlus(as.matrix(t(pwy_ex)),maxK=3,reps=50,pItem=0.8,pFeature=1,
#                               title=paste0(pwy,"_consensus_cluster"),
#                               clusterAlg="hc",distance="spearman",plot="png")
#icl = calcICL(results,plot="png",title=paste0(pwy,"_consensus_cluster"))
#pwy_genes_class=results[[2]]$consensusClass
#pwy_genes_class[pwy_genes_class==2]=-1
#plot(pwy_genes_class,t_val[pwy_genes]) #correlation correlation well with cluster
#pwy_score=apply(pwy_ex,2,function(x) sum(x*pwy_genes_class))
#cor(pwy_score,TGI[names(pwy_score)],method='spearman')
pwy_score=apply(pwy_ex,2,function(x) sum(x))
cor(pwy_score,TGI[names(pwy_score)],method='spearman')
pwy_genes_class2=ifelse(t_val[pwy_genes]>0,1,-1)
pwy_score=apply(pwy_ex,2,function(x) sum(x*pwy_genes_class2))
cor(pwy_score,TGI[names(pwy_score)],method='spearman')
pwy_score3=apply(pwy_ex,2,function(x) sum(x*t_val[names(x)]))
cor(pwy_score3,TGI[names(pwy_score)],method='spearman')
###mutation(pathway aggreated mutation analysis, especially cancer ralated signal pathway, filter pathway genes using oncogene/tumor suppress gene)
###variants filtering (confidence, dbsnp, 1000genome,...)
mut_dir='mutation'
files=list.files(mut_dir) 
csv_files=files[grep('table$',files)]
kwfs = gsub('(\\w+)\\-.*$','\\1',csv_files)
names(csv_files)=kwfs
mut_df=list()
snv_df=list()
#filter variants (functional important variants kept)
for(kwf in kwfs){
  mut_pdx=read.table(paste0(mut_dir,'/',csv_files[kwf]),h=T,sep="\t",stringsAsFactors = F)
  mut_pdx=subset(mut_pdx,CrownFIL=='PASS')
  mut_pdx$esp6500siv2_all=as.numeric(mut_pdx$esp6500siv2_all)
  mut_pdx$X1000g2014oct_all=as.numeric(mut_pdx$X1000g2014oct_all)
  mut_pdx=subset(mut_pdx,esp6500siv2_all<0.01 | is.na(esp6500siv2_all))
  mut_pdx=subset(mut_pdx,X1000g2014oct_all<0.01 | is.na(X1000g2014oct_all))
  #mut_pdx=subset(mut_pdx,!(is.na(esp6500siv2_all) & !is.na(X1000g2014oct_all) & X1000g2014oct_all>0.01))
  mut_pdx=subset(mut_pdx,!(ExonicFunc.ensGene %in%  c('synonymous SNV','.','unknown')))
  mut_pdx=subset(mut_pdx,!(ExonicFunc.ensGene=='nonsynonymous SNV' & (RadialSVM_pred!='D' | SIFT_pred!='D')))
  mut_pdx=subset(mut_pdx,!(ExonicFunc.ensGene %in% c('nonframeshift insertion','nonframeshift deletion') & End-Start<9))
  snv_pdx=subset(mut_pdx,ExonicFunc.ensGene == 'nonsynonymous SNV')
  mut_df[[kwf]]=unique(unlist(sapply(mut_pdx$GeneSymbol,strsplit,',')))
  snv_df[[kwf]]=unique(unlist(sapply(snv_pdx$GeneSymbol,strsplit,',')))
}  
save(list=c('mut_df','snv_df'),file='annovar_mut.Rdata')
all_mut=sort(table(unlist(mut_df)))
n_mut=sapply(mut_df,length)
n_nsm=sapply(snv_df,length)
all_genes=unique(unlist(mut_df))
all_mut_df=do.call(rbind,lapply(mut_df,function(x){v=numeric(length(all_genes));names(v)=all_genes;v[x]=1;v}))
res_all_mut_df=as.data.frame(cbind('Responder'=response[rownames(all_mut_df),'Responder'],all_mut_df))
#CTs=substr(rownames(res_all_mut_df),1,2)
#by(res_all_mut_df,CTs,function(x) table(x$Responder,x$KRAS))

res_mut_df=subset(res_all_mut_df,Responder != 2)
ae=apply(res_mut_df,2,function(x) length(unique(x))==1)
res_mut_df=res_mut_df[,!ae]
mut_pval=sapply(res_mut_df[,-1],function(x) chisq.test(x,res_mut_df$Responder)$p.value)
mut_pval=sort(mut_pval)

KRAS_mut=read.csv('KRAS_mutataion.csv',h=T,stringsAsFactors = F)
KRAS_G13D_G12V=dplyr::filter(KRAS_mut,Mutationtype=='Nonsynonymousmutation' & Source %in% c('Exome','PCR based')
                 & Amino_acid_change %in% c('G13D','G12V'))
KRAS_1213=dplyr::filter(KRAS_mut,Mutationtype=='Nonsynonymousmutation' & Source %in% c('Exome','PCR based')
                 & substr(Amino_acid_change,2,3) %in% c(12,13))
KRAS_mut_models = intersect(response$PDX,KRAS_1213$PDXmodel)
TGI[KRAS_mut_models]
response_KRAS= merge(response,KRAS_1213[,c('PDXmodel','Amino_acid_change')],by.x='PDX',by.y='PDXmodel',all.x=T)
response_KRAS=response_KRAS[!duplicated(response_KRAS),]
response_KRAS[is.na(response_KRAS$Amino_acid_change),'Amino_acid_change']='WT'
ggplot(response_KRAS)+aes(response,TGI,color=Amino_acid_change)+geom_jitter(width=0.2)+facet_wrap(~CT)+scale_color_brewer(palette='Set1')
write.csv(response_KRAS,'response_KRAS_mut.csv')
PA1=subset(response_KRAS,CT=='PA')

####Variants annotation from VEP
mut_dir='mutation/VEP/mutation'
files=list.files(mut_dir) 
csv_files=files[grep('table$',files)]
kwfs = gsub('(\\w+)P\\d+.*$','\\1',csv_files)
names(csv_files)=kwfs
driver_mut_df=list()
#filter variants (functional important variants kept)
#all_pdx_driver_mut=data.frame()
for(kwf in kwfs){
  mut_pdx=read.table(paste0(mut_dir,'/',csv_files[kwf]),h=T,sep="\t",stringsAsFactors = F,comment.char = "!")
  mut_pdx=subset(mut_pdx,DriverPrediction!='.')
  #mut_pdx=cbind('PDX'=kwf,mut_pdx)
  #all_pdx_driver_mut=rbind(all_pdx_driver_mut,mut_pdx)
  #mut_pdx=NULL
  driver_mut_df[[kwf]]=unique(mut_pdx$SYMBOL)
}  
#write.csv(all_pdx_driver_mut,'all_pdx_driver_mutation.csv')
all_driver_mut=rev(sort(table(unlist(driver_mut_df))))
n_driver_mut=sapply(driver_mut_df,length)
all_genes=unique(unlist(driver_mut_df))
driver_mut_df=do.call(rbind,lapply(driver_mut_df,function(x){v=numeric(length(all_genes));names(v)=all_genes;v[x]=1;v}))
res_driver_mut_df=as.data.frame(cbind('Responder'=response[rownames(driver_mut_df),'Responder'],driver_mut_df))
res_driver_mut_df=subset(res_driver_mut_df,Responder != 2)
ae=apply(res_driver_mut_df,2,function(x) length(unique(x))==1)
res_driver_mut_df=res_driver_mut_df[,!ae]
driver_mut_pval=sapply(res_driver_mut_df[,-1],function(x) chisq.test(x,res_driver_mut_df$Responder)$p.value)
driver_mut_pval=sort(driver_mut_pval)
response_mut=cbind(response,'TP53_mutation'=driver_mut_df[rownames(response),'TP53'])

#calculate number of genes carrying delterious mutation in each pathway
acsn_pwys = gmtPathways(gmt.file = "acsn_v1.1.gmt")
mut_genes=intersect(colnames(res_driver_mut_df)[-1],C2_cp[['KEGG_PATHWAYS_IN_CANCER']]@geneIds)
all_pw_mut_df=do.call(rbind,lapply(C2_cp,function(x) unlist(apply(driver_mut_df,1,function(v) sum(v[intersect(x@geneIds,mut_genes)])))))
rownames(all_pw_mut_df)=names(C2_cp)
response_mut=cbind(response_mut,'beta_catenin_pathway_mut'=all_pw_mut_df['ST_WNT_BETA_CATENIN_PATHWAY',rownames(response_mut)])
response_mut=cbind(response_mut,'Centromere'=unlist(pwys_df['REACTOME_DEPOSITION_OF_NEW_CENPA_CONTAINING_NUCLEOSOMES_AT_THE_CENTROMERE',rownames(response_mut)]))
ggplot(response_mut)+aes(factor(Responder),Centromere)+geom_boxplot(notch=T,varwidth = T)+
  geom_jitter(aes(color=factor(TP53_mutation),size=beta_catenin_pathway_mut),width=0.1)+facet_wrap(~substr(rownames(response_mut),1,2))+mytheme

ggplot(response_mut)+aes(Centromere,TGI)+geom_point(aes(color=factor(TP53_mutation),size=beta_catenin_pathway_mut))+
  geom_smooth(method='rlm',se=F)+facet_wrap(~substr(rownames(response_mut),1,2))+mytheme

pw_mut_df=do.call(rbind,lapply(C2_cp,function(x) unlist(apply(res_driver_mut_df[,-1],1,function(v) sum(v[intersect(x@geneIds,mut_genes)])))))
rownames(pw_mut_df)=names(C2_cp)
t.val=apply(pw_mut_df,1,function(x) t.test(x~res_driver_mut_df$Responder)$p.value)
pwy='ST_WNT_BETA_CATENIN_PATHWAY'
table(pw_mut_df[pwy,],res_driver_mut_df$Responder)
intersect(C2_cp[[pwy]]@geneIds,mut_genes)
pw_mut_df_binary=pw_mut_df #binary (either mutated or non-mutated)
pw_mut_df_binary[pw_mut_df_binary>0]=1
ae=apply(pw_mut_df_binary,1,function(x) length(unique(x))==1)
pw_mut_pval=sort(apply(pw_mut_df_binary[!ae,],1,function(x) chisq.test(x,res_driver_mut_df$Responder)$p.value))
table(pw_mut_df_binary['BIOCARTA_RNA_PATHWAY',],res_driver_mut_df$Responder)

acsn_pw_mut_df=do.call(rbind,lapply(acsn_pwys,function(x) unlist(apply(res_driver_mut_df[,-1],1,function(v) sum(v[intersect(x,mut_genes)])))))
t.val=apply(acsn_pw_mut_df,1,function(x) t.test(x~res_driver_mut_df$Responder)$p.value)
acsn_pw_mut_df[acsn_pw_mut_df>0]=1
ae=apply(acsn_pw_mut_df,1,function(x) length(unique(x))==1)
acsn_pw_mut_pval=sort(apply(acsn_pw_mut_df[!ae,],1,function(x) chisq.test(x,res_driver_mut_df$Responder)$p.value))
pwy='DESMOSOMES'
table(acsn_pw_mut_df[pwy,],res_driver_mut_df$Responder)
intersect(acsn_pwys[[pwy]],mut_genes) 
p53_pwys=sapply(acsn_pwys,function(x) 'TP53' %in% x)
acsn_pw_mut_pval[!p53_pwys[names(acsn_pw_mut_pval)]]
#the major problem is that mutation pattern is correlated with cancer type, and repsonse rate in different cancers are quite different
find_sig_ct=function(res_df){
  find_sig=function(res_df){
    ae=apply(res_df,2,function(x) length(unique(x))==1)
    res_df=res_df[,!ae]
    pval=sapply(res_df[,-1],function(x) chisq.test(x,res_df[,1])$p.value)
    pval[pval<0.05]
  }
  ct=substr(rownames(res_df),1,2)
  res_df_by_ct=by(res_df,ct,find_sig)
}  
gene_mut_ct=find_sig_ct(res_driver_mut_df)
#########################
CT=data.frame('ct'=substr(rownames(res_driver_mut_df),1,2))
x=dummyVars(~ct,data=CT)
ct_df=predict(x,CT)
sel_mut_genes=feature_sel(cbind(res_driver_mut_df,ct_df),sig=0.05,subsample = 1)
gene='PTPRB'
table(res_driver_mut_df$Responder,res_driver_mut_df[,gene],CT$ct)
filter(all_pdx_driver_mut,SYMBOL==gene)

###copy number variation
require(data.table)
cnv_all=fread('CN_gene_combined-2017-12-11.csv')
names(cnv_all)[1]='Symbol'
names(cnv_all)=gsub('P\\d+$','',names(cnv_all))
#require("biomaRt")
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",GRCh = 37)
#normal.chroms <- c(1:22, "X", "Y", "M")
#my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
#                    filters = c("hgnc_symbol", "chromosome_name"),
#                    values = list(hgnc_symbol=cnv_all$Symbol[1:5], chromosome_name=normal.chroms),
#                    mart = ensembl)
cnv_all=as.data.frame(cnv_all)
rownames(cnv_all)=cnv_all$Symbol
gene_info=cnv_all[,1:4]
cnv_info=cnv_all[,colnames(log_filteredExpr3)]
save(list=c('cnv_info','gene_info'),file='CNV_info.Rdata')
gene_cnv_ex_cor=sapply(intersect(rownames(cnv_info),rownames(log_filteredExpr3)),
                function(x) cor(unlist(cnv_info[x,]),unlist(log_filteredExpr3[x,]),method='spearman'))
cnv_genes=names(gene_cnv_ex_cor[gene_cnv_ex_cor>=0.3])
cnv_genes_info=cnv_info[cnv_genes,]
save(cnv_genes_info,file='cnv_genes_info.Rdata')
cor_all_cnv=cor(response[colnames(cnv_info),2:7],t(cnv_info[cnv_genes,]),method='spearman')
write.csv(t(cor_all_cnv),'cnv_eff_cor_all.csv')

cnv_genes_cg=round(cnv_genes_info) #CNV categorial 
cnv_genes_cg[cnv_genes_cg>4]=4
response_i=5-as.numeric(response[colnames(cnv_genes_cg),'response'])
cnv_cor_cg=sort(apply(cnv_genes_cg,1,cor,response_i))
gene='PTP4A1' #Protein tyrosine phosphatase which stimulates progression from G1 into S phase during mitosis;however,PTP4A2 not correlated
table(unlist(cnv_genes_cg[gene,]),response[colnames(cnv_genes_cg),'response'])





#fitting tumor growth model using double expontial model
dat2=dat %>% rename(name=Mouse,date=Day,size=TV)
out <- gdrate(dat2[,c('name','date','size')], 0.05, TRUE) #delayed kicked in/small TV change (not fitting well)
#fit linear mixed model 
require(nlme)

############exploratory analysis
tak_calc=read.csv('Takeda_calc.csv',h=T,stringsAsFactors = F)
tak_calc=subset(tak_calc,Model != '')
rownames(tak_calc)=tak_calc$Model
tgi_diff=TGI-tak_calc[names(TGI),'TGI']
tgi_df=data.frame('PDX'=names(TGI),'TGI'=TGI,'Takeda_TGI'=tak_calc[names(TGI),'TGI'],'diff'=tgi_diff)
tgi_df %>% arrange(desc(abs(tgi_diff)))


