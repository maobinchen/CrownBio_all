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

####gene expression 
response=read.csv('response_TAK931.csv',h=T)
#tak_calc=read.csv('Takeda_calc.csv',h=T)
#X=merge(response,tak_calc[,c('Model','GRI')],by.x = 'PDX',by.y='Model' )
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
ex_studied=ex_studied[,!colnames(ex_studied) %in% model_to_remove]

rows <- apply(ex_studied , 1 , function(x) quantile(x,0.8)>=1 )
filteredExpr1 <- ex_studied[ rows , ]
filteredExpr3=filteredExpr1[!(grepl("-",row.names(filteredExpr1))),] 
prior=1
log_filteredExpr3=log(filteredExpr3+prior,2)
boxplot(log_filteredExpr3)

save(list=c('dat1','log_filteredExpr3'),file='lmm.Rdata')

log_ex3=as.data.frame(t(log_filteredExpr3))
tgi_df=cbind('TGI'=response[rownames(log_ex3),'TGI'],log_ex3)
res_ex_df=cbind(response[rownames(log_ex3),],log_ex3)
res_df=cbind('Responder'=response[rownames(log_ex3),'Responder'],log_ex3)
res_df=subset(res_df,Responder!=2)
res_df$Responder=factor(ifelse(res_df$Responder==1,'Responder','Non_Responder'),levels=c('Responder','Non_Responder'))

cor_all=cor(response[rownames(log_ex3),2:7],log_ex3,method='spearman')

trainModel=function(df,nrep=5,fs=FALSE,ROC=TRUE){
  result=list()
  for(i in 1:nrep){
    train_idx <- createDataPartition(df[,1], p = 0.7, list = FALSE,  times = 1)
    train_df=df[train_idx,]
    test_df=df[-train_idx,]
    train_sel=train_df
    test_sel=test_df
    if(fs){
      filtered_genes=feature_sel(train_df,min_f=5,iter=200,nrep=10,m='spearman',subsample = 0.8,consistency = 0.8)
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
cf_models_genes_fsr=trainModel(res_df,nrep=10,fs=TRUE,ROC=TRUE)

require(GSVA)
library(GSEABase)
library(GSVAdata)
C2_cp=getGmt('c2.cp.v6.1.symbols.gmt',collectionType=BroadCollection(category="c2"),
             geneIdType=SymbolIdentifier())
C5_cc=getGmt('c5.cc.v6.1.symbols.gmt',collectionType=BroadCollection(category="c5"),
             geneIdType=SymbolIdentifier())
C5_mf=getGmt('c5.mf.v6.1.symbols.gmt',collectionType=BroadCollection(category="c5"),
             geneIdType=SymbolIdentifier())
acsn_pwys=getGmt('acsn_v1.1.gmt',
                 geneIdType=SymbolIdentifier())
all_gs=GeneSetCollection(c(C2_cp,C5_cc,acsn_pwys),geneIdType=SymbolIdentifier())
ss2=gsva(as.matrix(log_filteredExpr3),all_gs,mx.diff=TRUE,method="gsva",parallel.sz=4,kernel=TRUE)
pw_mat=ss2$es.obs
library(PerPAS)
require("biomaRt")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
gsInfo=getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),filters = 'hgnc_symbol', values = rownames(log_filteredExpr3),mart=ensembl)
eids=gsInfo[[2]]
names(eids)=gsInfo[[1]]
ex_studied_ens=log_filteredExpr3
rownames(ex_studied_ens)=make.names(eids[rownames(log_filteredExpr3)],unique = T)
norm.ex_studied=lcy.tORz.score(data=log2(ex_studied_ens+1),group=NULL,byrow=T,scale = T,method='median',type='zscore')
pas.score=lcy.pathway.scoring(data=norm.ex_studied,score.type = 'PDT')
pwys_df=as.data.frame(rbind(pas.score,pw_mat))
rownames(pwys_df)=make.names(rownames(pwys_df))
cor_all_pwys=as.data.frame(t(cor(response[colnames(pwys_df),2:7],t(pwys_df),method='spearman')))
write.csv(cor_all_pwys,'pwys_eff_cor_all.csv')

res_pwys_df=as.data.frame(cbind('Responder'=response[colnames(pas.score),'Responder'],t(pas.score),t(pw_mat)))
colnames(res_pwys_df)=make.names(colnames(res_pwys_df))
res_pwys_df=subset(res_pwys_df,Responder!=2)
res_pwys_df$Responder=factor(ifelse(res_pwys_df$Responder==1,'Responder','Non_Responder'),levels=c('Responder','Non_Responder'))
cf_pwys=feature_sel(res_pwys_df,sig=0.01,min_f=5,iter=200,nrep=10,subsample = 0.7)
cf_pwys_sel=res_pwys_df[,c('Responder',cf_pwys$sel_genes)]
featurePlot(x =cf_pwys_sel[,-1],y =cf_pwys_sel[,1], between = list(x = 1, y = 1),type = c("g", "p")) 
cf_models_pwys_fsr=trainModel(res_pwys_df,nrep=2,fs=TRUE,ROC=TRUE)
#cor_all_pwy=cor(response[colnames(pas.score),2:6],cbind(t(pw_mat),t(pas.score)),method='spearman')

res_genes_pwys_df=cbind(res_pwys_df,res_df[,-1])
cf_models_genes_pwys_fsr=trainModel(res_genes_pwys_df,nrep=10,fs=TRUE,ROC=TRUE)

###mutation
mut_dir='mutation/VEP/mutation'
files=list.files(mut_dir) 
csv_files=files[grep('table$',files)]
kwfs = gsub('(\\w+)P\\d+.*$','\\1',csv_files)
names(csv_files)=kwfs
driver_mut_df=list()
#filter variants (functional important variants kept)
for(kwf in kwfs){
  mut_pdx=read.table(paste0(mut_dir,'/',csv_files[kwf]),h=T,sep="\t",stringsAsFactors = F,comment.char = "!")
  mut_pdx=subset(mut_pdx,DriverPrediction!='.')
  driver_mut_df[[kwf]]=unique(mut_pdx$SYMBOL)
}  
all_genes=unique(unlist(driver_mut_df))
driver_mut_df=do.call(rbind,lapply(driver_mut_df,function(x){v=numeric(length(all_genes));names(v)=all_genes;v[x]=1;v}))
driver_mut_df=driver_mut_df[!rownames(driver_mut_df) %in% model_to_remove,]
#save(driver_mut_df,file='driver_mut_df.Rdata')
res_driver_mut_df=as.data.frame(cbind('Responder'=response[rownames(driver_mut_df),'Responder'],driver_mut_df))
res_driver_mut_df=subset(res_driver_mut_df,Responder != 2)
ae=apply(res_driver_mut_df,2,function(x) length(unique(x))==1)
res_driver_mut_df_nae=res_driver_mut_df[,!ae]
driver_mut_pval=sapply(res_driver_mut_df_nae[,-1],function(x) chisq.test(x,res_driver_mut_df_nae$Responder)$p.value)
driver_mut_pval=sort(driver_mut_pval)

cosmic=read.csv('Cosmic_consensus.csv',h=T,stringsAsFactors = F)
t1_genes=intersect(cosmic[cosmic$Tier==1,'Gene.Symbol'],all_genes) #573 cosmic tier1 cancer genes
geneset=all_gs
pw_mut_df=do.call(rbind,lapply(geneset,function(x) unlist(apply(res_driver_mut_df[,-1],1,function(v) sum(v[intersect(x@geneIds,t1_genes)])))))
rownames(pw_mut_df)=names(geneset)
pw_mut_df=pw_mut_df[complete.cases(pw_mut_df),]
t.val=sort(apply(pw_mut_df,1,function(x) t.test(x~factor(res_driver_mut_df$Responder))$p.value))
pwy='KEGG_REGULATION_OF_ACTIN_CYTOSKELETON'
pwy='PID_CDC42_PATHWAY'
pwy='GO_CELL_PROJECTION_MEMBRANE'
pwy='GO_CHROMOSOME_CENTROMERIC_REGION'
pwy='GO_DNA_HELICASE_ACTIVITY'
pwy='GO_CENTROSOME'
pwy='GO_MICROTUBULE_ORGANIZING_CENTER'
pwy='GO_CELL_LEADING_EDGE'
pwy='BIOCARTA_ALK_PATHWAY'
table(pw_mut_df[pwy,],res_driver_mut_df$Responder)
intersect(geneset[[pwy]]@geneIds,t1_genes)
pwys_mut_eff=t(sapply(names(t.val[t.val<0.01]),function(x) c(t.val[x],paste(intersect(geneset[[x]]@geneIds,t1_genes),collapse=","))))
colnames(pwys_mut_eff)=c('p.val','Cancer Driver Gene')
write.csv(pwys_mut_eff,'response_mutated_pwys.csv')

require(fgsea)
acsn_pwys = gmtPathways(gmt.file = "acsn_v1.1.gmt")
acsn_pw_mut_df=do.call(rbind,lapply(acsn_pwys,function(x) unlist(apply(res_driver_mut_df[,-1],1,function(v) sum(v[intersect(x,t1_genes)])))))
t.val=apply(acsn_pw_mut_df,1,function(x) t.test(x~res_driver_mut_df$Responder)$p.value)
acsn_pw_mut_df[acsn_pw_mut_df>0]=1
ae=apply(acsn_pw_mut_df,1,function(x) length(unique(x))==1)
acsn_pw_mut_pval=sort(apply(acsn_pw_mut_df[!ae,],1,function(x) chisq.test(x,res_driver_mut_df$Responder)$p.value))
pwy='G1_S_CHECKPOINT'
table(acsn_pw_mut_df[pwy,],res_driver_mut_df$Responder)
intersect(acsn_pwys[[pwy]],t1_genes)


####CNV data for ploidy information
cnv_info=read.csv('Absolute_ploidy.csv',h=T)
cnv_info$sample=gsub('P\\d+$','',cnv_info$sample)
ploidy_info=cnv_info$ploidy
names(ploidy_info)=cnv_info$sample

###MSI information
msi_info=read.csv('MSI_exome_seq.csv',h=T)
msi_info$Model=gsub('P\\d+$','',msi_info$Model)
msi_h=msi_info$Class
msi_s=msi_info$MSI.score
names(msi_h)=msi_info$Model
names(msi_s)=msi_info$Model



#sel_pwys='REACTOME_CTNNB1_PHOSPHORYLATION_CASCADE'

require(pheatmap)

sel_genes=union(c('TP53','CTNNB1','CHD8','CHD2','CHD3','CHD4'),names(driver_mut_pval[driver_mut_pval<0.05]))
sel_pwys=rownames(cor_all_pwys[order(cor_all_pwys$RTV_tc),][1:10,])
sel_pwys=sel_pwys[c(1,3,4,9)]
#sel_pwys=c(sel_pwys,'KEGG_HISTIDINE_METABOLISM')
sel_genes_mut=t(res_driver_mut_df)[sel_genes,]
sel_pwys_score=pwys_df[sel_pwys,rownames(res_driver_mut_df)]
sel_pwys_score=t(apply(sel_pwys_score,1,function(v) (v-min(v))/(max(v)-min(v))))
Cohort_class = ifelse(res_driver_mut_df$Responder==0,'Non_responder','Responder')
CT=factor(substr(rownames(res_driver_mut_df),1,2))

annotation_col = data.frame(
  Cancer_type = CT,
  Response=Cohort_class
)
ann_colors = list(
  Response = c('Non_responder'='firebrick3','Responder'='green'),
  Cancer_type=c('CR'='gray','LU'='cyan','PA'='yellow','OV'='purple')
)
rownames(annotation_col)=rownames(res_driver_mut_df)
pheatmap(rbind(sel_genes_mut,sel_pwys_score), show_colnames = T,
         clustering_distance_rows = 'minkowski',
         clustering_distance_cols = 'minkowski',
         cutree_cols = 3,
         annotation_col = annotation_col,annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.copy(png,paste0('Response','_mut_pwy.png'),width=1500,height=800, res=150)
dev.off()

###survival analysis
genSurvData=function(df){
  cutoff=600 #triple volume size as endpoint
  PDX=df[1,'PDX']
  treatment=as.character(df[1,'Treatment'])
  Mouse=df[1,"Mouse"]
  ind=which(df$TV>=cutoff)
  if(length(ind) > 0){
    i=ind[1]
    tri_t = df[i,'Day']
    c("PDX"=PDX,"Treatment"=treatment,"Mouse"=Mouse,"time"=tri_t,"censoring"=1)
  }else{
    c("PDX"=PDX,"Treatment"=treatment,"Mouse"=Mouse,"time"=df[nrow(df),'Day'],"censoring"=0)
  }
}
byMouse=by(dat,dat$Mouse,genSurvData)
surv=as.data.frame(do.call(rbind,byMouse))
require(survival)
require(survminer)
surv$time=as.numeric(as.character(surv$time))
surv <- within(surv, Treatment <- relevel(Treatment, ref = 'Vehicle'))
surv$censoring=as.numeric(as.character(surv$censoring))
surv$CT=substr(surv$PDX,1,2)
surv$Treatment=ifelse(surv$Treatment=='Vehicle',0,1)
s1=with(surv,Surv(time,censoring))
fit=survfit(s1~Treatment, data=surv)
coxph(s1~Treatment*CT,data=surv)
#plot(fit)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw())
surv_trt=subset(surv,Treatment==1)
s2=with(surv_trt,Surv(time,censoring))
fit2=survfit(s2~CT, data=surv_trt)
ggsurvplot(fit2,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw())  


#model building  

sel_genes=names(driver_mut_pval[driver_mut_pval<0.05])
sel_pwys=rownames(cor_all_pwys[order(cor_all_pwys$RTV_tc),][1:20,])
cf_pwys_mut=cbind(res_driver_mut_df[,c('Responder',sel_genes)],res_pwys_df[rownames(res_driver_mut_df),sel_pwys])
cf_pwys_mut$Responder=factor(ifelse(cf_pwys_mut$Responder==1,'Responder','Non_Responder'),levels=c('Responder','Non_Responder'))
CT=data.frame('ct'=substr(rownames(cf_pwys_mut),1,2))
x=dummyVars(~ct,data=CT)
ct_df=predict(x,CT)

all_pwys_mut=cbind(driver_mut_df[,sel_genes],t(pwys_df)[rownames(driver_mut_df),sel_pwys])
sel_all_df=as.data.frame(all_pwys_mut[,c('TP53',"G1_S_CHECKPOINT","REACTOME_DEPOSITION_OF_NEW_CENPA_CONTAINING_NUCLEOSOMES_AT_THE_CENTROMERE")])
sel_all_df$Non_TP53_mut_sum=apply(all_pwys_mut[,2:20],1,sum)
sel_all_df$CHD_mut=apply(driver_mut_df[rownames(sel_all_df),c('CHD1L','CHD8','CHD2','CHD3','CHD4')],1,sum)
sel_all_df=cbind('TGI'=response[rownames(sel_all_df),'TGI'],sel_all_df)
##models not used for biomarker selection
sel_test_df=sel_all_df[setdiff(rownames(sel_all_df),rownames(res_driver_mut_df)),] 
sel_test_df$Hypermutated=ifelse(sel_test_df$Non_TP53_mut_sum>=2,1,0)
write.csv(sel_test_df,'15models_test.csv')
names(sel_test_df)[4]='CENTROMERE'
require(corrgram)
corrgram(sel_test_df[,c(1,2,3,4,7)],lower.panel=panel.cor, upper.panel=panel.pts,
         diag.panel=panel.density,cor.method = 'spearman')
corrgram(sel_test_df[rownames(sel_test_df) != 'LU1380',c(1,2,3,4,7)],lower.panel=panel.cor, upper.panel=panel.pts,
         diag.panel=panel.density,cor.method = 'spearman')



require(caret)
fiveStats=function(...) c(twoClassSummary(...), defaultSummary(...))
newRF <- rfFuncs
newRF$summary <- fiveStats
index <- createMultiFolds(cf_pwys_mut[,1], times = 5)
#varSeq <- 1:(ncol(cf_pwys_mut)-1)
varSeq <- 1:20
ctrl <- rfeControl(method = "repeatedcv",repeats = 5,verbose = TRUE,functions = newRF,index = index)
rfRFE <- rfe(x = cbind(cf_pwys_mut[,-1],ct_df), y = cf_pwys_mut[,1],sizes = varSeq,metric = "Accuracy",rfeControl = ctrl,ntree = 1000)

cf_pwys=feature_sel(cf_pwys_mut,sig=0.05,min_f=5,iter=200,nrep=10,subsample = 0.7)
cf_sel=cf_pwys_mut[,c('Responder',cf_pwys$sel_genes)]
featurePlot(x =cf_sel[,-1],y =as.factor(cf_sel[,1]), between = list(x = 1, y = 1),type = c("g", "p"))


sel_df=cf_pwys_mut[,c('Responder','TP53',"G1_S_CHECKPOINT","REACTOME_DEPOSITION_OF_NEW_CENPA_CONTAINING_NUCLEOSOMES_AT_THE_CENTROMERE")]
sel_df$Non_TP53_mut_sum=apply(cf_pwys_mut[,3:21],1,sum)
sel_df$Centrosome_mut=pw_mut_df['GO_CENTROSOME',rownames(sel_df)]
sel_df$CHD_mut=apply(res_driver_mut_df[rownames(sel_df),c('CHD1L','CHD8','CHD2','CHD3','CHD4')],1,sum)
names(sel_df)[4]='CENTROMERE'
library(rpart)
fit=rpart(Responder~.,data=sel_df,method='class',control = rpart.control(cp = 0.01))
library(rattle)
fancyRpartPlot(fit)
fit2=rpart(Responder~.,data=sel_df[,-5],method='class',control = rpart.control(cp = 0.01))
fancyRpartPlot(fit2)
fit3=rpart(Responder~.,data=sel_df[,c(1,2,5)],method='class',control = rpart.control(minsplit = 10))
fancyRpartPlot(fit3)
#plot(fit,uniform=T)
#text(fit,use.n=T,all=T,cex=0.6)

sel_df$ploidy=ploidy_info[rownames(sel_df)]
boxplot(sel_df$ploidy~sel_df$TP53)
t.test(sel_df$ploidy~sel_df$TP53)
t.test(sel_df$ploidy~sel_df$Responder)

sel_df$MSI_H=msi_h[rownames(sel_df)]
sel_df$MSI_Score=msi_s[rownames(sel_df)]
table(sel_df$Responder,sel_df$MSI_H)
ggplot(sel_df)+aes(Responder,MSI_Score,col=MSI_H,shape=factor(TP53))+geom_jitter(width=0.1)+
  mytheme+facet_wrap(~substr(rownames(sel_df),1,2))
dev.copy(png,'MSI_H_profile.png',width=1200,height=1200, res=150)
dev.off()


CT=factor(substr(rownames(sel_df),1,2))

annotation_col = data.frame(
  Cancer_type = as.character(CT),
  Response=sel_df$Responder
)
ann_colors = list(
  Response = c('Non_Responder'='firebrick3','Responder'='green'),
  Cancer_type=c('CR'='gray','LU'='cyan','PA'='yellow','OV'='purple')
)
rownames(annotation_col)=rownames(sel_df)
sel_df_binary=sel_df[,c(2:4,9)]
names(sel_df_binary)[1]='TP53_mut'
sel_df_binary$Hypermutated=ifelse(sel_df$Non_TP53_mut_sum>=2,1,0)
sel_df_binary$CHD_mut=ifelse(sel_df$CHD_mut>=1,1,0)
min_max_norm=function(v) (v-min(v))/(max(v)-min(v))
sel_df_binary$G1_S_CHECKPOINT=min_max_norm(sel_df_binary$G1_S_CHECKPOINT)
sel_df_binary$CENTROMERE=min_max_norm(sel_df_binary$CENTROMERE)
sel_df_binary$MSI_H=ifelse(sel_df_binary$MSI_H=='non-MSI-H',0,1)
pheatmap(t(sel_df_binary), show_colnames = T,
         clustering_distance_rows = 'minkowski',
         clustering_distance_cols = 'minkowski',
         cutree_cols = 3,
         annotation_col = annotation_col,annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.copy(png,paste0('Response','_mut_pwy.png'),width=1500,height=800, res=150)
dev.off()




rfRFE <- rfe(x = sel_df[,-1], y = sel_df[,1],sizes = varSeq,metric = "Accuracy",rfeControl = ctrl,ntree = 1000)

sel_df$PDX=rownames(sel_df)


surv_sel=merge(surv_trt,sel_df,by='PDX')
s2=with(surv_sel,Surv(time,censoring))
fit2=survfit(s2~TP53, data=surv_sel)
fit3=survfit(s2~(Non_TP53_mut_sum>=2),data=surv_sel)
fit4=survfit(s2~(Centrosome_mut>=2),data=surv_sel)
surv_sel_tp53m=subset(surv_sel,TP53==1)
s3=with(surv_sel_tp53m,Surv(time,censoring))
fit5=survfit(s3~cut_number(G1_S_CHECKPOINT,2),data=surv_sel_tp53m)
fit5=survfit(s3~Centrosome_mut,data=surv_sel_tp53m)
ggsurvplot(fit5,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw())  

sel_test_df$PDX=rownames(sel_test_df)
surv_test_sel=merge(surv_trt,sel_test_df,by='PDX')
s4=with(surv_test_sel,Surv(time,censoring))
fit5=survfit(s4~cut_number(CENTROMERE,2),data=surv_test_sel)
fit5=survfit(s4~Hypermutated,data=surv_test_sel)
ggsurvplot(fit5,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw())  

##########identify secondary biomarker for TP53_mut
TP53_mut_pdxs=rownames(sel_df[sel_df$TP53==1,])
TP53m_df_g=log_filteredExpr3[,TP53_mut_pdxs] #genes
TP53m_df_p=pwys_df[,TP53_mut_pdxs] #pathways
TP53m_res=response[TP53_mut_pdxs,'Responder']
cor_all_genes=as.data.frame(t(cor(response[TP53_mut_pdxs,2:7],t(TP53m_df_g),method='spearman')))
write.csv(cor_all_genes,'TP53mut_cor_genes.csv')
cor_all_pwys=as.data.frame(t(cor(response[TP53_mut_pdxs,2:7],t(TP53m_df_p),method='spearman')))
write.csv(cor_all_pwys,'TP53mut_cor_pwys.csv')
TP53m_df_m=driver_mut_df[TP53_mut_pdxs,]
TP53m_driver_mut_df=as.data.frame(cbind('Responder'=response[TP53_mut_pdxs,'Responder'],TP53m_df_m))
ae=apply(TP53m_driver_mut_df,2,function(x) length(unique(x))==1)
TP53m_driver_mut_df_nae=TP53m_driver_mut_df[,!ae]
TP53m_driver_mut_pval=sapply(TP53m_driver_mut_df_nae[,-1],function(x) chisq.test(x,TP53m_driver_mut_df_nae$Responder)$p.value)
TP53m_driver_mut_pval=sort(TP53m_driver_mut_pval)
TP53m_res_driver_mut_df=as.data.frame(sapply(names(TP53m_driver_mut_pval),function(x) c(TP53m_driver_mut_pval[x],
                                                                                        as.vector(table(TP53m_driver_mut_df$Responder,TP53m_driver_mut_df[,x])))))
write.csv(t(TP53m_res_driver_mut_df),'TP53mut_cor_genes_mut.csv')
TP53m_pwy_mut_df=pw_mut_df[,TP53_mut_pdxs]
ae=apply(TP53m_pwy_mut_df,1,function(x) length(unique(x))==1)
TP53m_pwy_mut_df=TP53m_pwy_mut_df[!ae,]
t.val=sort(apply(TP53m_pwy_mut_df,1,function(x) t.test(x~factor(TP53m_res))$p.value))
pwy='GO_CENTROSOME'
table(TP53m_pwy_mut_df[pwy,],TP53m_res)
intersect(geneset[[pwy]]@geneIds,t1_genes)
pwys_mut_eff=t(sapply(names(t.val[t.val<0.01]),function(x) c(t.val[x],paste(intersect(geneset[[x]]@geneIds,t1_genes),collapse=","))))
colnames(pwys_mut_eff)=c('p.val','Cancer Driver Gene')
write.csv(pwys_mut_eff,'TP53mut_response_mutated_pwys.csv')


TP53m_res_g=as.data.frame(cbind('Responder'=TP53m_res,t(TP53m_df_g)))
TP53m_genes=feature_sel(TP53m_res_g,subsample = 1)
res_sel=TP53m_res_g[,c('Responder',TP53m_genes$sel_genes,'RDX','IDO1')]
featurePlot(x =res_sel[,-1],y =res_sel[,1], between = list(x = 1, y = 1),type = c("g", "p"))
table(TP53m_res_g$IDO1<0.1,TP53m_res_g$Responder) ##IDO1 seems to be a good predictor


TP53m_res_p=as.data.frame(cbind('Responder'=TP53m_res,t(TP53m_df_p)))
TP53m_pwys=feature_sel(TP53m_res_p,nrep=10)
res_sel=TP53m_res_p[,c('Responder',TP53m_pwys$sel_genes)]
featurePlot(x =res_sel[,-1],y =res_sel[,1], between = list(x = 1, y = 1),type = c("g", "p"))
table(TP53m_res_p$PID_VEGFR1_PATHWAY>0,TP53m_res_g$Responder) 
table(TP53m_res_p$REACTOME_AXON_GUIDANCE>0,TP53m_res_g$Responder) 
table(TP53m_res_p$PID_S1P_S1P1_PATHWAY>0,TP53m_res_g$Responder)
table(TP53m_res_p$KEGG_RIBOSOME>0,TP53m_res_g$Responder)

###DEG analysisREACTOME_AXON_GUIDANCE
deg=function(x,f){
  t=t.test(x~f)
  w=wilcox.test(x~f)
  c(t$estimate,t$statistic,'pval'=t$p.value,'pval_mwu'=w$p.value)
}
class=factor(TP53m_res)
#res_pdxs=response[response$response %in% c('CR','CPD'),'PDX']
#class=factor(response[res_pdxs,'Responder'])
Cohort_degTable=apply(TP53m_df_g,1,function(x) deg(x,class))
Cohort_degTable=as.data.frame(t(Cohort_degTable))
Cohort_degTable$diff=Cohort_degTable[,2]-Cohort_degTable[,1] #difference of log2(RSEM+1)
Cohort_degTable$gs=rownames(Cohort_degTable)
Cohort_degTable$adj_p = p.adjust(Cohort_degTable$pval_mwu,method='fdr')
write.csv(Cohort_degTable,'TP53mut_deg.csv')
Cohort_Pos_high=Cohort_degTable %>% filter(diff>0 & pval_mwu<0.05)
Cohort_Pos_low=Cohort_degTable %>% filter(diff<0 & pval_mwu<0.05)
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
gs2Reactome(Cohort_pos_up_gs,paste0('TP53mut_DEG','_POS_UP'),pval=0.05)
Cohort_pos_down_gs=Cohort_Pos_low$gs
gs2Reactome(Cohort_pos_down_gs,paste0('TP53mut_DEG','_POS_DOWN'))

Cohort_degTable=apply(TP53m_df_p,1,function(x) deg(x,class))
Cohort_degTable=as.data.frame(t(Cohort_degTable))
Cohort_degTable$diff=Cohort_degTable[,2]-Cohort_degTable[,1] #difference of log2(RSEM+1)
Cohort_degTable$adj_p = p.adjust(Cohort_degTable$pval_mwu,method='fdr')
write.csv(Cohort_degTable,'TP53mut_deg_pwys.csv')

#########
TP53m_sel_df=sel_df[TP53_mut_pdxs,c(1,3,4,5)]
TP53m_sel_df$IDO1=t(log_filteredExpr3)[TP53_mut_pdxs,'IDO1']

sel_pwys=c('KEGG_RIBOSOME','REACTOME_TRANSLATION',
           'REACTOME_ADAPTIVE_IMMUNE_SYSTEM','BIOCARTA_THELPER_PATHWAY','BIOCARTA_TCYTOTOXIC_PATHWAY',
           'PID_S1P_S1P1_PATHWAY','PID_ERBB1_RECEPTOR_PROXIMAL_PATHWAY','PID_VEGFR1_PATHWAY','BIOCARTA_HER2_PATHWAY','KEGG_MTOR_SIGNALING_PATHWAY','WIKI_._Angiogenesis_overview')
sel_pwys=TP53m_pwys$sel_genes
TP53m_pas_df=t(TP53m_df_p[sel_pwys,])
CT=factor(substr(rownames(TP53m_pas_df),1,2))

annotation_col = data.frame(
  Cancer_type = as.character(CT),
  Response=TP53m_sel_df$Responder,
  Hypermutated=ifelse(TP53m_sel_df$Non_TP53_mut_sum>=2,'Yes','No'),
  IDO1_expression=ifelse(TP53m_sel_df$IDO1>0.1,'Yes','No')
)
ann_colors = list(
  Response = c('Non_Responder'='firebrick3','Responder'='green'),
  Cancer_type=c('CR'='gray','LU'='cyan','PA'='yellow','OV'='purple'),
  Hypermutated=c('Yes'='red','No'='blue'),
  IDO1_expression=c('Yes'='red','No'='blue')
)
rownames(annotation_col)=rownames(TP53m_pas_df)
#TP53m_pas_df=t(apply(TP53m_pas_df,1,min_max_norm))
pheatmap(t(TP53m_pas_df), show_colnames = T,
         clustering_distance_rows = 'minkowski',
         clustering_distance_cols = 'minkowski',
         cutree_cols = 2,
         annotation_col = annotation_col,annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.copy(png,paste0('TP53mut_Response','_mut_pwy.png'),width=1500,height=800, res=150)
dev.off()

TP53m_sel_df_binary=TP53m_sel_df[,c(1,4,5)]
names(TP53m_sel_df_binary)[2:3]=c('Hypermutated','IDO1_expr')
TP53m_sel_df_binary$Hypermutated=ifelse(TP53m_sel_df_binary$Hypermutated>=2,1,0)
TP53m_sel_df_binary$IDO1_expr=ifelse(TP53m_sel_df_binary$IDO1_expr>=0.1,1,0)
TP53m_sel_df_binary=cbind(TP53m_sel_df_binary,TP53m_pas_df)
TP53m_sel_df_binary=TP53m_sel_df_binary[,-2]
require(rpart)
fit=rpart(Responder~.,data=TP53m_sel_df_binary,method='class',control = rpart.control(cp = 0.01))
require(rattle)
fancyRpartPlot(fit)
TP53m_fs=feature_sel(TP53m_sel_df_binary,nrep=10,min_f = 2)
TP53m_sel=TP53m_sel_df_binary[,c('Responder',TP53m_fs$sel_genes)]
featurePlot(x =TP53m_sel[,c(-1,-2)],y =as.factor(TP53m_sel[,1]), between = list(x = 1, y = 1),type = c("g", "p"))

cf_models_TP53m=trainModel(TP53m_sel_df_binary,nrep=10,fs=FALSE,ROC=TRUE)


#################deprecated################################
TP53_mut=factor(driver_mut_df[rownames(res_pwys_df),'TP53'])
res_pwys_df$TP53_mut=TP53_mut
ggplot(res_pwys_df)+aes(Responder,G1_S_CHECKPOINT)+
  geom_boxplot(notch=T)+geom_jitter(width=0.1,aes(col=TP53_mut))+
  facet_wrap(~substr(rownames(res_pwys_df),1,2))+mytheme
ggplot(res_pwys_df)+aes(Responder,REACTOME_CELL_CYCLE)+
  geom_boxplot(notch=T)+geom_jitter(width=0.1,aes(col=TP53_mut))+
  facet_wrap(~substr(rownames(res_pwys_df),1,2))+mytheme
ggplot(res_pwys_df)+aes(Responder,WNT_CANONICAL)+
  geom_boxplot(notch=T)+geom_jitter(width=0.1,aes(col=TP53_mut))+
  facet_wrap(~substr(rownames(res_pwys_df),1,2))+mytheme


g1s_ck_genes=acsn_pwys[['G1_S_CHECKPOINT']]
g1s_ck_ex=log_filteredExpr3[g1s_ck_genes,rownames(res_driver_mut_df)]
g1s_ck_ex=t(apply(g1s_ck_ex,1,function(v) (v-min(v))/(max(v)-min(v))))
pheatmap(rbind(sel_genes_mut['TP53',],g1s_ck_ex), show_colnames = T,
         clustering_distance_rows = 'minkowski',
         clustering_distance_cols = 'minkowski',
         cutree_cols = 3,
         annotation_col = annotation_col,annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))


