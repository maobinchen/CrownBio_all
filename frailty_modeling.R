require(survival)
require(survminer)
require(frailtypack)
genSurvData=function(df){
  cutoff=log(3) #triple volume size as endpoint
  df$diff=df$logTV-df[1,'logTV']
  PDX=df[1,'PDX']
  treatment=as.character(df[1,'Treatment'])
  Mouse=df[1,"Mouse"]
  ind = which(df$diff>=cutoff)
  if(length(ind) > 0){
    i=ind[1]
    tri_t = df[i-1,'Day']+(df[i,'Day']-df[i-1,'Day'])*(cutoff+df[1,'logTV']-df[i-1,'logTV'])/(df[i,'logTV']-df[i-1,'logTV']) #use linear interpolation at logTV scale to estimate triple time
    c("PDX"=PDX,"Treatment"=treatment,"Mouse"=Mouse,"time"=tri_t,"censoring"=1)
  }else{
    c("PDX"=PDX,"Treatment"=treatment,"Mouse"=Mouse,"time"=df[nrow(df),'Day'],"censoring"=0)
  }
}
byMouse=by(dat,dat$Mouse,genSurvData)
surv=as.data.frame(do.call(rbind,byMouse))
surv$time=as.numeric(as.character(surv$time))
surv <- within(surv, Treatment <- relevel(Treatment, ref = 'Vehicle'))
surv$censoring=as.numeric(as.character(surv$censoring))
surv$CT=substr(surv$PDX,1,2)
surv$Treatment=ifelse(surv$Treatment=='Vehicle',0,1)
frailty<-frailtyPenal(Surv(time,censoring)~cluster(PDX)+Treatment,data=surv,n.knots = 7,kappa=1000)
frailty
addFrail<-additivePenal(Surv(time,censoring)~cluster(PDX)+Treatment+slope(Treatment),data=surv,
                        correlation = T,n.knots = 8,kappa=10000)
addFrail

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
surv_trt=subset(surv,Treatment=='Treatment')
s2=with(surv_trt,Surv(time,censoring))
fit2=survfit(s2~CT, data=surv_trt)
ggsurvplot(fit2,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw())
frail2<-frailtyPenal(Surv(time,censoring)~cluster(PDX)+CT,data=surv,n.knots = 7,kappa=1000)

genes=c('TBX3','HIST2H4A','UBC')
gene_pdx_df=cbind(PDX=rownames(log_ex3),log_ex3[,genes])
surv_gene_df=merge(surv,gene_pdx_df,by='PDX')
#coxph(s1~Treatment*HIST2H4A,data=surv_gene_df)
addFrail<-additivePenal(Surv(time,censoring)~cluster(PDX)+Treatment*TBX3+slope(Treatment),data=surv_gene_df,
                        correlation = T,n.knots = 8,kappa=1000)
