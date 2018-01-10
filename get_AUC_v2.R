get_AUC=function(data,CI=FALSE){
	result = list()
	PDX_tr= with(data,by(data,Mouse,getAUC))
	AUC=do.call(rbind,PDX_tr)
	x=do.call(rbind,by(AUC,AUC$PDX,getMedianAUC))
	#auc_r=setNames(as.vector(x),names(x))
	result=list('AUC'=AUC,'AUC_r'=x[,1],'RTV'=x[,2],'wtmc'=x[,3])
	if(CI){
		AUC_simu=by(AUC,AUC$PDX,simulate_CI)
		AUC_r_CI=do.call(rbind,AUC_simu)
		AUC_r_CI$AUC=auc_r[rownames(AUC_r_CI)]
		AUC_r_CI$width=AUC_r_CI$Upper_bond-AUC_r_CI$Lower_bond
		result$AUC_r_CI=AUC_r_CI
	}
	result	
}

#from the last observation, finding the longest path with consistent slope
findPath=function(v){
	len_v=length(v)
	d=v[-1]-v[1:(len_v-1)]
	d=ifelse(d>0,1,-1)
	ladder=cumsum(rev(d))
	ind=which.max(abs(ladder))
	start=1
	if(max(abs(ladder))>=floor(len_v/2)) start=len_v-ind #the starting position for the longest consistent path; at least longer than half of the whole path
	start
}

#modified AUC (deal with V-shape(regrowth) or inverse V-shape(delayed kick in), goal: from the last observation, finding the longest path with consistent slope)
getAUC_old=function(df){
	df=as.data.frame(df)
	n=nrow(df)
	day_v=df$Day
	TV_v=log(df$TV+0.001)
	start=1
	if(length(TV_v)>3){
		start=findPath(TV_v)
	}
	day_v=day_v[start:n]
	TV_v=TV_v[start:n]
	n1=length(day_v)
	#TV_v=df$logTV
	AUC=0
	for(i in 2:n1){
		AUC = AUC+(day_v[i]-day_v[i-1])*(TV_v[i]+TV_v[i-1])/2
	}
	#print(AUC)
	day_span=day_v[n1]-day_v[1]
	AUC = AUC-TV_v[1]*day_span
	AUC = 2*AUC/(day_span*day_span) #tumor growth rate under exponential growth model
	data.frame('PDX'=df$PDX[1],'Treatment'=df$Treatment[1],'Mouse'=df$Mouse[1],'AUC'=AUC)
}
#calculate tumor growth rate on every time span, then get the weighted average tumor growth rate(assign more weights on latest observatons)
getAUC=function(df){
	df=as.data.frame(df)
	n=nrow(df)
	day_v=df$Day
	TV_v=log(df$TV+0.001)
	start=1
	if(length(TV_v)>3){
		start=findPath(TV_v)
	}
	day_v=day_v[start:n]
	TV_v=TV_v[start:n]
	n1=length(day_v)
	#TV_v=df$logTV
	AUC=0
	for(i in 2:n1){
		AUC = AUC+(day_v[i]-day_v[i-1])*(TV_v[i]+TV_v[i-1])/2
	}
	#print(AUC)
	day_span=day_v[n1]-day_v[1]
	AUC = AUC-TV_v[1]*day_span
	AUC = 2*AUC/(day_span*day_span) #tumor growth rate under exponential growth model
	data.frame('PDX'=df$PDX[1],'Treatment'=df$Treatment[1],'Mouse'=df$Mouse[1],'AUC'=AUC)
}

#calculate day 21 RTV and TVR, also calculate weight*K(t)-K(c), weight is the relative importantce we assigned to the treatment, adjusted by control
getMedianAUC=function(df){
	weight=2
	vehicle_v=df[df$Treatment=="Vehicle",]$AUC
	treatment_v=df[df$Treatment=="Treatment",]$AUC
	cbs=expand.grid(treatment_v,vehicle_v)
	AUC_r = cbs[,1]-cbs[,2]
	AUC_r=AUC_r[!(is.na(AUC_r) | is.nan(AUC_r))]
	d21_rtv=exp(mean(AUC_r)*21) #tumor size reduction relative to vehicle,normalize to 21 days 
	d21_tvr=exp(mean(treatment_v,na.rm=T)*21) #treatment tumor size ratio to inital tumor size
	wtmc=weight*mean(treatment_v,na.rm=T)-mean(vehicle_v,na.rm=T) #weighted kt minus kc
	c(d21_rtv,d21_tvr,wtmc)
}

#get simulated Confidence interval for medianAUCratio
simulate_CI=function(df,nrep=1000,CI=95){
	vehicle_v=df[df$Treatment=="Vehicle",]$AUC
	n_v=length(vehicle_v)
	treatment_v=df[df$Treatment=="Treatment",]$AUC
	n_t=length(treatment_v)
	median_auc_r_v = rep(NA,nrep)
	for(i in 1:nrep){
		v1=sample(vehicle_v,n_v,replace=T)
		v2=sample(treatment_v,n_t,replace=T)
		cbs=expand.grid(v2,v1)
		AUC_r = cbs[,1]-cbs[,2]
		AUC_r=AUC_r[!(is.na(AUC_r) | is.nan(AUC_r))]
		median_auc_r_v[i]=exp(median(AUC_r)*21)
	}
	data.frame('mean'=mean(median_auc_r_v),'median'=median(median_auc_r_v),'Lower_bond'=quantile(median_auc_r_v,0.05),'Upper_bond'=quantile(median_auc_r_v,0.95))
}