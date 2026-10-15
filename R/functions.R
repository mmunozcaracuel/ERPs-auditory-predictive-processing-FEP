sites <- c('Fz','FCz','Cz','FC3','FC4')
components <- c('CNV','MMN','PINV','N1','P2')
tier <- function(m)ifelse(m%in%c('CNV','MMN','PINV'),'Primary',ifelse(m%in%c('N1','P2'),'N1_P2','Other_exploratory'))
families <- list(Main=c('CAPE_TOT','SCIP_TOT','MASC_COR'),CAPE=c('CAPE_P','CAPE_N','CAPE_D'),
 SCIP=c('SCIP_AV_I','SCIP_MT','SCIP_FV','SCIP_AV_D','SCIP_VP'),MASC=c('MASC_NO','MASC_HIPO','MASC_HIPER'),
 PANSS=c('PANSS_P','PANSS_N','PANSS_PG'),SFS_total='SFS_TOT_T',
 SFS_subscales=c('SFS_AIS','SFS_COM','SFS_AE','SFS_AC','SFS_O','SFS_PRO','SFS_E'))

measures <- c(components, 'CNV_late','PINV_early','MMN_peak','MMN_latency','PINV_peak','PINV_latency')

# Validate participant-level invariance and a complete five-electrode design.
prepare_data <- function(input){
 a <- read.csv(input,check.names=FALSE,na.strings=c('','NA','NaN'))
 covariates<-c('GROUP','AGE','SEX','YEARS_EDUCATION','CPZ_EQ',unlist(families))
 need<-unique(c('Subject','Electrode',covariates,measures))
 absent<-setdiff(need,names(a))
 if(length(absent))stop('Missing columns: ',paste(absent,collapse=', '))
 a<-a[,need];a$Subject<-as.character(a$Subject)
 if(anyNA(a[c('Subject','Electrode','GROUP','SEX')]))stop('Missing participant, electrode or categorical values.')
 if(anyDuplicated(a[c('Subject','Electrode')]))stop('Duplicate participant-electrode keys.')
 if(!setequal(unique(a$Electrode),sites))stop('Input must contain only the five frontocentral electrodes; ROI is computed internally.')
 if(any(table(a$Subject)!=5))stop('Each participant must have all five electrodes.')
 if(!setequal(unique(a$GROUP),c('CON','FEP')))stop('GROUP must use CON and FEP.')
 if(length(unique(a$SEX))!=2L)stop('SEX must have two consistently coded levels.')
 for(v in covariates)if(any(vapply(split(a[[v]],a$Subject),function(x)length(unique(x)),0)!=1))stop('Participant-level variable varies across electrodes: ',v)
 pooled_outcomes<-unlist(families[1:4]);patient_outcomes<-unlist(families[5:7])
 numeric_all<-c('AGE','YEARS_EDUCATION',pooled_outcomes,measures)
 for(v in numeric_all)if(!is.numeric(a[[v]])||any(!is.finite(a[[v]])))stop('Expected complete finite numeric column: ',v)
 for(v in c('CPZ_EQ',patient_outcomes))if(!is.numeric(a[[v]])||any(!is.finite(a[a$GROUP=='FEP',v])))stop('Missing/non-numeric FEP measure: ',v)
 if(any(a$CPZ_EQ[a$GROUP=='FEP']<0))stop('Negative antipsychotic dose.')
 u<-a[!duplicated(a$Subject),]
 if(nrow(u)!=68L||any(table(u$GROUP)!=34L))warning('Sample differs from manuscript (34 per group); results will not reproduce the published sample.')
 a$GROUP<-factor(a$GROUP,levels=c('CON','FEP'));a$SEX<-factor(a$SEX);a$Subject<-factor(a$Subject)
 a$AGE_c<-a$AGE-mean(u$AGE);a$EDU_c<-a$YEARS_EDUCATION-mean(u$YEARS_EDUCATION)
 a$log_CPZ<-log1p(a$CPZ_EQ)
 roi<-a[a$Electrode=='FCz',]
 for(m in measures)roi[[m]]<-vapply(as.character(roi$Subject),function(s)mean(a[a$Subject==s,m]),0)
 roi$Electrode<-'ROI'
 rbind(a,roi)
}

bind_fill <- function(rows){
 fields<-unique(unlist(lapply(rows,names)))
 do.call(rbind,lapply(rows,function(x){for(nm in setdiff(fields,names(x)))x[[nm]]<-NA_real_;x[,fields]}))
}

# Estimation, standardization and Benjamini-Hochberg correction.
pool_sd<-function(y,g){z<-split(y,g);sqrt(sum(vapply(z,function(x)(length(x)-1)*var(x),0))/sum(lengths(z)-1))}
models<-list()
lm_result<-function(a,x,y,cov,key,group=FALSE,interaction=FALSE){
 a$Y<-a[[y]];if(!group)a$X<-a[[x]]
 needed<-c('Y',if(!group)'X',cov,if(group||interaction)'GROUP')
 a<-droplevels(a[complete.cases(a[,needed]),])
 f<-lm(reformulate(c(if(group)'GROUP' else if(interaction)'X*GROUP' else 'X',cov),'Y'),a)
 term<-if(group)'GROUPFEP' else if(interaction)'X:GROUPFEP' else 'X'
 M<-model.matrix(f);stopifnot(f$rank==ncol(M),nrow(a)>ncol(M)+5,!anyDuplicated(a$Subject))
 z<-coef(summary(f))[term,];ci<-confint(f,term)
 scale<-if(group)1/pool_sd(a$Y,a$GROUP) else sd(a$X)/sd(a$Y)
 data.frame(n=nrow(a),b=unname(z[1]),SE=unname(z[2]),df=df.residual(f),beta=unname(z[1])*scale,low=ci[1]*scale,high=ci[2]*scale,
 p=unname(z[4]))
}
bh<-function(z,keys,pcol='p',qcol='q'){
 key<-do.call(paste,c(z[keys],sep='|'));z[[qcol]]<-NA_real_;z$family_n<-NA_integer_
 for(k in unique(key)){j<-key==k;z[[qcol]][j]<-p.adjust(z[[pcol]][j],'BH');z$family_n[j]<-sum(j)};z
}

local_fdr<-function(x){z<-bh(x,c('Spec','Tier','Family','Site'),qcol='q_within_site');x$q_within_site<-z$q_within_site;x$local_family_n<-z$family_n;x}
