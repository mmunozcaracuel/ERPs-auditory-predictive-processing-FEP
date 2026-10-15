# Clinical and ERP group comparisons adjusted for age, sex and education.
group_rows<-list();int_rows<-list();clinical<-list()
{
 spec<-'Main';cov<-c('AGE_c','SEX','EDU_c')
 for(fam in names(families)[1:4])for(y in families[[fam]]){
  k<-paste('Clinical',spec,y,sep='__')
  clinical[[k]]<-cbind(Spec=spec,Family=fam,Outcome=y,lm_result(u,NULL,y,cov,k,TRUE))
 }
 for(m in measures){
  a<-p;a$Y<-a[[m]];a$Electrode<-factor(a$Electrode,levels=sites)
  f<-lmer(as.formula(paste('Y~GROUP*Electrode+',paste(cov,collapse='+'),'+(1|Subject)')),a,REML=TRUE)
  k<-paste('LMM',spec,m,sep='__');models[[k]]<-f
  e<-emmeans(f,~GROUP|Electrode,weights='equal',lmer.df='satterthwaite')
  z<-as.data.frame(summary(contrast(e,list(FEP_minus_HC=c(-1,1))),infer=c(TRUE,TRUE),adjust='none'))
  s<-vapply(as.character(z$Electrode),function(site)pool_sd(a$Y[a$Electrode==site],a$GROUP[a$Electrode==site]),0)
  group_rows[[k]]<-data.frame(Spec=spec,Tier=tier(m),ERP=m,Site=as.character(z$Electrode),Scope='Five_electrodes',n=length(unique(a$Subject)),
   b=z$estimate,SE=z$SE,df=z$df,beta=z$estimate/s,low=z$lower.CL/s,high=z$upper.CL/s,p=z$p.value)
  an<-anova(f,ddf='Satterthwaite')['GROUP:Electrode',]
  int_rows[[k]]<-data.frame(Spec=spec,Tier=tier(m),ERP=m,F=an$`F value`,df1=an$NumDF,df2=an$DenDF,p=an$`Pr(>F)`)
  k<-paste('ROI',spec,m,sep='__')
  group_rows[[k]]<-cbind(Spec=spec,Tier=tier(m),ERP=m,Site='ROI',Scope='ROI',lm_result(roi,NULL,m,cov,k,TRUE))
 }
}
gr<-bh(do.call(rbind,group_rows),c('Spec','Tier','Scope'))
gi<-bh(do.call(rbind,int_rows),c('Spec','Tier'))
cl<-bh(do.call(rbind,clinical),c('Spec','Family'))
cat('Group contrasts completed\n');print(subset(gr,Site%in%c('FCz','ROI')&ERP%in%components)[,c('Spec','ERP','Site','beta','p','q')])

# Separate ERP associations and ERP-by-group interaction models.
assoc<-list();slope_int<-list()
for(site in c(sites,'ROI'))for(fam in names(families))for(y in families[[fam]])for(m in components){
 spec<-'Main';patient<-fam%in%c('PANSS','SFS_total','SFS_subscales')
 a<-d[d$Electrode==site,];if(patient)a<-a[a$GROUP=='FEP',]
 cov<-c(if(!patient)'GROUP','AGE_c','SEX','EDU_c')
 k<-paste('Association',spec,site,fam,m,y,sep='__')
 meta<-data.frame(Spec=spec,Tier=tier(m),Family=fam,Sample=if(patient)'FEP' else 'Pooled',Site=site,Scope=if(site=='ROI')'ROI' else 'Five_electrodes',ERP=m,Outcome=y)
 assoc[[k]]<-cbind(meta,lm_result(a,m,y,cov,k))
 if(!patient){
  k<-sub('Association','ERP_group',k)
  slope_int[[k]]<-cbind(meta,lm_result(a,m,y,setdiff(cov,'GROUP'),k,interaction=TRUE))
 }
}
assoc<-bh(do.call(rbind,assoc),c('Spec','Tier','Family','Scope'))
si<-bh(do.call(rbind,slope_int),c('Spec','Tier','Family','Scope'))
assoc<-local_fdr(assoc);si<-local_fdr(si)

# Medication analyses reported in Table S6; dose is log1p transformed.
# No dose-phenotype associations, no extreme quartiles.
med<-list()
seven<-c(families$Main,families$PANSS,families$SFS_total)
{
 sample<-'Medicated31'
 spec<-'Main';a0<-d[d$GROUP=='FEP' & d$CPZ_EQ>0,]
 cov<-c('AGE_c','SEX','EDU_c')
 for(site in c(sites,'ROI'))for(m in components){
  a<-a0[a0$Electrode==site,];k<-paste('Dose_ERP',spec,sample,site,m,sep='__')
  med[[k]]<-cbind(Spec=spec,Sample=sample,Analysis='Dose_ERP',Tier=tier(m),Family='Dose_ERP',Site=site,Scope=if(site=='ROI')'ROI' else 'Five_electrodes',ERP=m,Outcome=m,lm_result(a,'log_CPZ',m,cov,k))
 }
 for(site in c(sites,'ROI'))for(m in components){
  a<-a0[a0$Electrode==site,]
  fs<-c(if(site=='FCz')list(Seven_outcomes=seven),families[c('PANSS','SFS_total','SFS_subscales')])
  for(fam in names(fs))for(y in fs[[fam]]){
   k<-paste('ERP_outcome_CPZ',spec,sample,site,m,fam,y,sep='__')
   med[[k]]<-cbind(Spec=spec,Sample=sample,Analysis='ERP_outcome_CPZ',Tier=tier(m),Family=fam,Site=site,Scope=if(fam=='Seven_outcomes')'FCz' else if(site=='ROI')'ROI' else 'Five_electrodes',ERP=m,Outcome=y,lm_result(a,m,y,c(cov,'log_CPZ'),k))
   no_dose<-lm_result(a,m,y,cov,paste0(k,'__no_dose'))
   med[[k]]$beta_no_dose<-no_dose$beta;med[[k]]$low_no_dose<-no_dose$low;med[[k]]$high_no_dose<-no_dose$high;med[[k]]$p_no_dose<-no_dose$p
   med[[k]]$delta_beta_dose<-med[[k]]$beta-no_dose$beta
  }
 }
}
med<-bh(bind_fill(med),c('Spec','Sample','Analysis','Tier','Family','Scope'))
med$q_within_site<-bh(med,c('Spec','Sample','Analysis','Tier','Family','Site'),qcol='q_local')$q_local

