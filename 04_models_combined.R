# Combined 04 model plotting script for full_manuscript.Rmd
# Source scripts: 04_models.R, 04_models_beta_select.R, 04_models_cross_validate.R
# This file is scoped to manuscript-included outputs only.

# ============================================================
# Section 1: DHARMa diagnostic plots used by the manuscript
# Outputs: plots/dharma_m*.png
# Original source: 04_models.R
# ============================================================
library(reshape2)
library(ggplot2)
library(mgcv)  
library(dplyr)
library(ggridges)
library(png)
library(PBSmodelling)
library(patchwork)
library(mgcv)
library(png)
library(grid)
library(DHARMa)
library(itsadug)
library(mgcViz)
library(directlabels)
outs<-read.csv("data/all_output.csv")
outs<-filter(outs,Year<2023)
alt_met<-read.csv("data/alt_metrics_calc.csv")
unc_mort<-read.csv("data/uncertainty_mort.csv")

#======================================
# full models:
# recruitment <- density + temperature + competitor + size
#==show two things
#=====1. some of the processes covary
#=====2. some of the processes can be explained by environmental variables

#==two figures: one for recruitment, one for mortality
#==first column are the estimates and the model fits
#==remaining columns represent the variables
#==only 'significant' variables are colored in
#==CHECK ON APPROPRIATE LAGS FOR RECRUITMENT
#==THINK ABOUT LARGE SCALE DRIVERS TOO

#==king crabs
use_stocks<-c("BBRKC","PIRKC","SMBKC","PIBKC")
rec_term<-NULL
mort_term<-NULL
out_plot_r<-NULL
out_plot_m<-NULL
dev_expl_m<-NULL
keep_AIC_r<-NULL
keep_AIC_m<-NULL
keep_AIC_m_big<-NULL
keep_conc<-list(list())
conc_cnt<-1
for(y in 1:length(use_stocks))
{
set1<-filter(outs,species==use_stocks[y])[,-c(1,6)]
set2<-filter(alt_met,stock==use_stocks[y])[,-1]
set3<-filter(unc_mort)
colnames(set2)[4]<-"species"

casted<-dcast(set1,Year~process,value.var="values")
colnames(casted)[4]<-"Other_mortality"
mod_dat<-merge(casted,set2,by="Year")
if(use_stocks[y]=="SMBKC")
  mod_dat<-mod_dat[-1,]
mod_dat$lag_temp<-c(NA,mod_dat$Temperature[-length(mod_dat$Temperature)])

#==recruits
base_r <- gam(Recruitment ~ 1, data = mod_dat,family=nb(link = "log"))
summary(base_r)
mod_r<-gam(data=mod_dat,Recruitment~s(Abundance,k=4)+s(Temperature,k=4),family=nb(link = "log"))
summary(mod_r)
plotted<-plot(mod_r,pages=1)

keep_AIC_r<-rbind(keep_AIC_r,c(AIC(base_r),AIC(mod_r)))

#ccf(mod_dat$Recruitment,mod_dat$Temperature,na.action=na.pass)
preds<-predict.gam(mod_r,type='response',se.fit=TRUE)

plo_gam_r<-data.frame(obs=mod_dat$Recruitment[!is.na(mod_dat$Recruitment)],
                      preds=preds$fit,
                      y_up=preds$fit+preds$se,
                      y_dn=preds$fit-preds$se,
                      year=mod_dat$Year[!is.na(mod_dat$Recruitment)],
                      stock=rep(use_stocks[y],length(preds$fit)))
out_plot_r<-rbind(out_plot_r,plo_gam_r)

for(x in 1:(length(plotted)))
{
  temp<-data.frame(x=plotted[[x]]$x,
                   y=plotted[[x]]$fit,
                   y_up=plotted[[x]]$fit+plotted[[x]]$se,
                   y_dn=plotted[[x]]$fit-plotted[[x]]$se,
                   covar=plotted[[x]]$xlab,
                   stock=use_stocks[y])
  rec_term<-rbind(rec_term,temp)
}


#==mortality
base_mod_m<-gam(data=mod_dat,Other_mortality~1,family=tw())
mod_m<-gam(data=mod_dat,Other_mortality~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3),family=tw())
mod_m_a<-gam(data=mod_dat,Other_mortality~s(Abundance,k=4),family=tw())
mod_m_t<-gam(data=mod_dat,Other_mortality~s(Temperature,k=4),family=tw())
mod_m_s<-gam(data=mod_dat,Other_mortality~s(Size,k=3),family=tw())

 keep_conc[[conc_cnt]]<- concurvity(mod_m,full=FALSE)$observed[-1,-1]
   conc_cnt<-conc_cnt+1
# gam.check(mod_m)
# acf(resid(mod_m))
# pacf(resid(mod_m))
simout<-simulateResiduals(mod_m,n=250)

png(paste("plots/dharma_m",use_stocks[y],".png",sep=''),height=6,width=8,res=350,units='in') 
plot(simout)
dev.off()


summary(mod_m)
summary(base_mod_m)
plot(mod_m,pages=1)
dev_expl_m<-c(dev_expl_m,round(summary(mod_m)$dev,2))
keep_AIC_m<-rbind(keep_AIC_m,c(AIC(base_mod_m),AIC(mod_m)))
keep_AIC_m_big<-rbind(keep_AIC_m_big,c(AIC(base_mod_m),AIC(mod_m),AIC(mod_m_a),AIC(mod_m_t),AIC(mod_m_s)))

mod_m1<-gam(data=mod_dat,Other_mortality~s(Abundance,k=4)+s(lag_temp,k=4)+s(Size,k=3),family=tw())
summary(mod_m1)
plotted<-plot(mod_m,pages=1)

# ccf(mod_dat$Other_mortality,mod_dat$Temperature)
# ccf(mod_dat$Other_mortality,mod_dat$Size)

preds<-predict.gam(mod_m,type='response',se.fit=TRUE)

plo_gam_m<-data.frame(obs=mod_dat$Other_mortality[!is.na(mod_dat$Other_mortality)],
                      preds=preds$fit,
                      y_up=preds$fit+preds$se,
                      y_dn=preds$fit-preds$se,
                      year=mod_dat$Year[!is.na(mod_dat$Other_mortality)],
                      stock=use_stocks[y])

out_plot_m<-rbind(out_plot_m,plo_gam_m)


for(x in 1:(length(plotted)))
{
  temp<-data.frame(x=plotted[[x]]$x,
                   y=plotted[[x]]$fit,
                   y_up=plotted[[x]]$fit+plotted[[x]]$se,
                   y_dn=plotted[[x]]$fit-plotted[[x]]$se,
                   covar=plotted[[x]]$xlab,
                   stock=use_stocks[y])
  mort_term<-rbind(mort_term,temp)
}

}

#==do a Pribs model
set1_rkc<-filter(outs,species%in%c("PIRKC"))[,-c(1,6)]
set2_rkc<-filter(alt_met,stock%in%c("PIRKC"))[,-1]
colnames(set2)[4]<-"species"

casted<-dcast(set1_rkc,Year~process,value.var="values")
colnames(casted)[4]<-"Other_mortality"
mod_dat_rkc<-merge(casted,set2_rkc,by="Year")

set1_bkc<-filter(outs,species%in%c("PIBKC"))[,-c(1,6)]
set2_bkc<-filter(alt_met,stock%in%c("PIBKC"))[,-1]
colnames(set2)[4]<-"species"

casted<-dcast(set1_bkc,Year~process,value.var="values")
colnames(casted)[4]<-"Other_mortality"
mod_dat_bkc<-merge(casted,set2_bkc,by="Year")

#==predict pirkc with pibkc abundance 
mod_dat<-mod_dat_rkc
mod_dat$bkc_n<-mod_dat_bkc$Abundance[-1]
mod_dat$bkc_ssb<-mod_dat_bkc$`Spawner abundance`[-1]
mod_m<-gam(data=mod_dat,Other_mortality~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3)+s(bkc_ssb,k=3)+s(bkc_n,k=3),family=tw())
summary(mod_m)
plot(mod_m,pages=1)

#==predict pibkc M with pirkc abundance 
mod_dat<-mod_dat_bkc[-1,]
mod_dat$rkc_n<-mod_dat_rkc$Abundance
mod_dat$rkc_ssb<-mod_dat_rkc$`Spawner abundance`
mod_m<-gam(data=mod_dat,Other_mortality~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3),family=tw())
summary(mod_m)
plot(mod_m,pages=1)

#==predict pibkc R with pirkc abundance 
mod_dat<-mod_dat_bkc[-1,]
mod_dat$rkc_n<-mod_dat_rkc$Abundance
mod_dat$rkc_ssb<-mod_dat_rkc$`Spawner abundance`
mod_m<-gam(data=mod_dat,Recruitment ~s(Abundance,k=3)+s(Temperature,k=3)+s(Size,k=3)+s(rkc_n,k=3),family=tw())
summary(mod_m)
plot(mod_m,pages=1)

#============================
# chionoecetes species
#============================
use_stocks<-c("Snow","Tanner")
for(y in 1:length(use_stocks))
{
  set1<-filter(outs,species==use_stocks[y])[,-c(1,6)]
  set2<-filter(alt_met,stock==use_stocks[y])[,-1]
  colnames(set2)[4]<-"species"
  
  casted<-dcast(set1,Year~process,value.var="values")
  colnames(casted)[4]<-"Other_mortality"
  mod_dat<-merge(casted,set2,by="Year")
  mod_dat$lag_temp<-c(NA,mod_dat$Temperature[-length(mod_dat$Temperature)])
  colnames(mod_dat)[colnames(mod_dat)=="Spawner abundance"]<-"Abundance"
  
  #==recruits
  mod_r<-gam(data=mod_dat,Recruitment~s(Abundance,k=4)+s(Temperature,k=4),family=nb(link = "log"))
  summary(mod_r)
  plotted<-plot(mod_r,pages=1)
  
  #ccf(mod_dat$Recruitment,mod_dat$Temperature,na.action=na.pass)
  preds<-predict.gam(mod_r,type='response',se.fit=TRUE)
  
  plo_gam_r<-data.frame(obs=mod_dat$Recruitment[!is.na(mod_dat$Recruitment)],
                        preds=preds$fit,
                        y_up=preds$fit+preds$se,
                        y_dn=preds$fit-preds$se,
                        year=mod_dat$Year[!is.na(mod_dat$Recruitment)],
                        stock=rep(use_stocks[y],length(preds$fit)))
  out_plot_r<-rbind(out_plot_r,plo_gam_r)
  
  for(x in 1:(length(plotted)))
  {
    temp<-data.frame(x=plotted[[x]]$x,
                     y=plotted[[x]]$fit,
                     y_up=plotted[[x]]$fit+plotted[[x]]$se,
                     y_dn=plotted[[x]]$fit-plotted[[x]]$se,
                     covar=plotted[[x]]$xlab,
                     stock=use_stocks[y])
    rec_term<-rbind(rec_term,temp)
  }
  
  
  #==immature mortality
  base_mod_imm_m<-gam(data=mod_dat,M_imm~1,family=tw())
  mod_imm_m<-gam(data=mod_dat,M_imm~s(N_imm,k=4)+s(Temperature,k=4)+s(Size,k=3),family=tw())
  keep_conc[[conc_cnt]]<- concurvity(mod_imm_m,full=FALSE)$observed[-1,-1]
  conc_cnt<-conc_cnt+1
  simout<-simulateResiduals(mod_imm_m,n=250)
  
  png(paste("plots/dharma_m",use_stocks[y],"_imm.png",sep=''),height=6,width=8,res=350,units='in') 
  plot(simout)
  dev.off()
  
  mod_imm_m_a<-gam(data=mod_dat,M_imm~s(N_imm,k=4),family=tw())
  mod_imm_m_t<-gam(data=mod_dat,M_imm~s(Temperature,k=4),family=tw())
  mod_imm_m_s<-gam(data=mod_dat,M_imm~s(Size,k=3),family=tw())
  
  summary(mod_imm_m)
  plot(mod_imm_m,pages=1)
  dev_expl_m<-c(dev_expl_m,round(summary(mod_imm_m)$dev,2))
  plotted<-plot(mod_imm_m,pages=1)
  
  keep_AIC_m<-rbind(keep_AIC_m,c(AIC(base_mod_imm_m),AIC(mod_imm_m)))
  keep_AIC_m_big<-rbind(keep_AIC_m_big,c(AIC(base_mod_imm_m),AIC(mod_imm_m),AIC(mod_imm_m_a),AIC(mod_imm_m_t),AIC(mod_imm_m_s)))
  
  preds<-predict.gam(mod_imm_m,type='response',se.fit=TRUE)
  
  plo_gam_m<-data.frame(obs=mod_dat$M_imm[!is.na(mod_dat$M_imm)],
                        preds=preds$fit,
                        y_up=preds$fit+preds$se,
                        y_dn=preds$fit-preds$se,
                        year=mod_dat$Year[!is.na(mod_dat$M_imm)],
                        stock=paste(use_stocks[y],"_imm",sep=''))
  
  out_plot_m<-rbind(out_plot_m,plo_gam_m)
  
  for(x in 1:(length(plotted)))
  {
    temp<-data.frame(x=plotted[[x]]$x,
                     y=plotted[[x]]$fit,
                     y_up=plotted[[x]]$fit+plotted[[x]]$se,
                     y_dn=plotted[[x]]$fit-plotted[[x]]$se,
                     covar=plotted[[x]]$xlab,
                     stock=paste(use_stocks[y],"_imm",sep=''))
    mort_term<-rbind(mort_term,temp)
  }
  
  #==mature mortality
  base_mod_imm_m<-gam(data=mod_dat,Other_mortality~1,family=tw())
  mod_imm_m<-gam(data=mod_dat,Other_mortality~s(Abundance ,k=4)+s(Temperature,k=4)+s(Size,k=3),family=tw())
  mod_imm_m_a<-gam(data=mod_dat,Other_mortality~s(Abundance,k=4),family=tw())
  mod_imm_m_t<-gam(data=mod_dat,Other_mortality~s(Temperature,k=4),family=tw())
  mod_imm_m_s<-gam(data=mod_dat,Other_mortality~s(Size,k=3),family=tw())
  keep_conc[[conc_cnt]]<- concurvity(mod_imm_m,full=FALSE)$observed[-1,-1]
  conc_cnt<-conc_cnt+1
  simout<-simulateResiduals(mod_imm_m,n=250)
  
  png(paste("plots/dharma_m",use_stocks[y],"_mat.png",sep=''),height=6,width=8,res=350,units='in') 
  plot(simout)
  dev.off()
  
  
  summary(mod_imm_m)
  plot(mod_imm_m,pages=1)
  dev_expl_m<-c(dev_expl_m,round(summary(mod_imm_m)$dev,2))
  plotted<-plot(mod_imm_m,pages=1)
  
  keep_AIC_m<-rbind(keep_AIC_m,c(AIC(base_mod_imm_m),AIC(mod_imm_m)))
  preds<-predict.gam(mod_imm_m,type='response',se.fit=TRUE)
  keep_AIC_m_big<-rbind(keep_AIC_m_big,c(AIC(base_mod_imm_m),AIC(mod_imm_m),AIC(mod_imm_m_a),AIC(mod_imm_m_t),AIC(mod_imm_m_s)))
  
  plo_gam_m<-data.frame(obs=mod_dat$Other_mortality[!is.na(mod_dat$Other_mortality)],
                        preds=preds$fit,
                        y_up=preds$fit+preds$se,
                        y_dn=preds$fit-preds$se,
                        year=mod_dat$Year[!is.na(mod_dat$Other_mortality)],
                        stock=paste(use_stocks[y],"_mat",sep=''))
  
  out_plot_m<-rbind(out_plot_m,plo_gam_m)
  
  for(x in 1:(length(plotted)))
  {
    temp<-data.frame(x=plotted[[x]]$x,
                     y=plotted[[x]]$fit,
                     y_up=plotted[[x]]$fit+plotted[[x]]$se,
                     y_dn=plotted[[x]]$fit-plotted[[x]]$se,
                     covar=plotted[[x]]$xlab,
                     stock=paste(use_stocks[y],"_mat",sep=''))
    mort_term<-rbind(mort_term,temp)
  }
  
}

# ============================================================
# Section 2: GAM model summary figures used by the manuscript
# Outputs: plots/fig_6.png and plots_select/all_mort_all3.png
# Original source: 04_models_beta_select.R
# ============================================================
library(reshape2)
library(ggplot2)
library(mgcv)  
library(dplyr)
library(ggridges)
library(png)
library(PBSmodelling)
library(patchwork)
library(png)
library(grid)
library(DHARMa)
library(ggplot2)
#library(itsadug)
#library(mgcViz)
#library(directlabels)
outs<-read.csv("data/all_output.csv")
outs<-filter(outs,Year<2026)
alt_met<-read.csv("data/alt_metrics_calc.csv")
unc_mort<-read.csv("data/uncertainty_mort.csv")
other_met<-read.csv("data/ice_extent.csv")

#======================================
# full models:
# recruitment <- density + temperature + competitor + size
#==show two things
#=====1. some of the processes covary
#=====2. some of the processes can be explained by environmental variables

#==two figures: one for recruitment, one for mortality
#==first column are the estimates and the model fits
#==remaining columns represent the variables
#==only 'significant' variables are colored in
#==CHECK ON APPROPRIATE LAGS FOR RECRUITMENT
#==THINK ABOUT LARGE SCALE DRIVERS TOO

#==king crabs
use_stocks<-c("BBRKC","PIRKC","SMBKC","PIBKC")
rec_term<-NULL
mort_term<-NULL
out_plot_r<-NULL
out_plot_m<-NULL
dev_expl_m<-NULL
keep_AIC_r<-NULL
keep_AIC_m<-NULL
keep_AIC_m_big<-NULL
keep_conc<-list(list())
conc_cnt<-1

for(y in 1:length(use_stocks))
{
set1<-filter(outs,species==use_stocks[y])[,-c(1,6)]
set2<-filter(alt_met,stock==use_stocks[y])[,-1]
set3<-filter(unc_mort,stock==use_stocks[y])[,c(2,6,7)]
colnames(set2)[4]<-"species"

casted<-dcast(set1,Year~process,value.var="values")
colnames(casted)[4]<-"Other_mortality"
mod_dat<-merge(casted,set2,by="Year")
mod_dat<-mod_dat[,-3]
mod_dat<-merge(mod_dat,other_met,by="Year",all=TRUE)
mod_dat<-merge(mod_dat,set3,by="Year",all=TRUE)
mod_dat$p_mort<-1-exp(-mod_dat$Other_mortality)
if(use_stocks[y]=="SMBKC")
  mod_dat<-mod_dat[-1,]
mod_dat$lag_temp<-c(NA,mod_dat$Temperature[-length(mod_dat$Temperature)])

#==p(mort)
mod_dat_base<-mod_dat[,c(2,6,7,11,14,15)]
mod_dat_base<-mod_dat_base[complete.cases(mod_dat_base),]

base_mod_m<-gam(data=mod_dat_base,p_mort~1,family = betar(link = "logit"),weights=1/sd)
mod_m<-gam(data=mod_dat,p_mort~s(Abundance,k=4,bs='ts')+s(Temperature,k=4,bs='ts')+s(Size,k=3,bs='ts')+s(Ice,k=3,bs='ts'),
           family = betar(link = "logit"),weights=1/sd,
           select=TRUE,method='REML')

#mod_m<-gam(data=mod_dat,p_mort~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3),family = betar(link = "logit"),weights=1/sd)

mod_m_a<-gam(data=mod_dat_base,p_mort~s(Abundance,k=4),family = betar(link = "logit"),weights=1/sd)
mod_m_t<-gam(data=mod_dat_base,p_mort~s(Temperature,k=4),family = betar(link = "logit"),weights=1/sd)
mod_m_s<-gam(data=mod_dat_base,p_mort~s(Size,k=3),family = betar(link = "logit"),weights=1/sd)
mod_m_i<-gam(data=mod_dat_base,p_mort~s(Ice,k=3),family = betar(link = "logit"),weights=1/sd)

keep_conc[[conc_cnt]]<- concurvity(mod_m,full=FALSE)$observed[-1,-1]
   conc_cnt<-conc_cnt+1

dev_expl_m<-c(dev_expl_m,round(summary(mod_m)$dev,2))
keep_AIC_m<-rbind(keep_AIC_m,c(AIC(base_mod_m),AIC(mod_m)))
keep_AIC_m_big<-rbind(keep_AIC_m_big,c(AIC(base_mod_m),AIC(mod_m),AIC(mod_m_a),AIC(mod_m_t),AIC(mod_m_s),AIC(mod_m_i)))

plotted<-plot(mod_m,pages=1)

preds<-predict.gam(mod_m,type='response',se.fit=TRUE)

plo_gam_m<-data.frame(obs=mod_dat$p_mort[as.numeric(names(preds$fit))],
                      preds=preds$fit,
                      y_up=preds$fit+2*preds$se,
                      y_dn=preds$fit-2*preds$se,
                      year=mod_dat$Year[as.numeric(names(preds$fit))],
                      stock=use_stocks[y])

out_plot_m<-rbind(out_plot_m,plo_gam_m)


for(x in 1:(length(plotted)))
{
  temp<-data.frame(x=plotted[[x]]$x,
                   y=plotted[[x]]$fit,
                   y_up=plotted[[x]]$fit+2*plotted[[x]]$se,
                   y_dn=plotted[[x]]$fit-2*plotted[[x]]$se,
                   covar=plotted[[x]]$xlab,
                   stock=use_stocks[y])
  mort_term<-rbind(mort_term,temp)
}

}
do_prib<-0
#==do a Pribs model
if(do_prib==1)
{
set1_rkc<-filter(outs,species%in%c("PIRKC"))[,-c(1,6)]
set2_rkc<-filter(alt_met,stock%in%c("PIRKC"))[,-1]
colnames(set2)[4]<-"species"

casted<-dcast(set1_rkc,Year~process,value.var="values")
colnames(casted)[4]<-"Other_mortality"
mod_dat_rkc<-merge(casted,set2_rkc,by="Year")

set1_bkc<-filter(outs,species%in%c("PIBKC"))[,-c(1,6)]
set2_bkc<-filter(alt_met,stock%in%c("PIBKC"))[,-1]
colnames(set2)[4]<-"species"

casted<-dcast(set1_bkc,Year~process,value.var="values")
colnames(casted)[4]<-"Other_mortality"
mod_dat_bkc<-merge(casted,set2_bkc,by="Year")

#==predict pirkc with pibkc abundance 
mod_dat<-mod_dat_rkc
mod_dat$bkc_n<-mod_dat_bkc$Abundance[-1]
mod_dat$bkc_ssb<-mod_dat_bkc$`Spawner abundance`[-1]
mod_m<-gam(data=mod_dat,Other_mortality~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3)+s(bkc_ssb,k=3)+s(bkc_n,k=3),family=tw())
summary(mod_m)
plot(mod_m,pages=1)

#==predict pibkc M with pirkc abundance 
mod_dat<-mod_dat_bkc[-1,]
mod_dat$rkc_n<-mod_dat_rkc$Abundance
mod_dat$rkc_ssb<-mod_dat_rkc$`Spawner abundance`
mod_m<-gam(data=mod_dat,Other_mortality~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3),family=tw())
summary(mod_m)
plot(mod_m,pages=1)

#==predict pibkc R with pirkc abundance 
mod_dat<-mod_dat_bkc[-1,]
mod_dat$rkc_n<-mod_dat_rkc$Abundance
mod_dat$rkc_ssb<-mod_dat_rkc$`Spawner abundance`
mod_m<-gam(data=mod_dat,Recruitment ~s(Abundance,k=3)+s(Temperature,k=3)+s(Size,k=3)+s(rkc_n,k=3),family=tw())
summary(mod_m)
plot(mod_m,pages=1)
}
#============================
# chionoecetes species
#============================
use_stocks<-c("Snow","Tanner")
for(y in 1:length(use_stocks))
{
  set1<-filter(outs,species==use_stocks[y])[,-c(1,6)]
  set2<-filter(alt_met,stock==use_stocks[y])[,-1]
  tmp_unc<-unc_mort[grep(use_stocks[y],unc_mort$stock),]
  imm_sd<-tmp_unc[grep('imm',tmp_unc$stock),c(6,7)]
  colnames(imm_sd)[1]<-'imm_sd'
  mat_sd<-tmp_unc[grep('mat',tmp_unc$stock),c(6,7)]
  colnames(mat_sd)[1]<-'mat_sd'

  casted<-dcast(set1,Year~process,value.var="values")
  colnames(casted)[4]<-"Other_mortality"
  mod_dat<-merge(casted,set2,by="Year")
  colnames(mod_dat)[colnames(mod_dat)=="Spawner abundance"]<-"Abundance"
  mod_dat<-mod_dat[,-2]
  mod_dat<-merge(mod_dat,other_met,by="Year",all=TRUE)
  mod_dat<-merge(mod_dat,imm_sd,by="Year",all=TRUE)
  mod_dat<-merge(mod_dat,mat_sd,by="Year",all=TRUE)
  mod_dat$p_mort_mat<-1-exp(-mod_dat$Other_mortality)
  mod_dat$p_mort_imm<-1-exp(-mod_dat$M_imm)    
  #write.csv(mod_dat,"data/AK_gam_dat_mat.csv")
  #==p(mort) mat
  mod_dat_base<-mod_dat[,c(6,7,8,12,14,15,16,17)]
  mod_dat_base<-mod_dat_base[complete.cases(mod_dat_base),]
  base_mod_m<-gam(data=mod_dat_base,p_mort_mat~1,family = betar(link = "logit"),weights=1/mat_sd)
  mod_mat_m<-gam(data=mod_dat,p_mort_mat~s(Abundance,k=4,bs='ts')+s(Temperature,k=4,bs='ts')+s(Size,k=3,bs='ts')+s(Ice,k=3,bs='ts'),
                 family = betar(link = "logit"),weights=1/mat_sd,
                 select=TRUE,method='REML')
  
  # plot(mod_mat_m,pages=1)
  # summary(mod_mat_m)
   #mod_mat_m<-gam(data=mod_dat,p_mort_mat~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3),family = betar(link = "logit"),weights=1/mat_sd)

  mod_mat_a<-gam(data=mod_dat_base,p_mort_mat~s(Abundance,k=4),family = betar(link = "logit"),weights=1/mat_sd)
  mod_mat_t<-gam(data=mod_dat_base,p_mort_mat~s(Temperature,k=4),family = betar(link = "logit"),weights=1/mat_sd)
  mod_mat_s<-gam(data=mod_dat_base,p_mort_mat~s(Size,k=3),family = betar(link = "logit"),weights=1/mat_sd)
  mod_mat_i<-gam(data=mod_dat_base,p_mort_mat~s(Ice,k=3),family = betar(link = "logit"),weights=1/mat_sd)

  # 
  keep_conc[[conc_cnt]]<- concurvity(mod_mat_m,full=FALSE)$observed[-1,-1]
  conc_cnt<-conc_cnt+1
  dev_expl_m<-c(dev_expl_m,round(summary(mod_mat_m)$dev,2))
  plotted<-plot(mod_mat_m,pages=1)
  
  keep_AIC_m<-rbind(keep_AIC_m,c(AIC(base_mod_m),AIC(mod_mat_m)))
  keep_AIC_m_big<-rbind(keep_AIC_m_big,c(AIC(base_mod_m),AIC(mod_mat_m),AIC(mod_mat_a),AIC(mod_mat_t),AIC(mod_mat_s),AIC(mod_mat_i)))
  
  preds<-predict.gam(mod_mat_m,type='response',se.fit=TRUE)
  
  plo_gam_m<-data.frame(obs=mod_dat$p_mort_mat[as.numeric(names(preds$fit))],
                        preds=preds$fit,
                        y_up=preds$fit+2*preds$se,
                        y_dn=preds$fit-2*preds$se,
                        year=mod_dat$Year[as.numeric(names(preds$fit))],
                        stock=paste(use_stocks[y],"_mat",sep=''))
  
  out_plot_m<-rbind(out_plot_m,plo_gam_m)
  
  for(x in 1:(length(plotted)))
  {
    temp<-data.frame(x=plotted[[x]]$x,
                     y=plotted[[x]]$fit,
                     y_up=plotted[[x]]$fit+2*plotted[[x]]$se,
                     y_dn=plotted[[x]]$fit-2*plotted[[x]]$se,
                     covar=plotted[[x]]$xlab,
                     stock=paste(use_stocks[y],"_mat",sep=''))
    mort_term<-rbind(mort_term,temp)
  }
  
  #==p(mort) imm
  base_mod_m<-gam(data=mod_dat_base,p_mort_imm~1,family = betar(link = "logit"),weights=1/imm_sd)
  mod_imm_m<-gam(data=mod_dat,p_mort_imm~s(Abundance,k=4,bs='ts')+s(Temperature,k=4,bs='ts')+s(Size,k=3,bs='ts')+s(Ice,k=3,bs='ts'),
                 family = betar(link = "logit"),weights=1/imm_sd,
                 select=TRUE,method='REML')
  # mod_imm_m<-gam(data=mod_dat,p_mort_imm~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3),family = betar(link = "logit"),weights=1/imm_sd)
  plot(mod_imm_m,pages=1)
  summary(mod_imm_m)
  mod_imm_a<-gam(data=mod_dat_base,p_mort_imm~s(Abundance,k=4),family = betar(link = "logit"),weights=1/imm_sd)
  mod_imm_t<-gam(data=mod_dat_base,p_mort_imm~s(Temperature,k=4),family = betar(link = "logit"),weights=1/imm_sd)
  mod_imm_s<-gam(data=mod_dat_base,p_mort_imm~s(Size,k=3),family = betar(link = "logit"),weights=1/imm_sd)
  mod_imm_i<-gam(data=mod_dat_base,p_mort_imm~s(Ice,k=3),family = betar(link = "logit"),weights=1/imm_sd)
  summary(mod_imm_s)
  
  keep_conc[[conc_cnt]]<- concurvity(mod_imm_m,full=FALSE)$observed[-1,-1]
  conc_cnt<-conc_cnt+1
  
  dev_expl_m<-c(dev_expl_m,round(summary(mod_imm_m)$dev,2))
  plotted<-plot(mod_imm_m,pages=1)
  
  keep_AIC_m<-rbind(keep_AIC_m,c(AIC(base_mod_m),AIC(mod_imm_m)))
  keep_AIC_m_big<-rbind(keep_AIC_m_big,c(AIC(base_mod_m),AIC(mod_imm_m),AIC(mod_imm_a),AIC(mod_imm_t),AIC(mod_imm_s),AIC(mod_imm_i)))
  
  preds<-predict.gam(mod_imm_m,type='response',se.fit=TRUE)
  
  plo_gam_m<-data.frame(obs=mod_dat$M_imm[as.numeric(names(preds$fit))],
                        preds=preds$fit,
                        y_up=preds$fit+2*preds$se,
                        y_dn=preds$fit-2*preds$se,
                        year=mod_dat$Year[as.numeric(names(preds$fit))],
                        stock=paste(use_stocks[y],"_imm",sep=''))
  
  out_plot_m<-rbind(out_plot_m,plo_gam_m)
  
  for(x in 1:(length(plotted)))
  {
    temp<-data.frame(x=plotted[[x]]$x,
                     y=plotted[[x]]$fit,
                     y_up=plotted[[x]]$fit+2*plotted[[x]]$se,
                     y_dn=plotted[[x]]$fit-2*plotted[[x]]$se,
                     covar=plotted[[x]]$xlab,
                     stock=paste(use_stocks[y],"_imm",sep=''))
    mort_term<-rbind(mort_term,temp)
  }
  
 
  
}

colnames(keep_AIC_m_big)<-c("No covars","All covars","Density","Temp","Size","Ice")

rownames(keep_AIC_m_big)<-c("BBRKC","PIRKC","SMBKC","PIBKC","Snow (mat)","Snow (imm)","Tanner (mat)","Tanner (imm)")

names(dev_expl_m)<-c("BBRKC","PIRKC","SMBKC","PIBKC","Snow (mat)","Snow (imm)","Tanner (mat)","Tanner (imm)")

plot_conc<-NULL
for(x in 1:length(keep_conc))
{
 tmp<-melt(keep_conc[[x]]) 
 #levels(tmp$Var1)<-c("s(Abundance)","s(Size)","s(Temperature)","s(Ice)")
 #levels(tmp$Var2)<-c("s(Abundance)","s(Size)","s(Temperature)","s(Ice)")
 tmp$stock<-names(dev_expl_m)[x]
 plot_conc<-rbind(plot_conc,tmp)
}

plot_conc[plot_conc==1]<-NA
p_conc<-ggplot(plot_conc)+
  geom_tile(aes(x=Var1,y=Var2,fill=value))+
  geom_text(aes(x=Var1,y=Var2,label=round(value,2)))+
  scale_fill_gradient(low='white',high='red',na.value='whitesmoke')+
  facet_wrap(~stock)+theme_bw()+
  theme(legend.position=c(.85,.15),
        axis.text.x = element_text(angle = 45, hjust = 1) )+
  #ggtitle("Estimated concurvity among variables by model")+
  ylab("")+xlab("")+guides(fill=guide_legend(title="Concurvity"))


alt_AIC<-keep_AIC_m_big
for(x in 1:nrow(alt_AIC))
  for(y in 2:ncol(alt_AIC))
    alt_AIC[x,y]<--1*(alt_AIC[x,1]-alt_AIC[x,y])

in_dat_t1<-as.data.frame(round(alt_AIC,2))

library(tidyr)
library(scales)
# Specify the columns to be colored
columns_to_color <- colnames(in_dat_t1[2:6])

# Reshape the data for ggplot2
data_long <- in_dat_t1 %>%
  mutate(RowID = rownames(in_dat_t1)) %>%
  pivot_longer(cols = all_of(columns_to_color), names_to = "Variable", values_to = "Value")

# Create the heatmap-like graphic
aic_table<-ggplot(data_long, aes(x = Variable, y = as.factor(RowID), fill = Value)) +
  geom_tile(color = "black") + # Add borders to cells
  scale_fill_gradient2(
    low = "forestgreen",
    mid = "white",
    high = "tomato3",
    midpoint = 0,
    space = "Lab",
    limits=c(-10,10),
    oob=squish
  ) +
  geom_text(aes(label = Value), color = "black") + # Add cell values as text
  labs( x = "",y='') +
  theme_minimal()+
  theme(legend.position="none")

use_dat<-melt(dev_expl_m)
colnames(use_dat)[1]<-"Deviance_explained"
use_dat$stock<-rownames(use_dat)

my_values <- dev_expl_m

# Convert the vector to a data frame, necessary for ggplot2
df <- data.frame(
  Value = my_values,
  Row = seq_along(my_values),  # Assign a row number for each value
  Row_Name = names(my_values)         # Add the row names to the data frame
)

# Create the ggplot
dev_ex_plot<-ggplot(df, aes(x = 1, y = as.factor(Row_Name), fill = Value)) + # Map x to a single column (1), y to Row, and fill to Value
  geom_tile(color = "black") +  # Create tiles, add a black border for clarity
  scale_fill_gradient(low = "white", high = "cornflowerblue") +  # Apply the gradient
  geom_text(aes(label = Value), color = "black") + # Add the value as text inside each tile
 # geom_text(aes(x = 0.5, label = Row_Name), hjust = 1, color = "black") + # Add row names to the left of the tiles
  labs( x = NULL, y = NULL) +  # Add title and remove axis labels
  theme_minimal() + 
  xlim(0,2)+
  theme(axis.text.x = element_blank(), # Remove x-axis text
        axis.ticks.x = element_blank(), # Remove x-axis ticks
        axis.text.y = element_blank(), # Remove x-axis text
        axis.ticks.y = element_blank(), # Remove x-axis ticks
        panel.grid = element_blank(),  # Remove grid lines
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, lineheight = 1.2) )  # Place the legend on the right



png("plots/fig_6.png",height=7,width=12,res=350,units='in') 

# 1. Ensure your original titles are attached to each individual plot
p_conc <- p_conc + 
  ggtitle("A. Estimated concurvity\namong variables by model")

aic_table <- aic_table + 
  labs(title = "B. Change in AIC\nfrom null model") +
  # Add the dummy facet to match structural height
  facet_wrap(~ "dummy") + 
  theme(
    strip.background = element_blank(), # Hides the gray box
    strip.text = element_blank()        # Hides the "dummy" text
  )

dev_ex_plot <- dev_ex_plot + 
  labs(title = "C. Deviance\nexplained") +
  # Add the dummy facet to match structural height
  facet_wrap(~ "dummy") + 
  theme(
    strip.background = element_blank(), # Hides the gray box
    strip.text = element_blank()        # Hides the "dummy" text
  )

# 2. Combine using your exact layout call
p_conc + aic_table + dev_ex_plot + plot_layout(widths = c(4, 2, 1))

dev.off()

in_col2<-c("#ff5050","#0034c377","#ff505077","#0034c3","#3da550","#3da55099","#ff8f38","#ff8f3899")
mort_term$covar[mort_term$covar=="N_imm"]<-"Abundance"

ughh<-out_plot_m[,c(1,5,6)]
colnames(ughh)<-c("est_m","year","stock'")

mod_mort<-unc_mort
mod_mort$up_m[mod_mort$up_m>max(mod_mort$est_m)]<-1.5*max(mod_mort$est_m)

mod_mort_cap <- mod_mort %>%
  group_by(stock) %>%
  mutate(up_m = pmin(up_m, 1.5 * max(est_m, na.rm = TRUE))) %>% # The `pmin` function is used here to take the minimum of the two values element-wise, effectively capping the up_m values.
  ungroup()

input<-as.data.frame(filter(unc_mort,Year<2026&Year!=2020))
input <- input %>%
  mutate(
    Year = as.numeric(as.character(Year)),
    est_m = as.numeric(as.character(est_m))
  )

input[input$stock=="Snow_mat",2]<-"Snow (mature)"
input[input$stock=="Snow_imm",2]<-"Snow (immature)"
input[input$stock=="Tanner_mat",2]<-"Tanner (mature)"
input[input$stock=="Tanner_imm",2]<-"Tanner (immature)"

out_plot_m[out_plot_m$stock=="Snow_mat",6]<-"Snow (mature)"
out_plot_m[out_plot_m$stock=="Snow_imm",6]<-"Snow (immature)"
out_plot_m[out_plot_m$stock=="Tanner_mat",6]<-"Tanner (mature)"
out_plot_m[out_plot_m$stock=="Tanner_imm",6]<-"Tanner (immature)"

mort_plot_alt_trans<-ggplot()+
  geom_point(data=input,aes(x=Year,y=1-exp(-est_m)),col='darkgrey')+
  geom_errorbar(data=input,aes(x=Year,ymin=1-exp(-dn_m),ymax=1-exp(-up_m)),col='grey')+
  geom_ribbon(data=out_plot_m,aes(x=year,ymin=y_dn,ymax=y_up,fill=stock),alpha=0.3,lwd=2)+
  geom_line(data=out_plot_m,aes(x=year,y=preds,col=stock),lwd=1.15)+  theme_bw()+ylab("p(mortality)")+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)+
  facet_wrap(~stock,ncol=1,scales='free_y')+
  theme(legend.position='none')+ylim(0,1)


mort_plot_trm<-ggplot(mort_term)+
  geom_line(aes(x=x,y=y,col=stock),lwd=2)+
  geom_ribbon(aes(x=x,ymin=y_dn,ymax=y_up,fill=stock),alpha=0.3,lwd=2)+
  theme_bw()+
  theme(legend.position='none',
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        legend.key.size=unit(.7, 'lines'))+
  geom_hline(yintercept=0,lty=2)+xlab("Observed value")+ylab("Smooth")+
  facet_grid(rows=vars(stock),cols=vars(covar),scales='free_x')+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)

mort_plot_1<-ggplot(mort_term)+
  geom_line(aes(x=x,y=y,col=stock),lwd=1.25)+
  geom_ribbon(aes(x=x,ymin=y_dn,ymax=y_up,fill=stock),alpha=0.2,lwd=2)+
  theme_bw()+
  theme(legend.position='none',
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        legend.key.size=unit(.7, 'lines'))+
  geom_hline(yintercept=0,lty=2)+xlab("Observed value")+ylab("Smooth")+
  facet_wrap(~covar,scales='free',ncol=1)+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)
mort_plot_2<-ggplot(mort_term)+
  geom_line(aes(x=x,y=y,col=stock),lwd=1.25)+
  geom_ribbon(aes(x=x,ymin=y_dn,ymax=y_up,fill=stock),alpha=0.2,lwd=2)+
  theme_bw()+
  theme(legend.position='none',
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        legend.key.size=unit(.7, 'lines'))+
  geom_hline(yintercept=0,lty=2)+xlab("Observed value")+ylab("Smooth")+
  facet_grid(stock~covar,scales='free')+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)

# mortality <- density + temperature + competitor + size
library(patchwork)
design<-"112
         112
         112"
png("plots_select/all_mort_all3.png",height=9,width=6,res=350,units='in') 
mort_plot_alt_trans +mort_plot_1 + plot_layout(nrow=2, design=design)
dev.off()

# ============================================================
# Section 3: Leave-one-out cross-validation figures
# Outputs: plots/model_crossval.png and plots/model_crossval_pred.png
# Original source: 04_models_cross_validate.R
# ============================================================
library(reshape2)
library(ggplot2)
library(mgcv)  
library(dplyr)
library(ggridges)
library(png)
library(PBSmodelling)
library(patchwork)
library(mgcv)
library(png)
library(grid)
library(DHARMa)
library(itsadug)
library(mgcViz)
library(directlabels)
outs<-read.csv("data/all_output.csv")
outs<-filter(outs,Year<2023)
alt_met<-read.csv("data/alt_metrics_calc.csv")
unc_mort<-read.csv("data/uncertainty_mort.csv")
other_met<-read.csv("data/coldpool_ice.csv")

big_keep_all<-NULL
out_plot_m_cv_all<-NULL

#==king crabs
use_stocks<-c("BBRKC","PIRKC","SMBKC","PIBKC")
for(y in 1:length(use_stocks))
{
  set1<-filter(outs,species==use_stocks[y])[,-c(1,6)]
  set2<-filter(alt_met,stock==use_stocks[y])[,-1]
  set3<-filter(unc_mort,stock==use_stocks[y])[,c(2,6,7)]
  colnames(set2)[4]<-"species"
  
  casted<-dcast(set1,Year~process,value.var="values")
  colnames(casted)[4]<-"Other_mortality"
  mod_dat<-merge(casted,set2,by="Year")
  mod_dat<-mod_dat[,-3]
  mod_dat<-merge(mod_dat,other_met,by="Year",all=TRUE)
  mod_dat<-merge(mod_dat,set3,by="Year",all=TRUE)
  mod_dat$p_mort<-1-exp(-mod_dat$Other_mortality)
  if(use_stocks[y]=="SMBKC")
    mod_dat<-mod_dat[-1,]
  mod_dat$lag_temp<-c(NA,mod_dat$Temperature[-length(mod_dat$Temperature)])
  
  #==p(mort)
  mod_dat_base<-mod_dat[,c(1,2,6,7,11,14,15)]
  mod_dat_base<-mod_dat_base[complete.cases(mod_dat_base),]
  big_keep<-NULL
  out_plot_m_cv<-NULL
  for(x in 1:nrow(mod_dat_base))
  {
    #=do LOOCV
   cv_dat<-mod_dat_base[-x,]
   mod_m<-gam(data=cv_dat,p_mort~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3)+s(ice,k=3),family = betar(link = "logit"),weights=1/sd)
   keeper<-as.data.frame(summary(mod_m)$s.table)
   keeper$stock<-use_stocks[y]
   keeper$smooth<-rownames(keeper)
   big_keep<-rbind(big_keep,keeper)
   
   #=record predictions
   preds<-predict.gam(mod_m,type='response',se.fit=TRUE)
   plo_gam_m_cv<-data.frame(obs=cv_dat$p_mort,
                         preds=preds$fit,
                         y_up=preds$fit+2*preds$se,
                         y_dn=preds$fit-2*preds$se,
                         year=cv_dat$Year,
                         stock=use_stocks[y],
                         group=x)
   out_plot_m_cv<-rbind(out_plot_m_cv,plo_gam_m_cv)
  }
  big_keep_all<-rbind(big_keep_all,big_keep)
  out_plot_m_cv_all<-rbind(out_plot_m_cv_all,out_plot_m_cv)
}
   
#==chionoecetes species
use_stocks<-c("Snow","Tanner")
for(y in 1:length(use_stocks))
{
  set1<-filter(outs,species==use_stocks[y])[,-c(1,6)]
  set2<-filter(alt_met,stock==use_stocks[y])[,-1]
  tmp_unc<-unc_mort[grep(use_stocks[y],unc_mort$stock),]
  imm_sd<-tmp_unc[grep('imm',tmp_unc$stock),c(6,7)]
  colnames(imm_sd)[1]<-'imm_sd'
  mat_sd<-tmp_unc[grep('mat',tmp_unc$stock),c(6,7)]
  colnames(mat_sd)[1]<-'mat_sd'
  
  casted<-dcast(set1,Year~process,value.var="values")
  colnames(casted)[4]<-"Other_mortality"
  mod_dat<-merge(casted,set2,by="Year")
  colnames(mod_dat)[colnames(mod_dat)=="Spawner abundance"]<-"Abundance"
  mod_dat<-mod_dat[,-2]
  mod_dat<-merge(mod_dat,other_met,by="Year",all=TRUE)
  mod_dat<-merge(mod_dat,imm_sd,by="Year",all=TRUE)
  mod_dat<-merge(mod_dat,mat_sd,by="Year",all=TRUE)
  mod_dat$p_mort_mat<-1-exp(-mod_dat$Other_mortality)
  mod_dat$p_mort_imm<-1-exp(-mod_dat$M_imm)    
  
  #==p(mort) mat
  mod_dat_base<-mod_dat[,c(1,6,7,8,12,14,15,16,17)]
  mod_dat_base<-mod_dat_base[complete.cases(mod_dat_base),]
  big_keep<-NULL
  out_plot_m_cv<-NULL
  for(x in 1:nrow(mod_dat_base))
  {
    #=do LOOCV
    cv_dat<-mod_dat_base[-x,]
    mod_m<-gam(data=cv_dat,p_mort_mat~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3)+s(ice,k=3),family = betar(link = "logit"),weights=1/mat_sd)
    keeper<-as.data.frame(summary(mod_m)$s.table)
    keeper$stock<-paste(use_stocks[y],"_mat",sep='')
    keeper$smooth<-rownames(keeper)
    big_keep<-rbind(big_keep,keeper)
    
    #=record predictions
    preds<-predict.gam(mod_m,type='response',se.fit=TRUE)
    plo_gam_m_cv<-data.frame(obs=cv_dat$p_mort_mat,
                             preds=preds$fit,
                             y_up=preds$fit+2*preds$se,
                             y_dn=preds$fit-2*preds$se,
                             year=cv_dat$Year,
                             stock=paste(use_stocks[y],"_mat",sep=''),
                             group=x)
    out_plot_m_cv<-rbind(out_plot_m_cv,plo_gam_m_cv)
  }
  big_keep_all<-rbind(big_keep_all,big_keep)
  out_plot_m_cv_all<-rbind(out_plot_m_cv_all,out_plot_m_cv)
  
  big_keep<-NULL
  out_plot_m_cv<-NULL
  #==p(mort) imm
  for(x in 1:nrow(mod_dat_base))
  {
    #=do LOOCV
    cv_dat<-mod_dat_base[-x,]
    mod_m<-gam(data=cv_dat,p_mort_imm~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3)+s(ice,k=3),family = betar(link = "logit"),weights=1/imm_sd)
    keeper<-as.data.frame(summary(mod_m)$s.table)
    keeper$stock<-paste(use_stocks[y],"_imm",sep='')
    keeper$smooth<-rownames(keeper)
    big_keep<-rbind(big_keep,keeper)
    
    #=record predictions
    preds<-predict.gam(mod_m,type='response',se.fit=TRUE)
    plo_gam_m_cv<-data.frame(obs=cv_dat$p_mort_imm,
                             preds=preds$fit,
                             y_up=preds$fit+2*preds$se,
                             y_dn=preds$fit-2*preds$se,
                             year=cv_dat$Year,
                             stock=paste(use_stocks[y],"_imm",sep=''),
                             group=x)
    out_plot_m_cv<-rbind(out_plot_m_cv,plo_gam_m_cv)
  }
  big_keep_all<-rbind(big_keep_all,big_keep)
  out_plot_m_cv_all<-rbind(out_plot_m_cv_all,out_plot_m_cv)
  


}

hist_in<-as.data.frame(big_keep_all[,c(4,5,6)])
dumbo<-data.frame(pval=big_keep_all[,4],
                  stock=big_keep_all[,5],
                  smooth=big_keep_all[,6])
png("plots/model_crossval.png",height=6,width=6,res=350,units='in') 
ggplot(dumbo)+
  geom_histogram(aes(x=pval,fill=smooth,group=smooth))+
  facet_wrap(~stock)+theme_bw()+
  #xlim(-0.01,.9)+
  geom_vline(xintercept=0.05,lty=2)+
  theme(legend.position=c(.85,.15))
dev.off()

in_col2<-c("#ff5050","#0034c377","#ff505077","#0034c3","#3da550","#3da55099","#ff8f38","#ff8f3899")


mort_plot_cv<-ggplot()+
  geom_point(data=filter(unc_mort,Year<2023&Year!=2020),aes(x=Year,y=1-exp(-est_m)),col='darkgrey')+
  geom_errorbar(data=filter(unc_mort,Year<2023&Year!=2020),aes(x=Year,ymin=1-exp(-dn_m),ymax=1-exp(-up_m)),col='grey')+
  geom_line(data=out_plot_m_cv_all,aes(x=year,y=preds,col=stock,group=group),lwd=1,alpha=0.7)+  theme_bw()+ylab("p(mortality)")+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)+
  facet_wrap(~stock,ncol=1,scales='free_y')+
  theme(legend.position='none')+ylim(0,1)

png("plots/model_crossval_pred.png",height=10,width=6,res=350,units='in') 
print(mort_plot_cv)
dev.off()
