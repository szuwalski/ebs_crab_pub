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

colnames(keep_AIC_m_big)<-c("No covars","All covars","Density","Temp","Size")

rownames(keep_AIC_m_big)<-c("BBRKC","PIRKC","SMBKC","PIBKC","Snow (imm)","Snow (mat)","Tanner (imm)","Tanner (mat)")

names(dev_expl_m)<-c("BBRKC","PIRKC","SMBKC","PIBKC","Snow (imm)","Snow (mat)","Tanner (imm)","Tanner (mat)")

plot_conc<-NULL
for(x in 1:length(keep_conc))
{
 tmp<-melt(keep_conc[[x]]) 
 levels(tmp$Var1)<-c("s(Abundance)","s(Size)","s(Temperature)")
 levels(tmp$Var2)<-c("s(Abundance)","s(Size)","s(Temperature)")
 tmp$stock<-names(dev_expl_m)[x]
 plot_conc<-rbind(plot_conc,tmp)
}

plot_conc[plot_conc==1]<-NA
png("plots/model_concurvity.png",height=7,width=7,res=350,units='in') 
ggplot(plot_conc)+
  geom_tile(aes(x=Var1,y=Var2,fill=value))+
  geom_text(aes(x=Var1,y=Var2,label=round(value,2)))+
  scale_fill_gradient(low='white',high='red',na.value='whitesmoke')+
  facet_wrap(~stock)+theme_bw()+
  theme(legend.position=c(.85,.15))+
  ggtitle("Estimated concurvity among variables by model")+
  ylab("")+xlab("")+guides(fill=guide_legend(title="Concurvity"))
dev.off()


alt_AIC<-keep_AIC_m_big
for(x in 1:nrow(alt_AIC))
  for(y in 2:ncol(alt_AIC))
    alt_AIC[x,y]<--1*(alt_AIC[x,1]-alt_AIC[x,y])

in_dat_t1<-as.data.frame(round(alt_AIC,2))

library(tidyr)
library(scales)
# Specify the columns to be colored
columns_to_color <- colnames(in_dat_t1[2:5])

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
  labs(title = "Change in AIC from null model", x = "",y='') +
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
  labs(title = "Deviance explained", x = NULL, y = NULL) +  # Add title and remove axis labels
  theme_minimal() + 
  xlim(0,2)+
  theme(axis.text.x = element_blank(), # Remove x-axis text
        axis.ticks.x = element_blank(), # Remove x-axis ticks
        axis.text.y = element_blank(), # Remove x-axis text
        axis.ticks.y = element_blank(), # Remove x-axis ticks
        panel.grid = element_blank(),  # Remove grid lines
        legend.position = "none")  # Place the legend on the right

png("plots/model_select.png",height=6,width=7,res=350,units='in') 
aic_table + dev_ex_plot +plot_layout(width=c(2,1))
dev.off()

in_col<-c("#ff5050","#0034c377","#ff505077","#0034c3")
in_col<-c("#ff5050","#0034c377","#ff505077","#0034c3","#3da550","#ff8f38")
in_col2<-c("#ff5050","#0034c377","#ff505077","#0034c3","#3da550","#3da55099","#ff8f38","#ff8f3899")
mort_term$covar[mort_term$covar=="N_imm"]<-"Abundance"
rec_plot<-ggplot(out_plot_r)+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=obs),col='black',lwd=2)+
  geom_ribbon(aes(x=year,ymin=y_dn,ymax=y_up,fill=stock),alpha=0.3,lwd=2)+
  geom_line(aes(x=year,y=preds,col=stock),lwd=2)+
  theme_bw()+ylab("Recruitment")+
  scale_color_manual(values=in_col)+
  scale_fill_manual(values=in_col)+
  facet_wrap(~stock,ncol=1)+
  theme(legend.position='none')

rec_plot_trm<-ggplot(rec_term)+
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
  scale_color_manual(values=in_col)+
  scale_fill_manual(values=in_col)
  
ughh<-out_plot_m[,c(1,5,6)]
colnames(ughh)<-c("est_m","year","stock'")

mod_mort<-unc_mort
mod_mort$up_m[mod_mort$up_m>max(mod_mort$est_m)]<-1.5*max(mod_mort$est_m)

mod_mort_cap <- mod_mort %>%
  group_by(stock) %>%
  mutate(up_m = pmin(up_m, 1.5 * max(est_m, na.rm = TRUE))) %>% # The `pmin` function is used here to take the minimum of the two values element-wise, effectively capping the up_m values.
  ungroup()


mort_plot<-ggplot()+
  geom_point(data=out_plot_m,aes(x=year,y=obs))+
  geom_errorbar(data=filter(mod_mort_cap,year<2023&year!=2020),aes(x=year,y=est_m,ymin=dn_m,ymax=up_m))+
  geom_ribbon(data=out_plot_m,aes(x=year,ymin=y_dn,ymax=y_up,fill=stock),alpha=0.3,lwd=2)+
  geom_line(data=out_plot_m,aes(x=year,y=preds,col=stock),lwd=1.5)+
  theme_bw()+ylab("Mortality")+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)+
  facet_wrap(~stock,ncol=1,scales='free_y')+
  theme(legend.position='none')+expand_limits(y=0)

mort_plot_alt<-ggplot()+
  geom_point(data=out_plot_m,aes(x=year,y=obs))+
  geom_errorbar(data=filter(mod_mort_cap,year<2023&year!=2020),aes(x=year,y=est_m,ymin=dn_m,ymax=up_m))+
  geom_ribbon(data=out_plot_m,aes(x=year,ymin=y_dn,ymax=y_up,fill=stock),alpha=0.3,lwd=2)+
  geom_line(data=out_plot_m,aes(x=year,y=preds,col=stock),lwd=1.5)+
  theme_bw()+ylab("Mortality")+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)+
  facet_wrap(~stock,ncol=1,scales='free_y')+
  theme(legend.position='none')+expand_limits(y=0)

mort_plot_alt_trans<-ggplot()+
  geom_point(data=out_plot_m,aes(x=year,y=1-exp(-obs)))+
  geom_errorbar(data=filter(unc_mort,year<2023&year!=2020),aes(x=year,y=est_m,ymin=1-exp(-dn_m),ymax=1-exp(-up_m)))+
  theme_bw()+ylab("p(mortality)")+
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
  geom_line(aes(x=x,y=y,col=stock),lwd=2)+
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
  facet_wrap(~covar,scales='free_x',ncol=1)+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)

# mortality <- density + temperature + competitor + size
library(patchwork)
png("plots/all_recruits_gam.png",height=8,width=7,res=350,units='in') 
rec_plot + rec_plot_trm
dev.off()
png("plots/all_mort_gam.png",height=8,width=8,res=350,units='in') 
mort_plot + mort_plot_trm
dev.off()
png("plots/all_mort_all.png",height=4,width=8,res=350,units='in') 
mort_plot_1
dev.off()

design<-"112
         112
         112"
png("plots/all_mort_all2.png",height=10,width=6,res=350,units='in') 
mort_plot +mort_plot_1 + plot_layout(nrow=2, design=design)
dev.off()



#==================================
# among pop cors
#========================
unique(outs$process)
library(GGally)
library(forecast)
#==by prGGally#==by process
#==put ccf in upper triangle
use_gg<-filter(outs,process=='Recruitment')
casted<-dcast(use_gg,Year~species,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

# Correlation matrix plot
p2 <- ggcorr(casted, label = TRUE, label_round = 2)
g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p<-ncol(casted)
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}
p1 = ggpairs(casted, lower = list(continuous = my_fn))
png("plots/all_rec_corr.png",height=10,width=10,res=350,units='in') 
print(p1)
dev.off()

#=========================
# other mortality

use_gg<-filter(outs,process%in%c('Other mortality',"M_mat"))
casted<-dcast(use_gg,Year~species,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

# Correlation matrix plot
p2 <- ggcorr(casted, label = TRUE, label_round = 2)
g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p<-ncol(casted)
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}

p1 = ggpairs(casted, lower = list(continuous = my_fn))
png("plots/all_mort_corr.png",height=10,width=10,res=350,units='in') 
print(p1)
dev.off()

# other mortality

use_gg<-filter(outs,process%in%c('Other mortality',"M_imm"))
casted<-dcast(use_gg,Year~species,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

# Correlation matrix plot
p2 <- ggcorr(casted, label = TRUE, label_round = 2)
g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p<-ncol(casted)
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}

p1 = ggpairs(casted, lower = list(continuous = my_fn))
png("plots/all_mort_imm_corr.png",height=10,width=10,res=350,units='in') 
print(p1)
dev.off()

#=========================
# abundance

use_gg<-filter(outs,process%in%c('Abundance'))
casted<-dcast(use_gg,Year~species,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

# Correlation matrix plot
p2 <- ggcorr(casted, label = TRUE, label_round = 2)
g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p<-ncol(casted)
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}

p1 = ggpairs(casted, lower = list(continuous = my_fn))
png("plots/all_abund_corr.png",height=10,width=10,res=350,units='in') 
print(p1)
dev.off()


((10*30) +(15*40)+(12*120)+(36*5)+40)/60




t_con<-read.csv("data/cody_ba.csv")
t_cod<-t_con%>%
  group_by(Year)%>%
  summarize(tot_cons=sum(cons,na.rm=T))
s_con<-read.csv("data/cody_op.csv")
s_cod<-s_con%>%
  group_by(Year)%>%
  summarize(tot_cons=sum(cons,na.rm=T))
t_M<-filter(outs,process%in%c("M_imm","Recruitment")&species=="Tanner")
s_M<-filter(outs,process%in%c("M_imm","Recruitment")&species=="Snow")

t_all<-merge(t_cod,t_M)
s_all<-merge(s_cod,s_M)
in_all<-rbind(t_all,s_all)

png("plots/reum_consumption.png",height=10,width=10,res=350,units='in') 
ggplot(in_all,aes(x=tot_cons,y=values))+
  geom_point()+
  geom_smooth()+
  facet_wrap(species~process)
dev.off()

ccf(filter(s_all,process=="Recruitment")$tot_cons,filter(s_all,process=="Recruitment")$values)
ccf(filter(t_all,process=="Recruitment")$tot_cons,filter(t_all,process=="Recruitment")$values)

ggplot(t_all)+
  geom_line(aes(x=Year,y=tot_cons))

ggplot(filter(t_all,process=="Recruitment"))+
  geom_point(aes(x=tot_cons,y=values))+
  geom_smooth(aes(x=tot_cons,y=values),method='lm')

ggplot(filter(s_all,process=="Recruitment"))+
  geom_point(aes(x=tot_cons,y=values))+
  geom_smooth(aes(x=tot_cons,y=values),method='lm')
