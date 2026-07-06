library(reshape2)
library(ggplot2)
library(mgcv)  
library(dplyr)
library(ggridges)
library(png)
library(PBSmodelling)
library(patchwork)

annotation_custom2 <-   function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, data){ layer(data = data, 
                                                                                                      stat = StatIdentity, 
                                                                                                      position = PositionIdentity,
                                                                                                      geom = ggplot2:::GeomCustomAnn,
                                                                                                      inherit.aes = TRUE, 
                                                                                                      params = list(grob = grob,xmin = xmin, xmax = xmax,
                                                                                                                    ymin = ymin, ymax = ymax))}


#==color coordinate the estimates with figure 1
rep_files<-c("/models/snow/snow_down.rep",
             "/models/tanner/tanner.rep")

species<-c("Snow", "Tanner")
outs_in<-list(list())

for(x in 1:length(rep_files))
  outs_in[[x]]<-readList(paste(getwd(),rep_files[x],sep=""))

#==get uncertainty in M estimates
cor_files<-c("/models/snow/snow_down.cor",
             "/models/tanner/tanner.cor")
keep_uncertainty_m_ch<-NULL

labs<-c("imm","mat")
for(x in 1:length(cor_files))
{
  ttt<-readLines(paste(getwd(),cor_files[x],sep=""))
  take_em<-grep('nat_m_dev',ttt)
  take_em_mat<-grep('nat_m_mat_dev',ttt)
  take_em_m<-grep('log_m_mu',ttt)
  yrs<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  for(z in 1:2)
  for(y in 1:length(take_em))
  {
    zzz<-unlist(strsplit(ttt[take_em[y]],split=' '))
    if(z==2)
     zzz<-unlist(strsplit(ttt[take_em_mat[y]],split=' '))
    m_mu<-unlist(strsplit(ttt[take_em_m[z]],split=' '))
    m_mu<-as.numeric(m_mu[nzchar(m_mu)][3])
    tmp<-data.frame(stock=paste(species[x],"_",labs[z],sep=""),
                    est_m=exp(m_mu+as.numeric(zzz[nzchar(zzz)][3])),
                    up_m=exp(m_mu+as.numeric(zzz[nzchar(zzz)][3])+as.numeric(zzz[nzchar(zzz)][4])*1.96),
                    dn_m=exp(m_mu+as.numeric(zzz[nzchar(zzz)][3])-as.numeric(zzz[nzchar(zzz)][4])*1.96),
                    sd=as.numeric(zzz[nzchar(zzz)][4]),
                    Year=yrs[y])
    keep_uncertainty_m_ch<-rbind(keep_uncertainty_m_ch,tmp)
  }
  
  take_em<-grep('total_population_n',ttt)
  yrs<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  for(y in 1:length(take_em))
  {
    zzz<-unlist(strsplit(ttt[take_em[y]],split=' '))
    tmp<-data.frame(stock=species[x],
                    tot_n=as.numeric(zzz[nzchar(zzz)][3]),
                    up_m=as.numeric(zzz[nzchar(zzz)][3])+as.numeric(zzz[nzchar(zzz)][4])*1.96,
                    dn_m=as.numeric(zzz[nzchar(zzz)][3])-as.numeric(zzz[nzchar(zzz)][4])*1.96,
                    sd=as.numeric(zzz[nzchar(zzz)][4]),
                    Year=yrs[y])
    keep_uncertainty_totn<-rbind(keep_uncertainty_totn,tmp)
  }
  
  take_em<-grep('fished_population_n',ttt)
  yrs<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  for(y in 1:length(take_em))
  {
    zzz<-unlist(strsplit(ttt[take_em[y]],split=' '))
    tmp<-data.frame(stock=species[x],
                    tot_n=as.numeric(zzz[nzchar(zzz)][3]),
                    up_m=as.numeric(zzz[nzchar(zzz)][3])+as.numeric(zzz[nzchar(zzz)][4])*1.96,
                    dn_m=as.numeric(zzz[nzchar(zzz)][3])-as.numeric(zzz[nzchar(zzz)][4])*1.96,
                    sd=as.numeric(zzz[nzchar(zzz)][4]),
                    Year=yrs[y])
    keep_uncertainty_fishn<-rbind(keep_uncertainty_fishn,tmp)
  }
}


#==need a flag for the type of mortality to split for snow + tanner?
#==plot snow and tanner together and then the kind crabs together?
#==snow and tanner with recruitment and immature mortality; fishing and mature mortality?
all_dat_imm<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  plot_dat<-data.frame(values=c(outs_in[[x]]$recruits/max(outs_in[[x]]$recruits),
                                outs_in[[x]]$`natural mortality`[,1]),
                       Year=c(years-5,rep(years,1)),
                       process=c(rep("Recruitment",length(years)),
                                 rep("Other mortality (imm)",length(years))))
  plot_dat$species<-species[x]
  all_dat_imm<-rbind(all_dat_imm,plot_dat)                
}


all_dat_imm$process<-as.character(all_dat_imm$process)
all_dat_imm$process<-factor(all_dat_imm$process, levels=c("Recruitment", "Other mortality (imm)"))

all_dat_mat<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  plot_dat<-data.frame(values=c(outs_in[[x]]$`mature natural mortality`[,1],
                                outs_in[[x]]$est_fishing_mort),
                       Year=c(rep(years,2)),
                       process=c(rep("Other mortality (mat)",length(years)),
                                 rep("Fishing mortality",length(years))))
  plot_dat$species<-species[x]
  all_dat_mat<-rbind(all_dat_mat,plot_dat)                
}
in_col<-c("#3da550","#ff8f38","#3da550","#ff8f38")
all_dat_mat$process<-as.character(all_dat_mat$process)
all_dat_mat$process<-factor(all_dat_mat$process, levels=c("Other mortality (mat)","Fishing mortality"))
all_dat_mat$values[all_dat_mat$process=='Fishing mortality' & all_dat_mat$values==0]<-NA
imm_proc<-ggplot()+
  geom_line(data=all_dat_imm,aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),cols=vars(species),scales='free_y')+
  theme_bw()+ylab("")+  scale_color_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

mat_proc<-ggplot()+
  geom_line(data=all_dat_mat,aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),cols=vars(species),scales='free_y')+
  theme_bw()+ylab("")+  scale_color_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

chion_proc_agg<-ggplot()+
  geom_line(data=rbind(all_dat_mat,all_dat_imm),aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),scales='free_y')+
  theme_bw()+ylab("")+
  scale_color_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# png("plots/fig_agg_chion.png",height=8,width=5,res=400,units='in')
# print(chion_proc_agg)
# dev.off()


#==need to put the CVs in the .REP files and pull here
#==immature indices
div_n<-c(max(outs_in[[1]]$imm_n_obs),max(outs_in[[2]]$imm_n_obs))
ind_dat_imm<-NULL

for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  df_1<-data.frame(pred=outs_in[[x]]$imm_numbers_pred/div_n[x],
                   obs=outs_in[[x]]$imm_n_obs/div_n[x],
                   year=years,
                   ci_dn=(outs_in[[x]]$imm_n_obs/div_n[x]) /  exp(1.96*sqrt(log(1+0.15^2))),
                   ci_up=(outs_in[[x]]$imm_n_obs/div_n[x]) *  exp(1.96*sqrt(log(1+0.15^2))))
  
  df_1$species<-species[x]
  df_1$color<-in_col[x]
  ind_dat_imm<-rbind(ind_dat_imm,df_1)
}


imm_abnd<-ggplot(data=ind_dat_imm)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred,col=color),lwd=1.5,alpha=.8)+
  theme_bw()+
  scale_color_manual(values=in_col)+
  ylab("Relative abundance")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  facet_grid(~species,scales='free_y')+xlim(1970,2023)+
  xlab("")+ guides(color="none")+ggtitle("IMMATURE")


#==mature indices
div_n<-c(max(outs_in[[1]]$mat_n_obs),max(outs_in[[2]]$mat_n_obs))
ind_dat_mat<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  df_1<-data.frame(pred=outs_in[[x]]$mat_numbers_pred/div_n[x],
                   obs=outs_in[[x]]$mat_n_obs/div_n[x],
                   year=years,
                   ci_dn=(outs_in[[x]]$mat_n_obs/div_n[x]) /  exp(1.96*sqrt(log(1+0.15^2))),
                   ci_up=(outs_in[[x]]$mat_n_obs/div_n[x]) *  exp(1.96*sqrt(log(1+0.15^2))))
  
  df_1$species<-species[x]
  df_1$color<-in_col[x]
  ind_dat_mat<-rbind(ind_dat_mat,df_1)
}

mat_abnd<-ggplot(data=ind_dat_mat)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred,col=color),lwd=1.5,alpha=.8)+
  theme_bw()+
  scale_color_manual(values=in_col)+
  ylab("Relative abundance")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),axis.line.y=element_blank())+
  facet_grid(~species,scales='free_y')+xlim(1970,2023)+
  xlab("")+ guides(color="none")+ggtitle("MATURE")

png("plots/fig_3.png",height=6,width=9,res=400,units='in')
(imm_abnd/imm_proc) | (mat_abnd/mat_proc)
dev.off()


#============================================
# correlations between time series
# need to pull spawning biomass in there too
#============================================
div_n_imm<-c(max(apply(outs_in[[1]]$pred_imm_pop_num,1,sum)),max(apply(outs_in[[2]]$pred_imm_pop_num,1,sum)))
div_n_mat<-c(max(apply(outs_in[[1]]$pred_mat_pop_num,1,sum)),max(apply(outs_in[[2]]$pred_mat_pop_num,1,sum)))

all_dat<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  plot_dat<-data.frame(values=c(outs_in[[x]]$recruits/max(outs_in[[x]]$recruits),
                                outs_in[[x]]$`natural mortality`[,1],
                                outs_in[[x]]$`mature natural mortality`[,1],
                                outs_in[[x]]$est_fishing_mort,
                                apply(outs_in[[x]]$pred_imm_pop_num,1,sum)/div_n_imm[x],
                                apply(outs_in[[x]]$pred_mat_pop_num,1,sum)/div_n_mat[x]),
                       Year=c(years-5,rep(years,5)),
                       process=c(rep("Recruitment",length(years)),
                                 rep("M_imm",length(years)),
                                 rep("M_mat",length(years)),
                                 rep("Fishing mortality",length(years)),
                                 rep("N_imm",length(years)),
                                 rep("Spawner abundance",length(years))))
  plot_dat$species<-species[x]
  all_dat<-rbind(all_dat,plot_dat)                
}

library(GGally)
casted<-dcast(all_dat,Year~species+process,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

p1 = ggpairs(casted, lower = list(continuous = my_fn))


#============================================
# SR relationsihp
#============================================

chion_srr_dat<-dcast(filter(all_dat,process%in%c("Recruitment","Spawner abundance")),species+Year~process,value.var="values")
chion_srr<-ggplot()+
  geom_point(data=chion_srr_dat,aes(x="Spawner abundance",y=Recruitment,col=species),size=2)+
  facet_wrap(~species)+
  theme_bw()+ylab("")+
  scale_color_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ylab("Recruits")+xlab("Relative spawner abundance")+
  expand_limits(x=0)

png("plots/chion_srr.png",height=4,width=6,res=400,units='in')
print(chion_srr)
dev.off()

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

png("plots/chion_cors.png",height=13,width=13,res=400,units='in')
print(p1)
dev.off()

all_dat$species_process<-paste(all_dat$species,"_",substring(all_dat$process,1,1),sep="")
out_dat<-rbind(all_dat_kc,all_dat)
write.csv(out_dat,"data/all_output.csv")

unc_m_est<-rbind(keep_uncertainty_m,keep_uncertainty_m_ch)
write.csv(unc_m_est,"data/uncertainty_mort.csv")

#==============================================
#==plot model fits
#=============================================
for(y in 1:length(species))
{
# imm
  years<-seq(outs_in[[y]]$"styr",outs_in[[y]]$"endyr")
  sizes<-outs_in[[y]]$"sizes"
  obs_comp<-outs_in[[y]]$"obs_imm_n_size"
  oioi<-sweep(obs_comp,1,apply(obs_comp,1,sum),FUN="/")
  rownames(oioi)<-years
  colnames(oioi)<-sizes
  df_1<-melt(oioi)
colnames(df_1)<-c("Year","Size","Proportion")
df_1$quant<-'Observed'

tmp_size<-outs_in[[y]]$'immature numbers at size'
rownames(tmp_size)<-years
colnames(tmp_size)<-sizes
df_3<-melt(tmp_size)
colnames(df_3)<-c("Year","Size","Proportion")
df_3$quant<-'Predicted'

input_size<-rbind(df_3,df_1)

imm_size_all<-ggplot(data=input_size,aes(x = (Size), y =Proportion,col=quant)) + 
  geom_line(lwd=1.1)+
  theme_bw()+
  ylab("Proportion")+theme(axis.title=element_text(size=11))+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Year)+ labs(col='') 

#==size comp residuals here
bubble_plot<-function(obs_sc,pred_sc)
{
  diffs<-(pred_sc-obs_sc)/(obs_sc) 
  diffs[diffs=="NaN"]<-0
  diffs[diffs=="Inf"]<-0
  diffs[diffs>5]<-5
  rownames(diffs)<-rownames(diffs)
  colnames(diffs)<-colnames(diffs)
  out_mat<-melt(diffs)
  out_mat$sign<-sign(out_mat$value)
  out_mat$value<-abs(out_mat$value)
  colnames(out_mat)<-c("Year","Size","value","sign")
  
  obs_mat<-melt(obs_sc)
  colnames(obs_mat)<-c("Year","Size","obs")
  
  in_mat<-merge(out_mat,obs_mat,by=c("Year","Size"))
  in_mat$obs<-in_mat$obs/median(in_mat$obs,na.rm=T)
  out_plot<-ggplot(in_mat)+
    geom_point(aes(y=Size,x=Year,size=value,col=as.factor(sign),alpha=obs))+
    theme_bw()+xlab("Year")+ylab("Carapace width (mm)")+
    labs(col="Error direction",size="Relative error")+guides(alpha='none')
  # out_plot<-ggplot(in_mat)+
  #   geom_point(aes(y=Size,x=Year,size=value,col=as.factor(sign)),alpha=.9)+
  #   theme_bw()+xlab("Year")+ylab("Carapace width (mm)")+
  #   labs(col="Residual direction",size="Residual size")
  return(out_plot)
}

png(paste('plots/',species[y],'imm_size_comp_resids.png',sep=''),height=5,width=5,res=350,units='in')
bubble_plot(obs_sc=oioi,pred_sc=tmp_size)
dev.off()

png(paste('plots/',species[y],'imm_size_comp_all.png',sep=''),height=8,width=8,res=350,units='in')
print(imm_size_all)
dev.off()

#==aggregate
df2<-data.frame(pred=apply(outs_in[[y]]$'immature numbers at size',2,median),
                Size=(outs_in[[y]]$sizes))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))

imm_size<-ggplot(data=df2,aes(x = as.numeric(as.character(Size)), y =pred)) + 
  geom_boxplot(data=df_1,aes(x = as.numeric(as.character(Size)), y =Proportion,group=as.numeric(as.character(Size)) ),fill='grey') +
  stat_summary(fun.y=mean, geom="line", aes(group=1),lwd=1.5,col='blue',alpha=.8)  + 
  theme_bw()+
  ylab("Proportion")+theme(axis.title=element_text(size=11))+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

#==mat
years<-seq(outs_in[[y]]$"styr",outs_in[[y]]$"endyr")
sizes<-outs_in[[y]]$"sizes"
obs_comp<-outs_in[[y]]$"obs_mat_n_size"
oioi<-sweep(obs_comp,1,apply(obs_comp,1,sum),FUN="/")
rownames(oioi)<-years
colnames(oioi)<-sizes
df_1<-melt(oioi)

colnames(df_1)<-c("Year","Size","Proportion")
df_1$quant<-'Observed'

tmp_size<-outs_in[[y]]$'mature numbers at size'
rownames(tmp_size)<-years
colnames(tmp_size)<-sizes
df_3<-melt(tmp_size)
colnames(df_3)<-c("Year","Size","Proportion")
df_3$quant<-'Predicted'

input_size<-rbind(df_3,df_1)

mat_size_all<-ggplot(data=input_size,aes(x = (Size), y =Proportion,col=quant)) + 
  geom_line(lwd=1.1)+
  theme_bw()+
  ylab("Size composition")+theme(axis.title=element_text(size=11))+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Year)+ labs(col='') 

png(paste('plots/',species[y],'mat_size_comp_all.png',sep=''),height=8,width=8,res=350,units='in')
print(mat_size_all)
dev.off()

png(paste('plots/',species[y],'mat_size_comp_resids.png',sep=''),height=5,width=5,res=350,units='in')
bubble_plot(obs_sc=oioi,pred_sc=tmp_size)
dev.off()

df2<-data.frame(pred=apply(outs_in[[y]]$'mature numbers at size',2,median),
                Size=(sizes))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))


mat_size<-ggplot(data=df2,aes(x = as.numeric(as.character(Size)), y =pred)) + 
  geom_boxplot(data=df_1,aes(x = as.numeric(as.character(Size)), y =Proportion,group=as.numeric(as.character(Size)) ),fill='grey') +
  stat_summary(fun.y=mean, geom="line", aes(group=1),lwd=1.5,col='blue',alpha=.8)  + 
  theme_bw()+
  ylab("")+
  scale_y_continuous(position = "right",limits=c(0,0.3))+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

#==fits index
# imm
div_n<-1000000000
df_1<-data.frame(pred=outs_in[[y]]$imm_numbers_pred/div_n,
                 obs=outs_in[[y]]$imm_n_obs/div_n,
                 year=years,
                 ci_dn=(outs_in[[y]]$imm_n_obs/div_n) /  exp(1.96*sqrt(log(1+outs_in[[y]]$imm_cv^2))),
                 ci_up=(outs_in[[y]]$imm_n_obs/div_n) *  exp(1.96*sqrt(log(1+outs_in[[y]]$imm_cv^2))),
                 recruits=(outs_in[[y]]$recruits)/div_n)

#df_1<-df_1[-32,]

imm_abnd<-ggplot()+
  geom_segment(data=df_1,aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(data=filter(df_1,year!=2020),aes(x=year,y=obs))+
  geom_line(data=df_1,aes(x=year,y=pred),lwd=1.5,col='blue',alpha=.8)+
  #geom_line(aes(x=year,y=recruits),col='purple',lwd=.7,lty=1,alpha=0.8)+
  theme_bw()+
  scale_x_continuous(position = "top",name='IMMATURE')+
  ylab("Abundance (billions)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# mature
df_1<-data.frame(pred=outs_in[[y]]$mat_numbers_pred/div_n,
                 obs=outs_in[[y]]$mat_n_obs/div_n,
                 year=years,
                 ci_dn=(outs_in[[y]]$mat_n_obs/div_n) /  exp(1.96*sqrt(log(1+outs_in[[y]]$mat_cv^2))),
                 ci_up=(outs_in[[y]]$mat_n_obs/div_n) *  exp(1.96*sqrt(log(1+outs_in[[y]]$mat_cv^2))))
#df_1<-df_1[-32,]

mat_abnd<-ggplot()+
  geom_segment(data=df_1,aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(data=filter(df_1,year!=2020),aes(x=year,y=obs))+
  geom_line(data=df_1,aes(x=year,y=pred),lwd=1.5,col='blue',alpha=.8)+
  theme_bw()+
  scale_x_continuous(position = "top",name='MATURE')+
  ylab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

png(paste('plots/',species[y],'ind_fits.png',sep=''),height=8,width=8,res=350,units='in')
print((imm_abnd/imm_size)|(mat_abnd/mat_size))
dev.off()


#==============================================
#==plot catch data
#=============================================
# retained
tmp_size<-outs_in[[y]]$'obs_retained_size_comp'
if(species[y]=="Snow")
  rownames(tmp_size)<-outs_in[[y]]$ret_cat_yrs
if(species[y]=="Tanner")
 rownames(tmp_size)<-outs_in[[y]]$ret_cat_size_yrs
colnames(tmp_size)<-sizes
df_1<-melt(tmp_size)
colnames(df_1)<-c("Year","Size","Proportion")
df_1$quant<-'Observed'

tmp_size1<-outs_in[[y]]$'pred_retained_size_comp'
rownames(tmp_size1)<-years[-length(years)]
colnames(tmp_size1)<-sizes
tmp_size1[tmp_size1=='nan']<-0
for(x in 1:ncol(tmp_size1))
  tmp_size1[,x]<-as.numeric(tmp_size1[,x])
df_3<-melt(tmp_size1)
colnames(df_3)<-c("Year","Size","Proportion")
df_3$quant<-'Predicted'

input_size<-rbind(df_3,df_1)
input_size$Proportion<-as.numeric(as.character(input_size$Proportion))

ret_size_all<-ggplot(data=input_size,aes(x = (Size), y =Proportion,col=quant)) + 
  geom_line(lwd=1.1)+
  theme_bw()+
  ylab("Proportion")+theme(axis.title=element_text(size=11))+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Year)+ labs(col='') 

# png(paste('plots/',species[y],'ret_size_comp_resids.png',sep=''),height=8,width=8,res=350,units='in')
# bub_dat<-tmp_size1[-which(tmp_size1[,1]=="0"),]
# so_annoy<-matrix(ncol=ncol(bub_dat),nrow=nrow(bub_dat))
# for(x in 1:ncol(so_annoy))
#   so_annoy[,x]<-as.numeric(bub_dat[,x])
# bubble_plot(obs_sc=tmp_size,pred_sc=so_annoy)
# dev.off()

png(paste('plots/',species[y],'ret_size_all.png',sep=''),height=8,width=8,res=350,units='in')
print(ret_size_all)
dev.off()

#==aggregate
df2<-filter(df_3,Proportion!=0)%>%
  group_by(Size)%>%
  summarize(pred=median(as.numeric(Proportion),na.rm=T))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))

ret_size<-ggplot() + 
  geom_boxplot(data=df_1,aes(x = as.numeric(as.character(Size)), y =Proportion,group=Size ),fill='grey') +
  geom_line(data=df2,aes(x = as.numeric(as.character(Size)), y =as.numeric(as.character(pred))),
            col='blue',lwd=2,alpha=.8)+
  theme_bw()+
  ylab("")+theme(axis.title=element_text(size=11))+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

names(outs_in[[y]])
# discards

if(species[y]=="Snow")
{
  tmp_size<-outs_in[[y]]$'obs_discard_size_comp'
  rownames(tmp_size)<-outs_in[[y]]$ret_cat_yrs
}
if(species[y]=="Tanner")
{
  tmp_size<-outs_in[[y]]$'obs_tot_size_comp'
  rownames(tmp_size)<-outs_in[[y]]$tot_cat_size_yrs
}
colnames(tmp_size)<-sizes
df_1<-melt(tmp_size)
colnames(df_1)<-c("Year","Size","Proportion")
df_1$quant<-'Observed'


if(species[y]=="Snow")
{
  tmp_size1<-outs_in[[y]]$'pred_discard_size_comp'
  rownames(tmp_size1)<-years[-length(years)]
  colnames(tmp_size1)<-sizes
}
if(species[y]=="Tanner")
{
  tmp_size1<-outs_in[[y]]$'pred_tot_size_comp'
  rownames(tmp_size1)<-years[-length(years)]
  colnames(tmp_size1)<-sizes
}


tmp_size1[tmp_size1=='nan']<-0
for(x in 1:ncol(tmp_size1))
  tmp_size1[,x]<-as.numeric(tmp_size1[,x])
df_3<-melt(tmp_size1)
colnames(df_3)<-c("Year","Size","Proportion")
df_3$quant<-'Predicted'

# png(paste('plots/',species[y],'disc_size_comp_resids.png',sep=''),height=8,width=8,res=350,units='in')
# bub_dat<-tmp_size1[-which(tmp_size1[,1]=="0"),]
# so_annoy<-matrix(ncol=ncol(bub_dat),nrow=nrow(bub_dat))
# for(x in 1:ncol(so_annoy))
#   so_annoy[,x]<-as.numeric(bub_dat[,x])
# bubble_plot(obs_sc=tmp_size,pred_sc=so_annoy)
# dev.off()

input_size<-rbind(df_3,df_1)
input_size$Proportion<-as.numeric(as.character(input_size$Proportion))
disc_size_all<-ggplot(data=input_size,aes(x = (Size), y =Proportion,col=quant)) + 
  geom_line(lwd=1.1)+
  theme_bw()+
  ylab("Proportion")+theme(axis.title=element_text(size=11))+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Year)+ labs(col='') 

png(paste('plots/',species[y],'disc_size_all.png',sep=''),height=8,width=8,res=350,units='in')
print(disc_size_all)
dev.off()

#==aggregate
if(species[y]=="Snow")
 df2<-data.frame(pred=apply(outs_in[[y]]$'pred_discard_size_comp',2,median),
                Size=(sizes))
if(species[y]=="Tanner")
  df2<-filter(df_3,Proportion!=0)%>%
  group_by(Size)%>%
  summarize(pred=median(as.numeric(Proportion),na.rm=T))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))

disc_size<-ggplot() + 
  geom_boxplot(data=df_1,aes(x = as.numeric(as.character(Size)), y =Proportion,group=Size ),fill='grey') +
  geom_line(data=df2,aes(x = as.numeric(as.character(Size)), y =as.numeric(as.character(pred))),
            col='blue',lwd=2,alpha=.8)+
  theme_bw()+
  ylab("")+theme(axis.title=element_text(size=11))+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(position = "right")

#==fits index
# retained
div_n<-1000000000
trans_ret<-rep(0,length(years)-1)
trans_ret[which(!is.na(match(years,outs_in[[y]]$ret_cat_yrs)))]<-outs_in[[y]]$ret_cat_numbers/div_n
df_1<-data.frame(pred=outs_in[[y]]$pred_retained_n/div_n,
                 obs=trans_ret,
                 year=years[-length(years)],
                 ci_dn=(trans_ret) /  exp(1.96*sqrt(log(1+0.05^2))),
                 ci_up=(trans_ret) *  exp(1.96*sqrt(log(1+0.05^2))))

ret_abnd<-ggplot(data=df_1)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred),lwd=1.5,col='blue',alpha=.8)+
  theme_bw()+
  scale_x_continuous(position = "top",name='RETAINED')+
  scale_y_continuous()+
  ylab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# discards
if(species[y]=="Snow")
{
trans_disc<-rep(0,length(years)-1)
trans_disc[which(!is.na(match(years,outs_in[[y]]$ret_cat_yrs)))]<-outs_in[[y]]$disc_cat_numbers/div_n
df_1<-data.frame(pred=outs_in[[y]]$pred_discard_n/div_n,
                 obs=trans_disc,
                 year=years[-length(years)],
                 ci_dn=(trans_disc) /  exp(1.96*sqrt(log(1+0.07^2))),
                 ci_up=(trans_disc) *  exp(1.96*sqrt(log(1+0.07^2))))

disc_abnd<-ggplot(data=df_1)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred),lwd=1.5,col='blue',alpha=.8)+
  theme_bw()+
  scale_x_continuous(position = "top",name='DISCARD')+
  ylab("Abundance (billions)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
}

if(species[y]=="Tanner")
{
  trans_disc<-rep(0,length(years)-1)
  trans_disc[which(!is.na(match(years,outs_in[[y]]$tot_cat_yrs)))]<-outs_in[[y]]$tot_cat_numbers/div_n
  df_1<-data.frame(pred=outs_in[[y]]$pred_tot_n/div_n,
                   obs=trans_disc,
                   year=years[-length(years)],
                   ci_dn=(trans_disc) /  exp(1.96*sqrt(log(1+0.07^2))),
                   ci_up=(trans_disc) *  exp(1.96*sqrt(log(1+0.07^2))))
  
  disc_abnd<-ggplot(data=df_1)+
    geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
    geom_point(aes(x=year,y=obs))+
    geom_line(aes(x=year,y=pred),lwd=1.5,col='blue',alpha=.8)+
    theme_bw()+
    scale_x_continuous(position = "top",name='DISCARD')+
    ylab("Abundance (billions)")+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
}
png(paste('plots/',species[y],'catch_fits.png',sep=''),height=8,width=8,res=350,units='in')
print((ret_abnd/ret_size)|(disc_abnd/disc_size))
dev.off()


##################################
# plot estimated and specified processes
###################################

dat_sel<-data.frame(value=c(outs_in[[y]]$'survey selectivity'[1,]),
                    sizes=rep(outs_in[[y]]$sizes),
                    est=c(rep("Estimated",length(outs_in[[y]]$sizes))),
                    era=c(rep("1982-present",length(outs_in[[y]]$sizes))))

if(y==2) # tanner
{
  dat_sel1<-data.frame(value=c(outs_in[[y]]$'survey selectivity'[1,]),
                      sizes=rep(outs_in[[y]]$sizes),
                      est=c(rep("Estimated",length(outs_in[[y]]$sizes))),
                      era=c(rep("1975-1981",length(outs_in[[y]]$sizes)))) 
  dat_sel2<-data.frame(value=c(outs_in[[y]]$'survey selectivity'[30,]),
                       sizes=rep(outs_in[[y]]$sizes),
                       est=c(rep("Estimated",length(outs_in[[y]]$sizes))),
                       era=c(rep("1982-present",length(outs_in[[y]]$sizes))))   
  dat_sel<-rbind(dat_sel1,dat_sel2)
}


fish_sel<-data.frame(value=c(outs_in[[y]]$ret_fish_sel[1,],outs_in[[y]]$total_fish_sel[1,]),
                     sizes=rep(outs_in[[y]]$sizes,2),
                     Fishery=c(rep("Retained",length(outs_in[[y]]$sizes)),rep("Total",length(outs_in[[y]]$sizes))))

t_molt<-outs_in[[y]]$'prob_term_molt'
colnames(t_molt)<-outs_in[[y]]$sizes
rownames(t_molt)<-seq(outs_in[[y]]$styr,outs_in[[y]]$endyr)
molt<-melt(t_molt)
colnames(molt)<-c("Year","Size","Probability")

s_sel<-ggplot()+
  geom_line(data=filter(dat_sel,est=="Estimated"),aes(x=sizes,y=value,col=era),lwd=2)+
  theme_bw()+ylab("Selectivity")+xlab("Carapace width (mm)")+
  annotate("text",x=75,y=0.9,label="Survey")+
  ylim(0,1)

molt_pl<-ggplot()+
  geom_line(data=filter(molt),aes(x=Size,y=Probability,col=Year,group=Year),lwd=1.2,alpha=.5)+
  theme_bw()+ylab("Molting probability")+xlab("Carapace width (mm)")+
  ylim(0,1)+theme(legend.position='none')


f_sel<-ggplot()+
  geom_line(data=fish_sel,aes(x=sizes,y=value,col=Fishery),lwd=2)+
  theme_bw()+xlab("Carapace width (mm)")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position=c(.25,.8))+ylim(0,1)+
  expand_limits(y=0)

indat<-outs_in[[y]]$size_trans
colnames(indat)<-outs_in[[y]]$sizes
rownames(indat)<-outs_in[[y]]$sizes

indat<-melt(indat)
indat1<-data.frame("Premolt"=as.numeric(indat$Var1),
                   "Postmolt"=as.numeric((indat$Var2)),
                   "Increment"=as.numeric(indat$value))

p <- ggplot(dat=indat1) 
p <- p + geom_density_ridges(aes(x=Premolt, y=Postmolt, height = Increment,
                                 group = Postmolt, 
                                 alpha=.9999),fill='blue',stat = "identity") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90)) +
  labs(x="Post-molt carapace width (mm)",y="Pre-molt carapace width (mm)") +
  xlim(min(outs_in[[y]]$sizes),max(outs_in[[y]]$sizes))


png(paste('plots/',species[y],'model_growth.png',sep=''),height=8,width=8,res=350,units='in')
print(p / (s_sel | f_sel |molt_pl) + plot_layout(nrow=2,heights=c(2,1)))
dev.off()




}

