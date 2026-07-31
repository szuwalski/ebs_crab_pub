# Trimmed duplicate of 03_figure_3_uncertainty.R for manuscript_publish
# Keeps only code needed for manuscript-included figures and required upstream data objects/files.

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
rep_files<-c("/models/snow/test/snow_down.rep",
             "/models/tanner/test/tanner.rep")

species<-c("Snow", "Tanner")
outs_in<-list(list())

for(x in 1:length(rep_files))
  outs_in[[x]]<-readList(paste(getwd(),rep_files[x],sep=""))

parse_cor_rows<-function(cor_lines, pattern, stock, years, year_offset=0)
{
  take_em<-grep(pattern,cor_lines)
  if(length(take_em)==0)
    return(data.frame(stock=character(),Year=numeric(),estimate=numeric(),sd=numeric()))
  n_rows<-min(length(take_em),length(years))
  out<-NULL
  for(y in 1:n_rows)
  {
    zzz<-unlist(strsplit(cor_lines[take_em[y]],split=' '))
    zzz<-zzz[nzchar(zzz)]
    tmp<-data.frame(stock=stock,
                    Year=years[y]+year_offset,
                    estimate=as.numeric(zzz[3]),
                    sd=as.numeric(zzz[4]))
    out<-rbind(out,tmp)
  }
  out
}

fit_uncertainty_ts<-function(fit_dat, stock, years, pred, cor_file, var_name)
{
  out<-fit_dat[fit_dat$stock==stock,]
  if(nrow(out)==0)
  {
    warning(paste("No",var_name,"rows found in",cor_file,"- plotting",stock,"fit without .cor uncertainty ribbon."))
    return(data.frame(Year=years,estimate=pred,sd=NA_real_))
  }
  out<-merge(data.frame(Year=years),out,by="Year",all.x=TRUE)
  if(any(is.na(out$sd)))
    warning(paste("Incomplete",var_name,"rows found in",cor_file,"- missing years will plot without fit uncertainty ribbon."))
  out
}

#==get uncertainty in M estimates
cor_files<-c("/models/snow/test/snow_down.cor",
             "/models/tanner/test/tanner.cor")
keep_uncertainty_m_ch<-NULL
keep_uncertainty_totn_ch<-NULL
keep_uncertainty_fishn_ch<-NULL
keep_uncertainty_rec_ch<-NULL
keep_uncertainty_imm_numbers_pred<-NULL
keep_uncertainty_mat_numbers_pred<-NULL
keep_uncertainty_f_ch<-NULL

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
  
  yrs<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  keep_uncertainty_totn_ch<-rbind(keep_uncertainty_totn_ch,
                               parse_cor_rows(ttt,'total_population_n',species[x],yrs))
  
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
    keep_uncertainty_fishn_ch<-rbind(keep_uncertainty_fishn_ch,tmp)
  }
  
  yrs<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  keep_uncertainty_rec_ch<-rbind(keep_uncertainty_rec_ch,
                                 parse_cor_rows(ttt,'recruits',species[x],yrs,year_offset=-5))
  keep_uncertainty_imm_numbers_pred<-rbind(keep_uncertainty_imm_numbers_pred,
                                           parse_cor_rows(ttt,'imm_numbers_pred',species[x],yrs))
  keep_uncertainty_mat_numbers_pred<-rbind(keep_uncertainty_mat_numbers_pred,
                                           parse_cor_rows(ttt,'mat_numbers_pred',species[x],yrs))
  
  take_em<-grep('f_dev',ttt)
  take_em_mu<-grep('log_f',ttt)
  log_f<-0
  log_f_sd<-0
  if(length(take_em_mu)>0)
  {
    zzz<-unlist(strsplit(ttt[take_em_mu[1]],split=' '))
    log_f<-as.numeric(zzz[nzchar(zzz)][3])
    log_f_sd<-as.numeric(zzz[nzchar(zzz)][4])
  }
  yrs_f<-yrs[outs_in[[x]]$est_fishing_mort>0]
  for(y in 1:min(length(take_em),length(yrs_f)))
  {
    zzz<-unlist(strsplit(ttt[take_em[y]],split=' '))
    tmp<-data.frame(stock=species[x],
                    est_f=exp(log_f+as.numeric(zzz[nzchar(zzz)][3])),
                    up_f=exp(log_f+as.numeric(zzz[nzchar(zzz)][3])+sqrt(log_f_sd^2+as.numeric(zzz[nzchar(zzz)][4])^2)*1.96),
                    dn_f=exp(log_f+as.numeric(zzz[nzchar(zzz)][3])-sqrt(log_f_sd^2+as.numeric(zzz[nzchar(zzz)][4])^2)*1.96),
                    sd=sqrt(log_f_sd^2+as.numeric(zzz[nzchar(zzz)][4])^2),
                    Year=yrs_f[y])
    keep_uncertainty_f_ch<-rbind(keep_uncertainty_f_ch,tmp)
  }
}


#==need a flag for the type of mortality to split for snow + tanner?
#==plot snow and tanner together and then the kind crabs together?
#==snow and tanner with recruitment and immature mortality; fishing and mature mortality?
all_dat_imm<-NULL
keep_rec_mle_ch<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  rec_ts<-filter(keep_uncertainty_rec_ch,stock==species[x])
  if(nrow(rec_ts)==0)
  {
    warning(paste("No recruits rows found in",cor_files[x],"- plotting",species[x],"recruitment without .cor uncertainty."))
    rec_ts<-data.frame(stock=species[x],
                       Year=years-5,
                       estimate=outs_in[[x]]$recruits,
                       sd=NA_real_)
  }
  rec_max<-max(rec_ts$estimate,na.rm=TRUE)
  keep_rec_mle_ch<-rbind(keep_rec_mle_ch,
                         data.frame(species=species[x],
                                    Year=rec_ts$Year,
                                    rec_mle=rec_ts$estimate,
                                    rec_max=rec_max))
  plot_dat<-data.frame(values=c(rec_ts$estimate/rec_max,
                                outs_in[[x]]$`natural mortality`[,1]),
                       Year=c(rec_ts$Year,rep(years,1)),
                       process=c(rep("Recruitment",length(rec_ts$Year)),
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
uncertainty_m_imm<-filter(keep_uncertainty_m_ch,grepl("_imm",stock))
uncertainty_m_imm<-data.frame(species=sub("_imm","",uncertainty_m_imm$stock),
                              Year=uncertainty_m_imm$Year,
                              process="Other mortality (imm)",
                              ci_dn=uncertainty_m_imm$dn_m,
                              ci_up=uncertainty_m_imm$up_m)
uncertainty_m_mat<-filter(keep_uncertainty_m_ch,grepl("_mat",stock))
uncertainty_m_mat<-data.frame(species=sub("_mat","",uncertainty_m_mat$stock),
                              Year=uncertainty_m_mat$Year,
                              process="Other mortality (mat)",
                              ci_dn=uncertainty_m_mat$dn_m,
                              ci_up=uncertainty_m_mat$up_m)
uncertainty_rec<-merge(keep_rec_mle_ch,
                       keep_uncertainty_rec_ch,
                       by.x=c("species","Year"),
                       by.y=c("stock","Year"))
uncertainty_rec<-data.frame(species=uncertainty_rec$species,
                            Year=uncertainty_rec$Year,
                            process="Recruitment",
                            ci_dn=pmax((uncertainty_rec$rec_mle-1.96*uncertainty_rec$sd)/uncertainty_rec$rec_max,0),
                            ci_up=(uncertainty_rec$rec_mle+1.96*uncertainty_rec$sd)/uncertainty_rec$rec_max)
uncertainty_f<-data.frame(species=keep_uncertainty_f_ch$stock,
                          Year=keep_uncertainty_f_ch$Year,
                          process="Fishing mortality",
                          ci_dn=keep_uncertainty_f_ch$dn_f,
                          ci_up=keep_uncertainty_f_ch$up_f)
uncertainty_imm<-rbind(uncertainty_rec,uncertainty_m_imm)
uncertainty_mat<-rbind(uncertainty_m_mat,uncertainty_f)
uncertainty_imm_cap<-aggregate(values~species+process,
                               data=filter(all_dat_imm,process%in%c("Recruitment","Other mortality (imm)")),
                               max,na.rm=TRUE)
uncertainty_imm_cap$cap<-uncertainty_imm_cap$values*1.25
uncertainty_imm<-merge(uncertainty_imm,uncertainty_imm_cap[,c("species","process","cap")],
                       by=c("species","process"),all.x=TRUE)
uncertainty_imm$ci_up<-pmin(uncertainty_imm$ci_up,uncertainty_imm$cap)
uncertainty_imm<-uncertainty_imm[order(uncertainty_imm$species,uncertainty_imm$process,uncertainty_imm$Year),]
uncertainty_imm$ribbon_group<-ave(uncertainty_imm$Year,
                                  uncertainty_imm$species,
                                  uncertainty_imm$process,
                                  FUN=function(x)cumsum(c(1,diff(x)>1)))
uncertainty_mat_cap<-aggregate(values~species+process,
                               data=filter(all_dat_mat,process%in%c("Other mortality (mat)","Fishing mortality")),
                               max,na.rm=TRUE)
uncertainty_mat_cap$cap<-uncertainty_mat_cap$values*1.25
uncertainty_mat<-merge(uncertainty_mat,uncertainty_mat_cap[,c("species","process","cap")],
                       by=c("species","process"),all.x=TRUE)
uncertainty_mat$ci_up<-pmin(uncertainty_mat$ci_up,uncertainty_mat$cap)
uncertainty_mat<-uncertainty_mat[order(uncertainty_mat$species,uncertainty_mat$process,uncertainty_mat$Year),]
uncertainty_mat$ribbon_group<-ave(uncertainty_mat$Year,
                                  uncertainty_mat$species,
                                  uncertainty_mat$process,
                                  FUN=function(x)cumsum(c(1,diff(x)>1)))
imm_proc<-ggplot()+
  geom_ribbon(data=uncertainty_imm,aes(x=Year,ymin=ci_dn,ymax=ci_up,fill=species,group=interaction(species,process,ribbon_group)),alpha=.2,colour=NA)+
  geom_line(data=all_dat_imm,aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),cols=vars(species),scales='free_y')+
  theme_bw()+ylab("")+  scale_color_manual(values=in_col)+
  scale_fill_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

mat_proc<-ggplot()+
  geom_ribbon(data=uncertainty_mat,aes(x=Year,ymin=ci_dn,ymax=ci_up,fill=species,group=interaction(species,process,ribbon_group)),alpha=.2,colour=NA)+
  geom_line(data=all_dat_mat,aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),cols=vars(species),scales='free_y')+
  theme_bw()+ylab("")+  scale_color_manual(values=in_col)+
  scale_fill_manual(values=in_col)+
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

#==need to put the CVs in the .REP files and pull here
#==immature indices
div_n<-c(max(outs_in[[1]]$imm_n_obs),max(outs_in[[2]]$imm_n_obs))
ind_dat_imm<-NULL

for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  fit_ts<-fit_uncertainty_ts(keep_uncertainty_imm_numbers_pred,species[x],years,
                             outs_in[[x]]$imm_numbers_pred,cor_files[x],"imm_numbers_pred")
  df_1<-data.frame(pred=outs_in[[x]]$imm_numbers_pred/div_n[x],
                   obs=outs_in[[x]]$imm_n_obs/div_n[x],
                   year=years,
                   ci_dn=(outs_in[[x]]$imm_n_obs/div_n[x]) /  exp(1.96*sqrt(log(1+0.15^2))),
                   ci_up=(outs_in[[x]]$imm_n_obs/div_n[x]) *  exp(1.96*sqrt(log(1+0.15^2))),
                   pred_ci_dn=pmax((fit_ts$estimate - 1.96 * fit_ts$sd) / div_n[x],0),
                   pred_ci_up=(fit_ts$estimate + 1.96 * fit_ts$sd) / div_n[x])
  
  df_1$species<-species[x]
  df_1$color<-in_col[x]
  ind_dat_imm<-rbind(ind_dat_imm,df_1)
}

ind_dat_imm$obs[ind_dat_imm$obs==0]<-NA
imm_abnd<-ggplot(data=ind_dat_imm)+
  geom_ribbon(aes(x=year,ymin=pred_ci_dn,ymax=pred_ci_up,fill=species),alpha=.18,colour=NA)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred,col=species),lwd=1.5,alpha=.8)+
  theme_bw()+
  scale_color_manual(values=in_col)+
  scale_fill_manual(values=in_col)+
  ylab("Relative abundance")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  facet_grid(~species,scales='free_y')+xlim(1970,2025)+
  xlab("")+ guides(color="none",fill="none")+ggtitle("IMMATURE")


#==mature indices
div_n<-c(max(outs_in[[1]]$mat_n_obs),max(outs_in[[2]]$mat_n_obs))
ind_dat_mat<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  fit_ts<-fit_uncertainty_ts(keep_uncertainty_mat_numbers_pred,species[x],years,
                             outs_in[[x]]$mat_numbers_pred,cor_files[x],"mat_numbers_pred")
  df_1<-data.frame(pred=outs_in[[x]]$mat_numbers_pred/div_n[x],
                   obs=outs_in[[x]]$mat_n_obs/div_n[x],
                   year=years,
                   ci_dn=(outs_in[[x]]$mat_n_obs/div_n[x]) /  exp(1.96*sqrt(log(1+0.15^2))),
                   ci_up=(outs_in[[x]]$mat_n_obs/div_n[x]) *  exp(1.96*sqrt(log(1+0.15^2))),
                   pred_ci_dn=pmax((fit_ts$estimate - 1.96 * fit_ts$sd) / div_n[x],0),
                   pred_ci_up=(fit_ts$estimate + 1.96 * fit_ts$sd) / div_n[x])
  
  df_1$species<-species[x]
  df_1$color<-in_col[x]
  ind_dat_mat<-rbind(ind_dat_mat,df_1)
}
#==make an upper limit for plotting
ind_dat_mat$pred_ci_up[ind_dat_mat$pred_ci_up>1.339]<-1.339
ind_dat_mat$obs[ind_dat_mat$obs==0]<-NA
mat_abnd<-ggplot(data=ind_dat_mat)+
  geom_ribbon(aes(x=year,ymin=pred_ci_dn,ymax=pred_ci_up,fill=species),alpha=.18,colour=NA)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred,col=species),lwd=1.5,alpha=.8)+
  theme_bw()+
  scale_color_manual(values=in_col)+
  scale_fill_manual(values=in_col)+
  ylab("Relative abundance")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),axis.line.y=element_blank())+
  facet_grid(~species,scales='free_y')+xlim(1970,2025)+
  xlab("")+ guides(color="none",fill="none")+ggtitle("MATURE")

png("plots/fig_3.png",height=6,width=9,res=400,units='in')
print((imm_abnd/imm_proc) | (mat_abnd/mat_proc))
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
  rec_ts<-filter(keep_uncertainty_rec_ch,stock==species[x])
  totn_ts<-filter(keep_uncertainty_totn_ch,stock==species[x])
  if(nrow(rec_ts)==0)
  {
    warning(paste("No recruits rows found in",cor_files[x],"- plotting",species[x],"recruitment without .cor uncertainty."))
    rec_ts<-data.frame(stock=species[x],
                       Year=years-5,
                       estimate=outs_in[[x]]$recruits,
                       sd=NA_real_)
  }
  rec_max<-max(rec_ts$estimate,na.rm=TRUE)
  totn_max<-max(totn_ts$estimate,na.rm=TRUE)
  plot_dat<-data.frame(values=c(rec_ts$estimate/rec_max,
                                outs_in[[x]]$`natural mortality`[,1],
                                outs_in[[x]]$`mature natural mortality`[,1],
                                outs_in[[x]]$est_fishing_mort,
                                totn_ts$estimate/totn_max,
                                apply(outs_in[[x]]$pred_imm_pop_num,1,sum)/div_n_imm[x],
                                apply(outs_in[[x]]$pred_mat_pop_num,1,sum)/div_n_mat[x]),
                       Year=c(rec_ts$Year,rep(years,3),totn_ts$Year,rep(years,2)),
                       process=c(rep("Recruitment",length(rec_ts$Year)),
                                 rep("M_imm",length(years)),
                                 rep("M_mat",length(years)),
                                 rep("Fishing mortality",length(years)),
                                 rep("Abundance",length(totn_ts$Year)),
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
unique(all_dat$process)
chion_srr_dat<-dcast(filter(all_dat,process%in%c("Recruitment","Spawner abundance")),species+Year~process,value.var="values")
chion_srr_dat$spbio<-as.numeric(chion_srr_dat$`Spawner abundance`)

chion_srr<-ggplot()+
  geom_point(data=chion_srr_dat,aes(x=spbio,y=Recruitment,col=species),size=2)+
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


# # Correlation matrix plot
 p2 <- ggcorr(casted, label = TRUE, label_round = 2)
 g2 <- ggplotGrob(p2)
 colors <- g2$grobs[[6]]$children[[3]]$gp$fill
# 
 ##THESE COLORS NEED TO BE FIXED
# # Change background color to tiles in the upper triangular matrix of plots 
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
# 
 png("plots/chion_cors.png",height=13,width=13,res=400,units='in')
 print(p1)
 dev.off()

all_dat$species_process<-paste(all_dat$species,"_",substring(all_dat$process,1,1),sep="")
out_dat<-rbind(all_dat_kc,all_dat)

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
