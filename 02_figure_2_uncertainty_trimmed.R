# Trimmed duplicate of 02_figure_2_uncertainty.R for manuscript_publish
# Keeps only code needed for manuscript-included figures and required upstream data objects/files.

#==pull in output from models
#==pull in pictures
#==plot time series with uncertainty
#==call out important periods in time serieslibrary(dplyr)
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


#==add PIBKC, PIRKC
rep_files<-c("/models/bbrkc/test/rkc.rep",
             "/models/pirkc/test/rkc.rep",
             "/models/smbkc/test/bkc.rep",
             "/models/pibkc/test/bkc.rep")

species<-c("BBRKC","PIRKC","SMBKC","PIBKC")
maturity<-c(120,120,105,120)
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

#==get uncertainty in M estimates
cor_files<-c("/models/bbrkc/test/rkc.cor",
             "/models/pirkc/test/rkc.cor",
             "/models/smbkc/test/bkc.cor",
             "/models/pibkc/test/bkc.cor")
input_m<-c(-1.4997,-1.229,-1.6225,-1.6225)
keep_uncertainty_m<-NULL
keep_uncertainty_totn<-NULL
keep_uncertainty_fishn<-NULL
keep_uncertainty_rec<-NULL
keep_uncertainty_numbers_pred<-NULL
keep_uncertainty_f<-NULL
for(x in 1:length(rep_files))
{
  ttt<-readLines(paste(getwd(),cor_files[x],sep=""))
  take_em<-grep('nat_m',ttt)
  yrs<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  for(y in 1:length(take_em))
  {
  zzz<-unlist(strsplit(ttt[take_em[y]],split=' '))
  tmp<-data.frame(stock=species[x],
                  est_m=exp(input_m[x]+as.numeric(zzz[nzchar(zzz)][3])),
                  up_m=exp(input_m[x]+as.numeric(zzz[nzchar(zzz)][3])+as.numeric(zzz[nzchar(zzz)][4])*1.96),
                  dn_m=exp(input_m[x]+as.numeric(zzz[nzchar(zzz)][3])-as.numeric(zzz[nzchar(zzz)][4])*1.96),
                  sd=as.numeric(zzz[nzchar(zzz)][4]),
                  Year=yrs[y])
  keep_uncertainty_m<-rbind(keep_uncertainty_m,tmp)
  }
  
  yrs<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  keep_uncertainty_totn<-rbind(keep_uncertainty_totn,
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
    keep_uncertainty_fishn<-rbind(keep_uncertainty_fishn,tmp)
  }  
  
  yrs<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  keep_uncertainty_rec<-rbind(keep_uncertainty_rec,
                              parse_cor_rows(ttt,'recruits',species[x],yrs,year_offset=-5))
  keep_uncertainty_numbers_pred<-rbind(keep_uncertainty_numbers_pred,
                                       parse_cor_rows(ttt,'numbers_pred',species[x],yrs))
  
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
    keep_uncertainty_f<-rbind(keep_uncertainty_f,tmp)
  }
  
}

#==need a flag for the type of mortality to split for snow + tanner?
#==plot snow and tanner together and then the kind crabs together?
#==snow and tanner with recruitment and immature mortality; fishing and mature mortality?
all_dat<-NULL
all_size_dat<-NULL
keep_rec_mle<-NULL
for(x in 1:length(outs_in))
{
years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
rec_ts<-filter(keep_uncertainty_rec,stock==species[x])
totn_ts<-filter(keep_uncertainty_totn,stock==species[x])
rec_max<-max(rec_ts$estimate,na.rm=TRUE)
totn_max<-max(totn_ts$estimate,na.rm=TRUE)
keep_rec_mle<-rbind(keep_rec_mle,
                    data.frame(species=species[x],
                               Year=rec_ts$Year,
                               rec_mle=rec_ts$estimate,
                               rec_max=rec_max))

#names(outs_in[[x]])
#outs_in[[x]]$pred_pop_num
tmp_ssn<-apply(outs_in[[x]]$pred_pop_num[,which(outs_in[[x]]$sizes>maturity[x])],1,sum)
tmp_ssn<-tmp_ssn/max(tmp_ssn)

cbind(apply(outs_in[[x]]$pred_pop_num,1,sum)/max(apply(outs_in[[x]]$pred_pop_num,1,sum),na.rm=T),
      tmp_ssn)
plot_dat<-data.frame(values=c(rec_ts$estimate/rec_max,
                              outs_in[[x]]$`natural mortality`[,1],
                              outs_in[[x]]$est_fishing_mort,
                              tmp_ssn,
                              totn_ts$estimate/totn_max),
                     Year=c(rec_ts$Year,rep(years,3),totn_ts$Year),
                     process=c(rep("Recruitment",length(rec_ts$Year)),
                               rep("Other mortality",length(years)),
                               rep("Fishing mortality",length(years)),
                               rep("Spawner abundance",length(years)),
                               rep("Abundance",length(totn_ts$Year))))
plot_dat$species<-species[x]
    all_dat<-rbind(all_dat,plot_dat) 
    
yerp<-data.frame(values=c(outs_in[[x]]$"survey selectivity"[1,],
                    outs_in[[x]]$"total_fish_sel",
                    outs_in[[x]]$"ret_fish_sel",
                    outs_in[[x]]$"in_prob_molt"),
           size=rep(outs_in[[x]]$sizes,4),
           process=c(rep("Survey selectivity",length(outs_in[[x]]$sizes)),
                     rep("Total fishery",length(outs_in[[x]]$sizes)),
                     rep("Retained selectivity",length(outs_in[[x]]$sizes)),
                     rep("Molting",length(outs_in[[x]]$sizes))),
           species=species[x])
all_size_dat<-rbind(all_size_dat,yerp) 
}

in_col<-c("#ff5050","#0034c377","#ff505077","#0034c3")
all_dat$values[all_dat$process=='Fishing mortality' & all_dat$values==0]<-NA
uncertainty_m<-data.frame(species=keep_uncertainty_m$stock,
                          Year=keep_uncertainty_m$Year,
                          process="Other mortality",
                          ci_dn=keep_uncertainty_m$dn_m,
                          ci_up=keep_uncertainty_m$up_m)
uncertainty_rec<-merge(keep_rec_mle,
                       keep_uncertainty_rec,
                       by.x=c("species","Year"),
                       by.y=c("stock","Year"))
uncertainty_rec<-data.frame(species=uncertainty_rec$species,
                            Year=uncertainty_rec$Year,
                            process="Recruitment",
                            ci_dn=pmax((uncertainty_rec$rec_mle-1.96*uncertainty_rec$sd)/uncertainty_rec$rec_max,0),
                            ci_up=(uncertainty_rec$rec_mle+1.96*uncertainty_rec$sd)/uncertainty_rec$rec_max)
uncertainty_totn<-merge(keep_uncertainty_totn,
                        aggregate(estimate~stock,data=keep_uncertainty_totn,max,na.rm=TRUE),
                        by="stock",suffixes=c("","_max"))
uncertainty_totn<-data.frame(species=uncertainty_totn$stock,
                             Year=uncertainty_totn$Year,
                             process="Abundance",
                             ci_dn=pmax((uncertainty_totn$estimate-1.96*uncertainty_totn$sd)/uncertainty_totn$estimate_max,0),
                             ci_up=(uncertainty_totn$estimate+1.96*uncertainty_totn$sd)/uncertainty_totn$estimate_max)
uncertainty_f<-data.frame(species=keep_uncertainty_f$stock,
                          Year=keep_uncertainty_f$Year,
                          process="Fishing mortality",
                          ci_dn=keep_uncertainty_f$dn_f,
                          ci_up=keep_uncertainty_f$up_f)
uncertainty_dat<-rbind(uncertainty_m,uncertainty_rec,uncertainty_f,uncertainty_totn)
uncertainty_cap<-aggregate(values~species+process,
                           data=filter(all_dat,process%in%c("Recruitment","Other mortality","Fishing mortality")),
                           max,na.rm=TRUE)
uncertainty_cap$cap<-uncertainty_cap$values*1.25
uncertainty_dat<-merge(uncertainty_dat,uncertainty_cap[,c("species","process","cap")],
                       by=c("species","process"),all.x=TRUE)
uncertainty_dat$ci_up<-pmin(uncertainty_dat$ci_up,uncertainty_dat$cap)
uncertainty_dat<-uncertainty_dat[order(uncertainty_dat$species,uncertainty_dat$process,uncertainty_dat$Year),]
uncertainty_dat$ribbon_group<-ave(uncertainty_dat$Year,
                                  uncertainty_dat$species,
                                  uncertainty_dat$process,
                                  FUN=function(x)cumsum(c(1,diff(x)>1)))
king_proc<-ggplot()+
  geom_ribbon(data=filter(uncertainty_dat,!process%in%c("Spawner abundance","Abundance")),aes(x=Year,ymin=ci_dn,ymax=ci_up,fill=species,group=interaction(species,process,ribbon_group)),alpha=.2,colour=NA)+
  geom_line(data=filter(all_dat,!process%in%c("Spawner abundance","Abundance")),aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),cols=vars(species),scales='free_y')+
  theme_bw()+ylab("")+
  scale_color_manual(values=in_col)+
  scale_fill_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

king_proc_agg<-ggplot()+
  geom_line(data=filter(all_dat,process!="Abundance"),aes(x=Year,y=values,col=species),lwd=1.2)+
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

king_proc_agg_size<-ggplot()+
  geom_line(data=all_size_dat,aes(x=size,y=values,col=species),lwd=1.2)+
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
div_n<-c(max(outs_in[[1]]$n_obs),max(outs_in[[2]]$n_obs),max(outs_in[[3]]$n_obs),max(outs_in[[4]]$n_obs))
ind_dat<-NULL

for(x in 1:length(outs_in))
{
years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
fit_ts<-filter(keep_uncertainty_numbers_pred,stock==species[x])
df_1<-data.frame(pred=outs_in[[x]]$numbers_pred/div_n[x],
                 obs=outs_in[[x]]$n_obs/div_n[x],
                 year=years,
                 ci_dn=(outs_in[[x]]$n_obs/div_n[x]) /  exp(1.96*sqrt(log(1+outs_in[[x]]$sigma_numbers^2))),
                 ci_up=(outs_in[[x]]$n_obs/div_n[x]) *  exp(1.96*sqrt(log(1+outs_in[[x]]$sigma_numbers^2))),
                 pred_ci_dn=pmax((fit_ts$estimate - 1.96 * fit_ts$sd) / div_n[x],0),
                 pred_ci_up=(fit_ts$estimate + 1.96 * fit_ts$sd) / div_n[x])

df_1$species<-species[x]
df_1$color<-in_col[x]
ind_dat<-rbind(ind_dat,df_1)
}

ind_dat$ci_up[ind_dat$ci_up>1.5]<-1.5

king_abnd<-ggplot(data=ind_dat)+
  geom_ribbon(aes(x=year,ymin=pred_ci_dn,ymax=pred_ci_up,fill=species),alpha=.18,colour=NA)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred,col=species),lwd=1.5)+
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
  xlab("")+ guides(color="none",fill="none")


png("plots/fig_2.png",height=6,width=8,res=400,units='in')
print(king_abnd/king_proc + plot_layout(heights = c(2, 3)))
dev.off()


#=========================
# plot spawner abundance vs recruits
#===============================
srr_dat<-filter(all_dat,process%in%c("Recruitment","Spawner abundance"))
in_dat<-data.frame(dcast(data=srr_dat,formula=Year+species~process , value.var = "values"))

king_srr<-ggplot()+
  geom_point(data=in_dat,aes(x=Spawner.abundance,y=Recruitment,col=species),size=2)+
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

library(GGally)
all_dat$species_process<-paste(all_dat$species,"_",substring(all_dat$process,1,1),sep="")
casted<-dcast(all_dat,Year~species_process,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

p1 = ggpairs(casted, lower = list(continuous = my_fn))


# Correlation matrix plot
 p2 <- ggcorr(casted, label = TRUE, label_round = 2)
 g2 <- ggplotGrob(p2)
 colors <- g2$grobs[[6]]$children[[3]]$gp$fill
# 
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
 png("plots/rkc_cors.png",height=13,width=13,res=400,units='in')
 print(p1)
 dev.off()

all_dat_kc<-all_dat

casted<-dcast(filter(all_dat,process%in%c("Abundance","Recruitment")),Year~species_process,value.var="values")
#p1<-ggpairs(casted[,-1])+ theme_grey(base_size = 8)
#p2 <- ggcorr(casted[,-1], label = TRUE)

abund<-filter(all_dat,process%in%c("Abundance"))[,c(1,2,4)]
colnames(abund)<-c("Abundance","Year","Species")
rec<-filter(all_dat,process%in%c("Recruitment"))[,c(1,2,4)]
colnames(rec)<-c("Recruitment","Year","Species")
plot_srr<-merge(rec,abund,by=c("Year"))
plot_srr <- rec %>% right_join(abund, by=c("Year","Species"))
ggplot(plot_srr)+
  geom_point(aes(x=Abundance,y=Recruitment))+
  geom_smooth(aes(x=Abundance,y=Recruitment),span=1)+
  theme_bw()+facet_wrap(~Species)



#===============================================
# get the fits and diagnostic plots for all stocks
#=================================================

for(y in 1:length(species))
{

dat_sel<-data.frame(value=c(outs_in[[y]]$'survey selectivity'[1,]),
                    sizes=rep(outs_in[[y]]$sizes),
                    est=c(rep("Estimated",length(outs_in[[y]]$sizes))))
fish_sel<-data.frame(value=c(outs_in[[y]]$ret_fish_sel,outs_in[[y]]$total_fish_sel),
                     sizes=rep(outs_in[[y]]$sizes,2),
                     Fishery=c(rep("Retained",length(outs_in[[y]]$sizes)),rep("Total",length(outs_in[[y]]$sizes))))
molt<-data.frame(value=c(outs_in[[y]]$'in_prob_molt'),
                 sizes=rep(outs_in[[y]]$sizes),
                 est=c(rep("Estimated",length(outs_in[[y]]$sizes))))

s_sel<-ggplot()+
  geom_line(data=filter(dat_sel,est=="Estimated"),aes(x=sizes,y=value),col=2,lwd=2)+
  theme_bw()+ylab("Selectivity")+xlab("Carapace length (mm)")+
  annotate("text",x=75,y=0.9,label="Survey")+
  ylim(0,1)

molt_pl<-ggplot()+
  geom_line(data=filter(molt,est=="Estimated"),aes(x=sizes,y=value),col=2,lwd=2)+
  theme_bw()+ylab("Molting probability")+xlab("Carapace length (mm)")+
  ylim(0,1)


f_sel<-ggplot()+
  geom_line(data=fish_sel,aes(x=sizes,y=value,col=Fishery),lwd=2)+
  theme_bw()+xlab("Carapace length (mm)")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position=c(.25,.8),
        legend.background = element_blank())+ylim(0,1)+
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
  labs(x="Post-molt carapace length (mm)",y="Pre-molt carapace length (mm)") +
  xlim(min(outs_in[[y]]$sizes),max(outs_in[[y]]$sizes))


png(paste('plots/',species[y],'model_growth.png',sep=''),height=8,width=8,res=350,units='in')
print(p / (s_sel | f_sel |molt_pl) + plot_layout(nrow=2,heights=c(2,1)))
dev.off()

#==============================================
#==plot survey data
#=============================================
years<-seq(outs_in[[y]]$"styr",outs_in[[y]]$"endyr")
sizes<-outs_in[[y]]$"sizes"
obs_comp<-outs_in[[y]]$"n_obs_len"
rownames(obs_comp)<-years
colnames(obs_comp)<-sizes
df_1<-melt(obs_comp)
colnames(df_1)<-c("Year","Size","Proportion")
df_1$quant<-'Observed'

tmp_size<-outs_in[[y]]$'pred numbers at size'
rownames(tmp_size)<-years
colnames(tmp_size)<-sizes
df_3<-melt(tmp_size)
colnames(df_3)<-c("Year","Size","Proportion")
df_3$quant<-'Predicted'

input_size<-rbind(df_3,df_1)

n_size_all<-ggplot(data=input_size,aes(x = (Size), y =Proportion,col=quant)) + 
  geom_line(lwd=1.1)+
  theme_bw()+
  ylab("Proportion")+theme(axis.title=element_text(size=11))+
  xlab("Carapace length (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Year)+ labs(col='') 

png(paste("plots/",species[y],"n_size_comp_all.png",sep=""),height=8,width=8,res=350,units='in') 
print(n_size_all)
dev.off()


#==aggregate
df2<-data.frame(pred=apply(outs_in[[y]]$'pred numbers at size',2,median),
                Size=(sizes))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))

all_size<-ggplot(data=df2,aes(x = factor(Size), y =pred)) + 
  geom_boxplot(data=df_1,aes(x = Size, y =Proportion ),fill='grey') +
  stat_summary(fun.y=mean, geom="line", aes(group=1),lwd=1.5,col='blue',alpha=.8)  + 
  theme_bw()+
  ylab("Proportion")+theme(axis.title=element_text(size=11))+
  xlab("Carapace length (mm)")+
  scale_x_discrete(breaks=seq(32.5,132.5,10))+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


#==fits index
# names(outs_in[[1]])
div_n<-1000000000
df_1<-data.frame(pred=outs_in[[y]]$numbers_pred/div_n,
                 obs=outs_in[[y]]$n_obs/div_n,
                 year=years,
                 ci_dn=(outs_in[[y]]$n_obs/div_n) /  exp(1.96*sqrt(log(1+outs_in[[y]]$surv_n_cv^2))),
                 ci_up=(outs_in[[y]]$n_obs/div_n) *  exp(1.96*sqrt(log(1+outs_in[[y]]$surv_n_cv^2))),
                 recruits=(outs_in[[y]]$recruits)/div_n)

plt_abnd<-ggplot(data=df_1)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred),lwd=1.5,col='blue',alpha=.8)+
  #geom_line(aes(x=year,y=recruits),col='purple',lwd=.7,lty=1,alpha=0.8)+
  theme_bw()+
  scale_x_continuous(position = "top",name='TOTAL')+
  ylab("Abundance (billions)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



png(paste("plots/",species[y],"ind_fits.png",sep=""),height=8,width=5,res=350,units='in') 
print(plt_abnd/all_size)
dev.off()

#==============================================
#==plot catch data
#=============================================
#retained
input<-paste(getwd(),"/models/",species[y],"/test/catch_dat.DAT",sep='')
cat_dat<-readLines(input)
ind<-grep("which years of retained",cat_dat)[1]
indn<-grep("number of years of retained",cat_dat)[1]
n_yr<-scan(input,skip=indn,n=1,quiet=T)
ret_yrs<-scan(input,skip=ind,n=n_yr,quiet=T)

#==fits catch
div_n<-1000000000
trans_ret<-rep(0,length(years)-1)
trans_ret[which(!is.na(match(years,ret_yrs)))]<-outs_in[[y]]$ret_cat_numbers/div_n
df_1<-data.frame(pred=outs_in[[y]]$pred_retained_n/div_n,
                 obs=trans_ret,
                 year=years[-length(years)],
                 ci_dn=(trans_ret) /  exp(1.96*sqrt(log(1+0.05^2))),
                 ci_up=(trans_ret) *  exp(1.96*sqrt(log(1+0.05^2))))

ret_abnd<-ggplot(data=df_1)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred),lwd=1.5,col='blue',alpha=.8)+
  #geom_line(aes(x=year,y=recruits),col='purple',lwd=.7,lty=1,alpha=0.8)+
  theme_bw()+
  scale_x_continuous(position = "top",name='RETAINED')+
  scale_y_continuous()+
  ylab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



if(species[y]%in%c("BBRKC","SMBKC"))
{
  ind<-grep("which years of retained fishery size comp data",cat_dat)[1]
  indn<-grep("number of years of retained fishery size comp data",cat_dat)[1]
  n_yr<-scan(input,skip=indn,n=1,quiet=T)
  ret_yrs<-scan(input,skip=ind,n=n_yr,quiet=T)
tmp_size<-outs_in[[y]]$'obs_retained_size_comp'
rownames(tmp_size)<-ret_yrs
colnames(tmp_size)<-sizes
df_1<-melt(tmp_size)
colnames(df_1)<-c("Year","Size","Proportion")
df_1$quant<-'Observed'

tmp_size<-outs_in[[y]]$'pred_retained_size_comp'
rownames(tmp_size)<-years[-length(years)]
colnames(tmp_size)<-sizes
df_3<-melt(tmp_size)
colnames(df_3)<-c("Year","Size","Proportion")
df_3$quant<-'Predicted'

input_size<-rbind(df_3,df_1)
input_size<-filter(input_size,Proportion!='nan')
input_size$Proportion<-as.numeric(input_size$Proportion)

ret_size_all<-ggplot(data=input_size,aes(x = (Size), y =Proportion,col=quant)) +
  geom_line(lwd=1.1)+
  theme_bw()+
  ylab("Proportion")+theme(axis.title=element_text(size=11))+
  xlab("Carapace length (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Year)+ labs(col='')

png(paste("plots/",species[y],"ret_size_all.png",sep=""),height=8,width=8,res=350,units='in')
print(ret_size_all)
dev.off()

#==aggregate
derp<-outs_in[[y]]$'pred_retained_size_comp'
derp2<-(derp[-which(derp[,1]=='nan'),])
derp3<-matrix(as.numeric(derp2),ncol=ncol(derp))
df2<-data.frame(pred=apply(derp3,2,median),
                Size=(sizes))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))

ret_size<-ggplot() +
  geom_boxplot(data=df_1,aes(x = as.numeric(as.character(Size)), y =Proportion,group=Size ),fill='grey') +
  geom_line(data=df2,aes(x = as.numeric(as.character(Size)), y =as.numeric(as.character(pred))),
            col='blue',lwd=2,alpha=.8)+
  theme_bw()+
  ylab("")+theme(axis.title=element_text(size=11))+
  xlab("Carapace length (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


ind<-grep("total catch years",cat_dat)[1]
indn<-grep("number of years of discard",cat_dat)[1]
n_yr<-scan(input,skip=indn,n=1,quiet=T)
tot_yrs<-scan(input,skip=ind,n=n_yr,quiet=T)

#==fits total catch
div_n<-1000000000
trans_ret<-rep(0,length(years)-1)
trans_ret[which(!is.na(match(years,tot_yrs)))]<-outs_in[[y]]$tot_cat_numbers/div_n
df_1<-data.frame(pred=outs_in[[y]]$pred_tot_n/div_n,
                 obs=trans_ret,
                 year=years[-length(years)],
                 ci_dn=(trans_ret) /  exp(1.96*sqrt(log(1+0.05^2))),
                 ci_up=(trans_ret) *  exp(1.96*sqrt(log(1+0.05^2))))

tot_abnd<-ggplot(data=df_1)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred),lwd=1.5,col='blue',alpha=.8)+
  #geom_line(aes(x=year,y=recruits),col='purple',lwd=.7,lty=1,alpha=0.8)+
  theme_bw()+
  scale_x_continuous(position = "top",name='TOTAL')+
  scale_y_continuous()+
  ylab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# total
tmp_size<-outs_in[[y]]$'obs_tot_size_comp'
rownames(tmp_size)<-tot_yrs
colnames(tmp_size)<-sizes
df_1<-melt(tmp_size)
colnames(df_1)<-c("Year","Size","Proportion")
df_1$quant<-'Observed'

tmp_size<-outs_in[[y]]$'pred_tot_size_comp'
rownames(tmp_size)<-years[-length(years)]
colnames(tmp_size)<-sizes
df_3<-melt(tmp_size)
colnames(df_3)<-c("Year","Size","Proportion")
df_3$quant<-'Predicted'

input_size<-rbind(df_3,df_1)
input_size<-filter(input_size,Proportion!='nan')
input_size$Proportion<-as.numeric(input_size$Proportion)

tot_size_all<-ggplot(data=input_size,aes(x = (Size), y =Proportion,col=quant)) +
  geom_line(lwd=1.1)+
  theme_bw()+
  ylab("Proportion")+theme(axis.title=element_text(size=11))+
  xlab("Carapace length (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Year)+ labs(col='')

png(paste("plots/",species[y],"tot_size_all.png",sep=""),height=8,width=8,res=350,units='in')
print(tot_size_all)
dev.off()

#==aggregate
derp<-outs_in[[y]]$'pred_tot_size_comp'
derp2<-(derp[-which(derp[,1]=='nan'),])
derp3<-matrix(as.numeric(derp2),ncol=ncol(derp))
df2<-data.frame(pred=apply(derp3,2,median),
                Size=(sizes))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))

tot_size<-ggplot() +
  geom_boxplot(data=df_1,aes(x = as.numeric(as.character(Size)), y =Proportion,group=Size ),fill='grey') +
  geom_line(data=df2,aes(x = as.numeric(as.character(Size)), y =as.numeric(as.character(pred))),
            col='blue',lwd=2,alpha=.8)+
  theme_bw()+
  ylab("")+theme(axis.title=element_text(size=11))+
  xlab("Carapace length (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


png(paste("plots/",species[y],"catch_fits.png",sep=""),height=8,width=8,res=350,units='in')
print((ret_abnd/ret_size)|(tot_abnd/tot_size))
dev.off()
}else
{
png(paste("plots/",species[y],"catch_fits.png",sep=""),height=8,width=8,res=350,units='in')
print(ret_abnd)
dev.off()
}
}
