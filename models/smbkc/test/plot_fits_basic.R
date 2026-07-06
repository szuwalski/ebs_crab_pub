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


#==add PIBKC, PIRKC
rep_files<-c("/models/smbkc/test/bkc.rep")
species<-"SMBKC"
outs_in<-list(list())

for(x in 1:length(rep_files))
  outs_in[[x]]<-readList(paste(getwd(),rep_files[x],sep=""))


#==need a flag for the type of mortality to split for snow + tanner?
#==plot snow and tanner together and then the kind crabs together?
#==snow and tanner with recruitment and immature mortality; fishing and mature mortality?
all_dat<-NULL
all_size_dat<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  
  plot_dat<-data.frame(values=c(outs_in[[x]]$recruits/max(outs_in[[x]]$recruits),
                                outs_in[[x]]$`natural mortality`[,1],
                                outs_in[[x]]$est_fishing_mort,
                                apply(outs_in[[x]]$pred_pop_num,1,sum)/max(apply(outs_in[[x]]$pred_pop_num,1,sum),na.rm=T)),
                       Year=c(years-5,rep(years,3)),
                       process=c(rep("Recruitment",length(years)),
                                 rep("Other mortality",length(years)),
                                 rep("Fishing mortality",length(years)),
                                 rep("Abundance",length(years))))
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
king_proc<-ggplot()+
  geom_line(data=filter(all_dat,process!="Spawner abundance"),aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),cols=vars(species),scales='free_y')+
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

tmp<-filter(all_dat,process=="Other mortality")
plot(tmp$values~tmp$Year,type='l',xlim=c(1975,2022))

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
  
  n_size_all2<-ggplot(data=input_size,aes(x = (Size), y =Proportion,col=quant)) + 
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
    facet_wrap(~Year)+ labs(col='')+ylim(0,0.2)
  

  if(species[y]=="PIRKC")
  {
    n_size_all_s<-ggplot(data=filter(input_size,Year>1990),aes(x = (Size), y =Proportion,col=quant)) + 
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

  }
  
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
    
    

  }
}

print(king_proc)
print(p / (s_sel | f_sel |molt_pl) + plot_layout(nrow=2,heights=c(2,1)))
print(n_size_all)
print(n_size_all_s)
print(tot_size_all)
print(ret_size_all)
print(plt_abnd/all_size)
print(ret_abnd)
