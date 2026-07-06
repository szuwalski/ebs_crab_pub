library(reshape2)
library(ggplot2)
library(mgcv)  
library(dplyr)
library(ggridges)
library(png)
library(PBSmodelling)
library(patchwork)




norm_tot_n<-keep_uncertainty_totn
unq_st<-unique(norm_tot_n$stock)
for(x in 1:length(unq_st))
{
 indi<-which(norm_tot_n$stock==unq_st[x])  
 norm_tot_n[indi,2:4]<- norm_tot_n[indi,2:4]/max(norm_tot_n[indi,2])
}

rep_files<-c("/models/bbrkc/test/rkc.rep",
             "/models/pirkc/test/rkc.rep",
             "/models/smbkc/test/bkc.rep",
             "/models/pibkc/test/bkc.rep")

species<-c("BBRKC","PIRKC","SMBKC","PIBKC")
outs_in<-list(list())
tot_abn<-data.frame(year=NULL,
                    stock=NULL,
                    abund=NULL)
for(x in 1:length(rep_files))
{
  outs_in[[x]]<-readList(paste(getwd(),rep_files[x],sep=""))
  fishable_n<-apply(sweep(outs_in[[x]]$pred_pop_num,2,outs_in[[x]]$ret_fish_sel,FUN="*"),1,sum)
  tmp_df<-data.frame(year=seq(outs_in[[x]]$styr,outs_in[[x]]$endyr),
                     stock=species[x],
                     abund=apply(outs_in[[x]]$pred_pop_num,1,sum),
                     other_m=outs_in[[x]]$"natural mortality"[,1],
                     fish_m=outs_in[[x]]$"est_fishing_mort",
                     fish_n=fishable_n,
                     norm_fish_n=fishable_n/max(fishable_n))
  tot_abn<-rbind(tot_abn,tmp_df)
}
                     
rep_files<-c("/models/snow/test/snow_down.rep",
             "/models/tanner/test/tanner.rep")
species<-c("Snow","Tanner")
outs_in<-list(list())
for(x in 1:length(rep_files))
{
  outs_in[[x]]<-readList(paste(getwd(),rep_files[x],sep=""))
  fishable_n<-apply(outs_in[[x]]$pred_mat_pop_num*outs_in[[x]]$ret_fish_sel,1,sum)
  tmp_df<-data.frame(year=seq(outs_in[[x]]$styr,outs_in[[x]]$endyr),
                     stock=species[x],
                     abund=apply(outs_in[[x]]$pred_mat_pop_num + outs_in[[x]]$pred_imm_pop_num ,1,sum),
                     other_m=outs_in[[x]]$"mature natural mortality"[,1],
                     fish_m=c(outs_in[[x]]$"est_fishing_mort"),
                     fish_n=fishable_n,
                     norm_fish_n=fishable_n/max(fishable_n))
  tot_abn<-rbind(tot_abn,tmp_df)
}


tot_abn$prop_fish<-tot_abn$fish_n/tot_abn$abund


plot(outs_in[[x]]$ret_fish_sel[1,],type='b')
lines(outs_in[[x]]$ret_fish_sel[31,])

filter(tot_abn,stock=="Snow")
tot_abn$Fishing<-NA
tot_abn$Fishing[tot_abn$stock%in%c("Snow","BBRKC","Tanner")]<-"Regularly fished"
tot_abn$Fishing[tot_abn$stock%in%c("PIRKC","PIBKC","SMBKC")]<-"Rarely fished"
#==pull total abundance for each species
ggplot(tot_abn)+
  geom_line(aes(x=year,y=abund))+
  theme_bw()+
  facet_wrap(~stock,scales='free_y')


ggplot(tot_abn)+
  geom_line(aes(x=year,y=fish_n))+
  theme_bw()+
  facet_wrap(~stock,scales='free_y')


in_col<-c("#ff5050","#0034c377","#ff505077","#0034c3","#3da550","#ff8f38")

prop_fish<-ggplot(tot_abn)+
  geom_line(aes(x=year,y=prop_fish,group=stock,col=stock),lwd=1.2)+
  theme_bw()+
  facet_wrap(~Fishing,scales='free_y',ncol=1)+
  scale_color_manual(values=in_col)+ylab("Proportion of population abundance vulnerable to fishery")+
  theme(legend.position='none')+ylim(0,1)

prop_fish2<-ggplot(tot_abn)+
  geom_line(aes(x=year,y=prop_fish,group=stock,col=stock),lwd=1.2)+
  theme_bw()+
  facet_wrap(~Fishing,scales='free_y',ncol=1)+
  scale_color_manual(values=in_col)+ylab("Proportion of population abundance vulnerable to fishery")+
  theme(legend.position=c(.5,.3))+ylim(0,1)


unq_st<-unique(keep_uncertainty_fishn$stock)
stoopid<-keep_uncertainty_fishn
for(x in 1:length(unq_st))
 stoopid[which(keep_uncertainty_fishn$stock==unq_st[x]),2:4]<- stoopid[which(keep_uncertainty_fishn$stock==unq_st[x]),2:4]/max(stoopid[which(keep_uncertainty_fishn$stock==unq_st[x]),2]) 

stoopid$Fishing<-NA
stoopid$Fishing[stoopid$stock%in%c("Snow","BBRKC","Tanner")]<-"Regularly fished"
stoopid$Fishing[stoopid$stock%in%c("PIRKC","PIBKC","SMBKC")]<-"Rarely fished" 

stoopid$up_m[stoopid$up_m>1.2]<-1.2
stoopid$dn_m[stoopid$dn_m<0.01]<-0.01

fishable_n<-ggplot()+
  geom_ribbon(data=stoopid,aes(x=Year,ymin=dn_m,ymax=up_m,fill=stock),alpha=.3)+
  geom_line(data=stoopid,aes(x=Year,y=tot_n,col=stock),alpha=.7)+
  geom_point(data=stoopid,aes(x=Year,y=tot_n,col=stock),alpha=.7,size=.51)+
  theme_bw()+
  facet_wrap(~Fishing,scales='free_y',ncol=1)+
  scale_color_manual(values=in_col)+ylab("Relative abundance")+ 
  scale_fill_manual(values=in_col)

png("plots/rel_abn_prop_fish.png",height=4,width=7,res=400,units='in')
prop_fish + fishable_n
dev.off()

#==fished_n_rel comes from 06_ref_points
png("plots/rel_abn_prop_fish2.png",height=4,width=7,res=400,units='in')
prop_fish + fished_n_rel
dev.off()

median(filter(tot_abn,year==2025)$norm_fish_n)

tot_abn2 <- tot_abn %>%
  dplyr::group_by(stock)%>%
  dplyr::arrange(year)%>%
  dplyr::mutate(norm_abund=abund/max(abund),
          norm_f=fish_m/max(fish_m),
          norm_m=other_m/max(other_m),
         perc_change=(norm_abund-dplyr::lag(norm_abund))/dplyr::lag(norm_abund))%>%
  ungroup()


tot_abn3 <- norm_tot_n %>%
  dplyr::group_by(stock)%>%
  dplyr::arrange(Year)%>%
  dplyr::mutate(perc_change=(tot_n-dplyr::lag(tot_n))/dplyr::lag(tot_n))%>%
  ungroup()


count_big<-tot_abn2%>%
  group_by(stock)%>%
  summarize(ups=sum(perc_change>.40,na.rm=T),
            downs=sum(perc_change< -.40,na.rm=T))%>%
  ungroup()

count_big2<-tot_abn3%>%
  group_by(stock)%>%
  summarize(ups=sum(perc_change>.40,na.rm=T),
            downs=sum(perc_change< -.40,na.rm=T))%>%
  ungroup()

ggplot(tot_abn3,aes(x=Year,y=tot_n))+
  geom_rect(data=filter(tot_abn2,(perc_change)>.30),
            aes(xmin=Year-0.5,xmax=Year+0.5,ymin=-Inf,ymax=Inf),
            fill='green',alpha=0.3,inherit.aes=FALSE)+
  geom_rect(data=filter(tot_abn2,(perc_change)< -.30),
            aes(xmin=Year-0.5,xmax=Year+0.5,ymin=-Inf,ymax=Inf),
            fill='red',alpha=0.3,inherit.aes=FALSE)+
  geom_line()+
  geom_point()+
  facet_wrap(~stock)+
  labs(x="Year",y="Relative abundance")+
  theme_bw()

in_col<-c("#ff5050","#0034c377","#ff505077","#0034c3","#3da550","#ff8f38")

ts_drop<-ggplot(tot_abn2,aes(x=year,y=norm_abund))+
  geom_rect(data=filter(tot_abn2,(perc_change)< -.30),
            aes(xmin=year-0.5,xmax=year+0.5,ymin=-Inf,ymax=Inf),
            fill='grey',alpha=0.8,inherit.aes=FALSE)+
  geom_line(aes(col=stock))+
  geom_point(aes(col=stock))+
  scale_color_manual(values=in_col)+
  facet_wrap(~stock)+
  labs(x="Year",y="Relative abundance")+
  theme_bw()+theme(legend.position='none')+xlab("")

ts_drop2<-ggplot(tot_abn2,aes(x=year,y=norm_abund))+
  geom_rect(data=filter(tot_abn2,(perc_change)< -.30),
            aes(xmin=year-0.5,xmax=year+0.5,ymin=-Inf,ymax=Inf),
            fill='grey',alpha=0.8,inherit.aes=FALSE)+
  geom_line(aes(col=stock))+
  geom_point(aes(col=stock))+
  scale_color_manual(values=in_col)+
  facet_wrap(~stock)+
  labs(x="Year",y="Relative abundance")+
  theme_bw()+theme(legend.position='none',
                   axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank())+xlab("")


drops<-tot_abn2
drops$perc_change[drops$perc_change> -.30]<-0

sub_drop<-drops[,c(1,2,ncol(drops))]
sub_drop$year<-sub_drop$year-1
colnames(sub_drop)[3]<-"perc_change_lag"
drops<-merge(drops,sub_drop,by=c('year','stock'))

png("plots/perc_drops.png",height=4,width=6,res=400,units='in')
collect_drop<-ggplot(drops)+
  geom_bar(aes(x=year,y=perc_change_lag*100,fill=stock),position="stack", stat="identity")+
  theme_bw() + scale_fill_manual(values=in_col)+
  ylab("Percent \nchange")+theme(legend.position='none',axis.line = element_line(colour = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.border = element_blank(),
                                  panel.background = element_blank())+
  xlab("")
print(collect_drop)
dev.off()


drops3<-tot_abn2
sub_drop3<-drops3[,c(1,2,ncol(drops3))]
sub_drop3$year<-sub_drop3$year-1
colnames(sub_drop3)[3]<-"perc_change_lag"
drops3<-merge(drops3,sub_drop3,by=c('year','stock'))
collect_drop_all<-ggplot(drops3)+
  geom_bar(aes(x=year,y=perc_change_lag*100,fill=stock),position="stack", stat="identity")+
  theme_bw() + scale_fill_manual(values=in_col)+
  ylab("Percent \nchange")+theme(legend.position='none',axis.line = element_line(colour = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank())+ylim(-200,0)

png("plots/perc_change.png",height=4,width=6,res=400,units='in')
print(collect_drop_all)
dev.off()


all_drop<-filter(drops3,perc_change_lag<0)%>%
  group_by(year)%>%
  summarize(sum_ch=sum(perc_change_lag))

# colnames(all_drop)[1]<-"Year"
# yoink<-merge(all_drop,temp_mu)
# ggplot(yoink,aes(x=mean_tmp,y=sum_ch))+
#   geom_point()+
#   geom_smooth()+
#   theme_bw()
# modd<-gam(data=filter(yoink,Year>1982),sum_ch~s(mean_tmp))
# summary(modd)
# plot(modd)

ggplot(drops)+
  geom_area(aes(x=year,y=abund,fill=stock))+ 
  scale_fill_manual(values=in_col)+theme_bw()
ggplot(drops)+
  geom_area(aes(x=year,y=norm_abund,fill=stock))+ 
  scale_fill_manual(values=in_col)+theme_bw()

png("plots/big_change_all.png",height=6,width=6,res=400,units='in')
ts_drop/collect_drop+ plot_layout(nrow=2, heights = c(3, 1))
dev.off()

write.csv(filter(drops,perc_change_lag!=0&!is.na(perc_change_lag)),'drops.csv')

data <- drops %>%
  filter(perc_change_lag!=0&!is.na(perc_change_lag)) %>%
  arrange(stock, year) %>%
  group_by(stock) %>%
  mutate(group_id = as.factor(cumsum(year != lag(year, default = first(year)) + 1)) )

 data2 <- data %>%
  dplyr::group_by(stock, group_id) %>%
  dplyr::summarize(start_year = min(year), 
                   end_year = max(year), 
                   avg_perc_chng = mean(perc_change_lag),
                   fishing=mean(norm_f),
                   mortality=mean(norm_m)) %>%
  ungroup()
library(tidyr)
 # Reshape data to long format for plotting
 df_long <- data2 %>%
   pivot_longer(cols = c( mortality, fishing), 
                names_to = "variable", 
                values_to = "value")
 
 df_long$comb<-paste(df_long$stock," (",df_long$start_year,"-",df_long$end_year,")",sep="")
 df_long$yr_range<-paste(" ",df_long$start_year,"-",df_long$end_year,")",sep="")
 df_long$variable[df_long$variable=="fishing"]<-"Fishing\n mortality"
 df_long$variable[df_long$variable=="mortality"]<-"Other\n mortality"
 
 prev_chk<-NULL
 df_long$comb2<-NA
 df_long$comb2[1]<-paste(df_long$stock[1]," (",df_long$start_year[1],"-",df_long$end_year[1],")",sep="")
 for(x in 2:nrow(df_long))
 {
 df_long$comb2[x]<-paste(df_long$stock[x]," (",df_long$start_year[x],"-",df_long$end_year[x],")",sep="")
 if(df_long$stock[x]==df_long$stock[x-1])
   df_long$comb2[x]<-paste("(",df_long$start_year[x],"-",df_long$end_year[x],")",sep="")
 }
 
 # Plot
 reasons<-ggplot(df_long) +
   geom_tile(aes(x=variable,y=comb,fill=value)) +  # Create heatmap tiles
   geom_text(aes(x=variable,y=comb,label = round(value, 2)), color = "black", size = 5) +  # Overlay numbers
   #scale_fill_distiller(palette = "Spectral")+
   scale_fill_gradient2( low = ("blue"),
                         mid = "white",
                         high = ("red"),
                         midpoint = 0.5,)+
   labs(x = "Variable", y = "Stock", fill = "Value") +
   theme_minimal() +  # Clean theme
   theme(axis.text.x = element_text(angle = 45, hjust = 1),
         legend.position='none') +
   xlab("")+
   scale_y_discrete(position='right')
 
 design <- "
  111113
  111113
  111113
  111113
  222223
"
 png("plots/figure2_alt.png",height=7,width=9,res=400,units='in')
 ts_drop2+collect_drop +reasons+ plot_layout(nrow=2, design=design)+ plot_annotation(tag_levels = 'A')
 dev.off()
 
 
 #==try to put 'reasons' on the individual figures using geom_tiles
 #==or just change the color of the rectangles given 

 df_long$col<-1
 df_long$col[which(df_long$variable=="Fishing\n mortality" & df_long$value >0.5)]<-2
 ts_drop3<-ggplot()+
   geom_rect(data=df_long,aes(xmin=start_year +0.5,xmax=end_year +1.5,ymin=-Inf,ymax=Inf,fill=as.factor(col)),
             alpha=0.6,inherit.aes=FALSE)+
   geom_line(data=tot_abn2,aes(x=year,y=norm_abund,col=stock))+
   geom_point(data=tot_abn2,aes(x=year,y=norm_abund,col=stock))+
   scale_color_manual(values=in_col)+
   facet_wrap(~stock)+
   guides(colour = "none")+
   labs(x="Year",y="Relative abundance")+
   theme_bw()+theme(legend.position=c(.91,.07),
                    axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    legend.title=element_blank(),
                    legend.background = element_blank(),
                    legend.box.background = element_blank(),
                    legend.key = element_blank())+xlab("")+
 scale_fill_manual(values=c('grey','red'),
                   labels = c('1' = "F < Average", '2' = "F > Average"))
 
 in_col2<-c('grey','red',"#ff5050","#0034c377","#ff505077","#0034c3","#3da550","#ff8f38")
 
 mod_input<-norm_tot_n
 mod_input$up_m[mod_input$up_m>1.2]<-1.2
 mod_input$dn_m[mod_input$dn_m<0.01]<-0.01
 
 ts_drop3<-ggplot()+
   geom_rect(data=df_long,aes(xmin=start_year +0.5,xmax=end_year +1.5,ymin=-Inf,ymax=Inf,fill=as.factor(col)),
             alpha=0.6,inherit.aes=FALSE)+
   geom_ribbon(data=mod_input,aes(x=Year,ymin=dn_m,ymax=up_m,fill=stock),alpha=.3)+
   geom_line(data=mod_input,aes(x=Year,y=tot_n,col=stock),alpha=.7)+
   geom_point(data=mod_input,aes(x=Year,y=tot_n,col=stock),alpha=.7,size=.51)+
   scale_color_manual(values=in_col)+
   facet_wrap(~stock)+
   guides(colour = "none")+
   labs(x="Year",y="Relative abundance")+
   theme_bw()+theme(legend.position='none',
                    axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    legend.title=element_blank(),
                    legend.background = element_blank(),
                    legend.box.background = element_blank(),
                    legend.key = element_blank(),
                    plot.margin=unit(c(.1,.1,.1,.1),'cm'))+xlab("")+
    scale_fill_manual(values=in_col2)
 filter(mod_input,stock=='PIRKC')
 
 design <- "
  1
  1
  1
  1
  1
  2
"
 png("plots/figure2alt.png",height=5,width=8,res=400,units='in')
 ts_drop3+collect_drop + plot_layout(nrow=2, design=design)+ plot_annotation(tag_levels = 'A')
 dev.off()
 