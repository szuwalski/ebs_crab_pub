#==read in par and input
library(reshape2)
library(ggplot2)
library(mgcv)  
library(dplyr)
library(ggridges)
library(png)
library(PBSmodelling)
library(patchwork)


#==QUESTION
#==Can I improve performance when mortality events are expected with different shaped control rules?
#==or some other algorithmic approach?

#==need an example that links mortality to density

#==run 02_figure_2 and 03_figure_3 first to get keep_uncertainty_fishn
norm_fish_n<-keep_uncertainty_fishn
unq_st<-unique(norm_fish_n$stock)
for(x in 1:length(unq_st))
{
  indi<-which(norm_fish_n$stock==unq_st[x])  
  norm_fish_n[indi,2:4]<- norm_fish_n[indi,2:4]/max(norm_fish_n[indi,2])
}

#==add PIBKC, PIRKC
rep_files<-c("/models/bbrkc/test/rkc.rep")
maturity<-c(120)
outs_in<-list(list())

for(x in 1:length(rep_files))
  outs_in[[x]]<-readList(paste(getwd(),rep_files[x],sep=""))

names(outs_in[[1]])
in_outs<-outs_in[[1]]
pop_dy<-function(in_outs,
                 f_mort,
                 proj_yr,
                 constant_m=TRUE,
                 sample_m=FALSE,
                 disc_surv=0.80,
                 steepness,
                 srr_r0,
                 srr_s0,
                 n0,
                 dens_m=FALSE,
                 min_m=0.1,
                 max_m=0.5)
{
  n_at_size<-matrix(ncol=length(in_outs$sizes),nrow=proj_yr)
  n_at_size[1,]<-in_outs$"numbers_pred"[1]*in_outs$"pred numbers at size"[1,]
  MMB<-rep(0,proj_yr)
  recruits<-rep(0,proj_yr)
  ret_catch<-matrix(ncol=length(in_outs$sizes),nrow=proj_yr)
  
  #==fill nat M
  if(sample_m==FALSE)
    nat_m<-matrix(median(in_outs$`natural mortality`),ncol=length(in_outs$sizes),nrow=proj_yr) 
  
  if(sample_m==TRUE)
  {
    nat_m<-matrix(ncol=length(in_outs$sizes),nrow=proj_yr)
    for(x in 1:proj_yr)
     nat_m[x,]<-sample(in_outs$`natural mortality`,size=1)
  }
 
  if(dens_m==TRUE)
    nat_m<-matrix(median(in_outs$`natural mortality`),ncol=length(in_outs$sizes),nrow=proj_yr)   
  
  for(x in 1:(proj_yr-1))
  {
   if(dens_m==TRUE&x>10) # make M a function of density 
   {
    slope<-(max_m-min_m)/n0
    nat_m[x,]<-  min_m + slope*sum(n_at_size[x,])
   }

   tmp_n <- n_at_size[x,]*exp(-0.17*nat_m[x,]) 
 
   # fishery
   tmp_catch <- tmp_n * (1-exp(-f_mort*in_outs$total_fish_sel))
   ret_catch[x,] <- tmp_catch*in_outs$ret_fish_sel
   tmp_n <- tmp_n *exp(-f_mort*in_outs$total_fish_sel)
   tmp_n <- tmp_n + (tmp_catch*(1-in_outs$ret_fish_sel))*disc_surv
   
   # growth
   trans_n <- c(in_outs$size_trans %*% (tmp_n*in_outs$prob_molt))
   tmp_n <- trans_n + tmp_n*(1-in_outs$prob_molt)
   
   # recruits
   MMB[x]<-sum(tmp_n[in_outs$sizes>120])
   num   <- 4*srr_r0*steepness*MMB[x]
   denom <- (1-steepness)*srr_r0*(srr_s0/srr_r0)+(5*steepness-1)*MMB[x]
   recruits[x]<- num/denom 
   
   tmp_n[1]<-0.3*recruits[x]
   tmp_n[2]<-0.5*recruits[x]
   tmp_n[3]<-0.2*recruits[x]
   
   n_at_size[x+1,] <- tmp_n*exp(-0.83*nat_m[x,]) 
  }
  list(MMB=MMB,recruits=recruits,yield=apply(ret_catch,1,sum),nat_m=nat_m,abund=apply(n_at_size,1,sum))
  }


proj_yr<-200
no_fish<-pop_dy(in_outs=outs_in[[1]],
           f_mort=0,
           proj_yr=proj_yr,
           sample_m=FALSE,
           disc_surv=0.80,
           steepness=1,
           srr_r0=100000,
           srr_s0=100)

b0<-no_fish$MMB[proj_yr-1]
n0<-no_fish$abund[proj_yr-1]

#==loop over f_mort
f_morts<-seq(0,1,0.01)
pop_out<-list(list())
yields1<-rep(0,length(f_morts))
mmbs1<-rep(0,length(f_morts))
for(y in 1:length(f_morts))
{
  pop_out[[y]]<-pop_dy(in_outs=outs_in[[1]],
           f_mort=f_morts[y],
           proj_yr=proj_yr,
           sample_m=FALSE,
           disc_surv=0.80,
           steepness=0.5,
           srr_r0=100000,
           srr_s0=b0,
           n0=n0)
  yields1[y]<-pop_out[[y]]$yield[proj_yr-1]
  mmbs1[y]<-pop_out[[y]]$MMB[proj_yr-1]
}

plot(yields1~f_morts)
plot(yields1~mmbs1,type='b')
plot(mmbs1~f_morts)
eq_df<-data.frame(f_mort=f_morts,yield=yields1,MMB=mmbs1)
eq_df$f_mort[which(eq_df$yield==max(eq_df$yield))]

#==density dependence in M
#==NEED TO CALCULATE N0 AND B0 BASED ON DENSITY DEPENDENT MORTALITY
pop_out<-list(list())
yields2<-rep(0,length(f_morts))
mmbs2<-rep(0,length(f_morts))
out_m<-rep(0,length(f_morts))
for(y in 1:length(f_morts))
{
  pop_out[[y]]<-pop_dy(in_outs=outs_in[[1]],
                       f_mort=f_morts[y],
                       proj_yr=proj_yr,
                       sample_m=FALSE,
                       disc_surv=0.80,
                       steepness=0.5,
                       srr_r0=100000,
                       srr_s0=b0,
                       n0=n0,
                       dens_m=TRUE)
  yields2[y]<-pop_out[[y]]$yield[proj_yr-1]
  mmbs2[y]<-pop_out[[y]]$MMB[proj_yr-1]
  out_m[y]<-pop_out[[y]]$nat_m[proj_yr-1]
}

plot(out_m~f_morts)
plot(yields2~f_morts)
plot(yields1~mmbs1,type='b')
plot(mmbs1~f_morts)
eq_df<-data.frame(f_mort=f_morts,yield=yields1,MMB=mmbs1)
eq_df$f_mort[which(eq_df$yield==max(eq_df$yield))]




#==stochastic now
#==loop over f_mort
nsims<-1000
yields<-matrix(ncol=length(f_morts),nrow=nsims)
mmbs<-matrix(ncol=length(f_morts),nrow=nsims)
for(x in 1:nsims)
for(y in 1:length(f_morts))
{
  pop_out<-pop_dy(in_outs=outs_in[[1]],
                       f_mort=f_morts[y],
                       proj_yr=proj_yr,
                       sample_m=TRUE,
                       disc_surv=0.80,
                       steepness=0.5,
                       srr_r0=100000,
                       srr_s0=b0,
                       n0=n0)
  yields[x,y]<-pop_out$yield[proj_yr-1]
  mmbs[x,y]<-pop_out$MMB[proj_yr-1]
}

colnames(yields)<-f_morts
rownames(yields)<-seq(1,nsims)
melted<-melt(yields)
colnames(melted)<-c("Sim","f_mort","yield")
ggplot()+
  geom_boxplot(data=melted,aes(x=f_mort,y=yield,group=f_mort),fill='grey')+
  geom_line(data=eq_df,aes(x=f_mort,y=yield),col='red',lwd=1.25)+
  theme_bw()

ggplot()+
  geom_point(data=melted,aes(x=f_mort,y=yield,group=f_mort),fill='grey')+
  geom_line(data=eq_df,aes(x=f_mort,y=yield),col='red',lwd=1.25)+
  theme_bw()

colnames(mmbs)<-f_morts
rownames(mmbs)<-seq(1,nsims)
melted_mmb<-melt(mmbs)
colnames(melted_mmb)<-c("Sim","f_mort","MMB")

ggplot()+
  geom_boxplot(data=melted_mmb,aes(x=f_mort,y=MMB,group=f_mort),fill='grey')+
  geom_line(data=eq_df,aes(x=f_mort,y=MMB),col='red',lwd=1.25)+
  theme_bw()

merged<-merge(melted,melted_mmb)
ggplot()+
  geom_point(data=filter(merged,yield>0.01),aes(x=MMB,y=yield))+
  geom_line(data=eq_df,aes(x=MMB,y=yield),col='red',lwd=1.25)+
  theme_bw()

ggplot()+
  geom_point(data=merged,aes(x=f_mort,y=yield,col=MMB))+
  geom_line(data=eq_df,aes(x=f_mort,y=yield),col='red',lwd=1.25)+
  theme_bw()

# png("plots/yield_curv.png",height=5,width=7,res=350,units='in') 
# ggplot()+
#   geom_boxplot(data=merged,aes(x=f_mort,y=yield,group=f_mort),fill='grey')+
#   geom_line(data=eq_df,aes(x=f_mort,y=yield),col='red',lwd=1.25)+
#   theme_bw()+xlab("Fishing mortality")+
#   ylab("Yield")+xlim(0,1.5)
# dev.off()

ggplot()+
  geom_boxplot(data=merged,aes(x=f_mort,y=MMB,group=f_mort),fill='grey')+
  geom_line(data=eq_df,aes(x=f_mort,y=MMB),col='red',lwd=1.25)+
  theme_bw()+xlab("Fishing mortality")+
  ylab("MMB")+xlim(0,1.5)


save_b0<-rep(0,nsims)
for(x in 1:nsims)
  save_b0[x]<-pop_dy(in_outs=outs_in[[1]],
           f_mort=0,
           proj_yr=proj_yr,
           sample_m=TRUE,
           disc_surv=0.80,
           steepness=1,
           srr_r0=100000,
           srr_s0=100,n0=n0)$MMB[proj_yr-1]
other_b0<-median(save_b0)

yields2<-matrix(ncol=length(f_morts),nrow=nsims)
mmbs2<-matrix(ncol=length(f_morts),nrow=nsims)
for(x in 1:nsims)
  for(y in 1:length(f_morts))
  {
    pop_out<-pop_dy(in_outs=outs_in[[1]],
                    f_mort=f_morts[y],
                    proj_yr=proj_yr,
                    sample_m=TRUE,
                    disc_surv=0.80,
                    steepness=0.5,
                    srr_r0=100000,
                    srr_s0=other_b0,n0=n0)
    yields2[x,y]<-pop_out$yield[proj_yr-1]
    mmbs2[x,y]<-pop_out$MMB[proj_yr-1]
  }

colnames(yields2)<-f_morts
rownames(yields2)<-seq(1,nsims)
melted_y2<-melt(yields2)
colnames(melted_y2)<-c("Sim","f_mort","yield")

png("plots/yield_curv2.png",height=5,width=7,res=350,units='in') 
f_comp<-ggplot()+
  geom_boxplot(data=melted_y2,aes(x=f_mort,y=yield,group=f_mort),fill='grey')+
  geom_line(data=eq_df,aes(x=f_mort,y=yield),col='red',lwd=1.25)+
  theme_bw()+xlab("Fishing mortality")+
  ylab("Yield")+xlim(0,0.6)
dev.off()

eq_df$f_mort[which(eq_df$yield==max(eq_df$yield))]


meds<-melted_y2%>%
  dplyr::group_by(f_mort)%>%
  dplyr::summarize(med_yld=median(yield))
library(mgcv)
print(meds,n=60)
mod<-gam(data=meds,med_yld~s(f_mort))
summary(mod)
plot(mod)
meds$preds<-predict(mod)
meds$f_mort[which(meds$med_yld==max(meds$med_yld))]

f_comp<-ggplot()+
  geom_boxplot(data=melted_y2,aes(x=f_mort,y=yield,group=f_mort),fill='grey')+
  geom_line(data=eq_df,aes(x=f_mort,y=yield),col='red',lwd=1.25)+
  geom_line(data=meds,aes(x=f_mort,y=preds),col='black',lwd=1.25)+
  theme_bw()+xlab("Fishing mortality")+
  ylab("Yield")+xlim(0,0.6)

colnames(mmbs2)<-f_morts
rownames(mmbs2)<-seq(1,nsims)
melted_mmb2<-melt(mmbs2)
colnames(melted_mmb2)<-c("Sim","f_mort","MMB")

bio_comp<-ggplot()+
  geom_boxplot(data=melted_mmb2,aes(x=f_mort,y=MMB,group=f_mort),fill='grey')+
  geom_line(data=eq_df,aes(x=f_mort,y=MMB),col='red',lwd=1.25)+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+xlim(0,0.6)


mod_input<-norm_fish_n
mod_input$up_m[mod_input$up_m>1.2]<-1.2
mod_input$dn_m[mod_input$dn_m<0.01]<-0.01
mod_input$fish<-"Regularly fished"
mod_input$fish[which(mod_input$stock%in%c("PIBKC","PIRKC","SMBKC"))]<-"Rarely fished"
in_col<-c("#ff5050","#0034c377","#ff505077","#0034c3","#3da550","#ff8f38")

fished_n_rel<-ggplot()+
  geom_ribbon(data=mod_input,aes(x=Year,ymin=dn_m,ymax=up_m,fill=stock),alpha=.3)+
  geom_line(data=mod_input,aes(x=Year,y=tot_n,col=stock),alpha=.7)+
  geom_point(data=mod_input,aes(x=Year,y=tot_n,col=stock),alpha=.7,size=.51)+
  scale_color_manual(values=in_col)+
  facet_wrap(~fish,ncol=1)+
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
  scale_fill_manual(values=in_col)





png("plots/yield_curv3.png",height=8,width=4,res=350,units='in') 
bio_comp / f_comp
dev.off()

png("plots/yield_curv3_stock.png",height=7,width=9,res=350,units='in') 
bio_comp / f_comp | fishable_n 
dev.off()

#==need 05_size_changes ran first
png("plots/yield_curv6_stock.png",height=7,width=11,res=350,units='in') 
prop_fish2 | fished_n_rel |bio_comp / f_comp 
dev.off()

merged2<-merge(melted_y2,melted_mmb2)
ggplot()+
  geom_point(data=filter(merged2,yield>0.01),aes(x=MMB,y=yield))+
  geom_line(data=eq_df,aes(x=MMB,y=yield),col='red',lwd=1.25)+
  theme_bw()

ggplot()+
  geom_point(data=merged2,aes(x=f_mort,y=yield,col=MMB))+
  geom_line(data=eq_df,aes(x=f_mort,y=yield),col='red',lwd=1.25)+
  theme_bw()
