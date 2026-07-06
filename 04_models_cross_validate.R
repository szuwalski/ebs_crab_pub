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
png("plots/model_concurvity.png",height=8,width=8,res=350,units='in') 
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

mort_plot_alt_trans<-ggplot()+
  geom_point(data=filter(unc_mort,Year<2023&Year!=2020),aes(x=Year,y=1-exp(-est_m)),col='darkgrey')+
  geom_errorbar(data=filter(unc_mort,Year<2023&Year!=2020),aes(x=Year,ymin=1-exp(-dn_m),ymax=1-exp(-up_m)),col='grey')+
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

# mortality <- density + temperature + competitor + size
library(patchwork)
design<-"112
         112
         112"
png("plots/all_mort_all3.png",height=9,width=6,res=350,units='in') 
mort_plot_alt_trans +mort_plot_1 + plot_layout(nrow=2, design=design)
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
