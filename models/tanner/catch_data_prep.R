ret_dat<-read.csv("C:\\Users\\cody.szuwalski\\Work\\ebs_crab\\data\\tanner catch\\RetainedCatchZCs.CrabFisheries.Scaled.Binned.csv")
tot_dat<-read.csv("C:\\Users\\cody.szuwalski\\Work\\ebs_crab\\data\\tanner catch\\TotalCatchZCs.CrabFisheries.Scaled.Binned.csv")

ret_dat<-read.csv("C:\\Users\\cody.szuwalski\\Work\\ebs_crab\\data\\tanner catch\\RetainedCatch.ZCs.AllAsTCF.ScaledAndBinned.csv")
tot_dat<-read.csv("C:\\Users\\cody.szuwalski\\Work\\ebs_crab\\data\\tanner catch\\TotalCatch.ZCs.CrabFisheries.Scaled.Binned.csv")

library(dplyr)
library(ggplot2)
library(reshape2)

#==retained data
ret_m<-filter(ret_dat,sex=='male')%>%
  group_by(year,size)%>%
  summarize(tot_n=sum(abundance))

ret_m_tot<-ret_m%>%
  group_by(year)%>%
  summarize(tot_n=sum(tot_n))

out_ret_sc<-dcast(ret_m,year~size,value.var='tot_n')
yrm<-out_ret_sc[,-1]
ret_size_comp<-sweep(yrm,1,ret_m_tot$tot_n,FUN="/")
rownames(ret_size_comp)<-out_ret_sc[,1]

#==total data
tot_m<-filter(tot_dat,sex=='male'&area=='all EBS')%>%
  group_by(year,size)%>%
  summarize(tot_n=sum(abundance))

tot_m_tot<-tot_m%>%
  group_by(year)%>%
  summarize(tot_n=sum(tot_n))

out_tot_sc<-dcast(tot_m,year~size,value.var='tot_n')
yrm<-out_tot_sc[,-1]
tot_size_comp<-sweep(yrm,1,tot_m_tot$tot_n,FUN="/")
rownames(tot_size_comp)<-out_tot_sc[,1]

write.csv(ret_m_tot,'data/tanner catch/1_retained_tanner1.csv')
write.csv(ret_size_comp,'data/tanner catch/1_retained_tanner_sc1.csv')
write.csv(tot_m_tot,'data/tanner catch/1_total_tanner1.csv')
write.csv(tot_size_comp,'data/tanner catch/1_total_tanner_sc1.csv')

new<-ggplot()+
  geom_point(data=tot_m_tot,aes(x=year,y=tot_n),col='blue')+
  geom_point(data=ret_m_tot,aes(x=year,y=tot_n),col='red')+
  geom_line(data=tot_m_tot,aes(x=year,y=tot_n),col='blue')+
  geom_line(data=ret_m_tot,aes(x=year,y=tot_n),col='red')+
  theme_bw()


old<-ggplot()+
  geom_point(data=tot_m_tot,aes(x=year,y=tot_n),col='blue')+
  geom_point(data=ret_m_tot,aes(x=year,y=tot_n),col='red')+
  geom_line(data=tot_m_tot,aes(x=year,y=tot_n),col='blue')+
  geom_line(data=ret_m_tot,aes(x=year,y=tot_n),col='red')+
  theme_bw()
