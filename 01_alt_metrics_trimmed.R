# Trimmed duplicate of 01_alt_metrics.R for manuscript_publish
# Keeps only code needed for manuscript-included figures and required upstream data objects/files.


library(maps)
library("rnaturalearth")
library(interp)
library(RColorBrewer)
library(reshape2) # for melt
library(mgcv)  
library(PBSmapping)
library(mapdata)    #some additional hires data
#library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(ggridges)
library(crabpack)
library(mgcv)
library(dplyr)
library(reshape2)
library(ggplot2)


#===pull data
endyr<-2025
channel <- "API"
in_region<-"EBS"

#==opilio
opilio <- crabpack::get_specimen_data(species = "SNOW",
                                             region = in_region,
                                             years = c(1982:endyr),
                                             channel = channel)

opie_cpue <- crabpack::calc_cpue(crab_data = opilio,
                            species = "SNOW",
                            region = "EBS",
                            district = "ALL",
                            sex="male",
                            size_min=25,
                            size_max=135,
                            bin_1mm=TRUE)

use_opie<-merge(opie_cpue,select(opilio$haul,"YEAR","STATION_ID","GEAR_TEMPERATURE"))

#==tanner

tanner <- crabpack::get_specimen_data(species = "TANNER",
                                      region = in_region,
                                      years = c(1975:endyr),
                                      channel = channel)

tan_cpue <- crabpack::calc_cpue(crab_data = tanner,
                                 species = "TANNER",
                                 region = "EBS",
                                 district = "ALL",
                                 sex="male",
                                 size_min=25,
                                 size_max=185,
                                 bin_1mm=TRUE)

use_tan<-merge(tan_cpue,select(tanner$haul,"YEAR","STATION_ID","GEAR_TEMPERATURE"))

#==red king crab
rkc <- crabpack::get_specimen_data(species = "RKC",
                                      region = in_region,
                                      years = c(1975:endyr),
                                      channel = channel)

rkc_cpue <- crabpack::calc_cpue(crab_data = rkc,
                                species = "RKC",
                                region = "EBS",
                                district = "ALL",
                                sex="male",
                                size_min=45,
                                size_max=220,
                                bin_1mm=TRUE)

use_rkc<-merge(rkc_cpue,select(rkc$haul,"YEAR","STATION_ID","GEAR_TEMPERATURE"))


#==blue king crab
bkc <- crabpack::get_specimen_data(species = "BKC",
                                     region = in_region,
                                     years = c(1975:endyr),
                                     channel = channel)

bkc_cpue <- crabpack::calc_cpue(crab_data = bkc,
                                species = "BKC",
                                region = "EBS",
                                district = "ALL",
                                sex="male",
                                size_min=45,
                                size_max=200,
                                bin_1mm=TRUE)

use_bkc<-merge(bkc_cpue,select(bkc$haul,"YEAR","STATION_ID","GEAR_TEMPERATURE"))

#============================================
#==calculate temperature occupied + mean size
#============================================
#==temperature occupied opilio
use_df_op<-filter(use_opie,nchar(STATION_ID)<5)
op_occ_temp_yr<-use_df_op%>%
  group_by(YEAR)%>%
  summarize(Temperature=weighted.mean(GEAR_TEMPERATURE, w = CPUE, na.rm = T),
            avg_size=weighted.mean(SIZE_1MM,w=CPUE,na.rm=T))
colnames(op_occ_temp_yr)<-c("Year","Temperature","Size")
op_occ_temp_yr$stock<-"Snow"

#==temperature occupied tanner
use_df_tn<-filter(use_tan,nchar(STATION_ID)<5)
tn_occ_temp_yr<-use_df_tn%>%
  group_by(YEAR)%>%
  summarize(Temperature=weighted.mean(GEAR_TEMPERATURE, w = CPUE, na.rm = T),
            avg_size=weighted.mean(SIZE_1MM,w=CPUE,na.rm=T))
colnames(tn_occ_temp_yr)<-c("Year","Temperature","Size")
tn_occ_temp_yr$stock<-"Tanner"

#==W of 168 for prib dist, S of 58.39 lat
#==bbrkc
use_df_bbrkc<-filter(use_rkc,nchar(STATION_ID)<5& LONGITUDE >-168)
bbrkc_occ_temp_yr<-use_df_bbrkc%>%
  group_by(YEAR)%>%
  summarize(Temperature=weighted.mean(GEAR_TEMPERATURE, w = CPUE, na.rm = T),
            avg_size=weighted.mean(SIZE_1MM,w=CPUE,na.rm=T))
colnames(bbrkc_occ_temp_yr)<-c("Year","Temperature","Size")
bbrkc_occ_temp_yr$stock<-"BBRKC"

#==pirck
use_df_bbrkc<-filter(use_rkc,nchar(STATION_ID)<5& LONGITUDE <=-168)
pirkc_occ_temp_yr<-use_df_bbrkc%>%
  group_by(YEAR)%>%
  summarize(Temperature=weighted.mean(GEAR_TEMPERATURE, w = CPUE, na.rm = T),
            avg_size=weighted.mean(SIZE_1MM,w=CPUE,na.rm=T))
colnames(pirkc_occ_temp_yr)<-c("Year","Temperature","Size")
pirkc_occ_temp_yr$stock<-"PIRKC"

#==smbkc
use_df_smbkc<-filter(use_bkc,nchar(STATION_ID)<5& LATITUDE >58.39)
smbkc_occ_temp_yr<-use_df_smbkc%>%
  group_by(YEAR)%>%
  summarize(Temperature=weighted.mean(GEAR_TEMPERATURE, w = CPUE, na.rm = T),
            avg_size=weighted.mean(SIZE_1MM,w=CPUE,na.rm=T))
colnames(smbkc_occ_temp_yr)<-c("Year","Temperature","Size")
smbkc_occ_temp_yr$stock<-"SMBKC"


#==pibkc
use_df_smbkc<-filter(use_bkc,nchar(STATION_ID)<5& LATITUDE <58.39)
pibkc_occ_temp_yr<-use_df_smbkc%>%
  group_by(YEAR)%>%
  summarize(Temperature=weighted.mean(GEAR_TEMPERATURE, w = CPUE, na.rm = T),
            avg_size=weighted.mean(SIZE_1MM,w=CPUE,na.rm=T))
colnames(pibkc_occ_temp_yr)<-c("Year","Temperature","Size")
pibkc_occ_temp_yr$stock<-"PIBKC"

plotters<-rbind(op_occ_temp_yr,tn_occ_temp_yr,bbrkc_occ_temp_yr,
                pirkc_occ_temp_yr,smbkc_occ_temp_yr,pibkc_occ_temp_yr)


png("plots/fig_4.png",height=6,width=8,res=400,units='in')
p<-ggplot(plotters)+
  geom_line(aes(x=Year,y=Temperature,col=stock),lwd=2)+
  scale_color_manual(values=c("#ff5050","#0034c377","#ff505077","#0034c3","#3da550","#ff8f38"))+
  theme_bw()+ylab("Temperature occupied (C)")

q<-ggplot(plotters)+
  geom_line(aes(x=Year,y=Size,col=stock),lwd=2)+
  scale_color_manual(values=c("#ff5050","#0034c377","#ff505077","#0034c3","#3da550","#ff8f38"))+
  theme_bw()+ylab("Average size (mm)")

p/q+plot_layout(guides='collect')
dev.off()


write.csv(plotters,"data/alt_metrics_calc.csv")
