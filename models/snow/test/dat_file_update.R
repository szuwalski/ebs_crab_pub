#==script for automatically updating the .DAT file for EBS snow crab
library(crabpack)
library(mgcv)
library(dplyr)
library(reshape2)
library(ggplot2)

#==define size ranges of interest
sizes<-seq(27.5,132.5,5)
styr<-1982
endyr<-2025
years<-seq(styr,endyr)

channel <- "API"
specimen_data <- crabpack::get_specimen_data(species = "SNOW",
                                             region = "EBS",
                                             years = c(styr:endyr),
                                             channel = channel)

#==get indices (do this for comparison, but ultimately do not need this...except for CVs)
indices <- crabpack::calc_bioabund(crab_data = specimen_data,
                                       species = "SNOW",
                                       region = "EBS",
                                       years =c(styr:endyr),
                                       sex = "male",
                                       size_min=25,
                                       shell_condition = "all_categories")

#==get n at carapace width by shell condition
n_at_size <- crabpack::calc_bioabund(crab_data = specimen_data,
                                   species = "SNOW",
                                   region = "EBS",
                                   years = c(styr:endyr),
                                   sex = "male",
                                   size_min=25,
                                   shell_condition = "all_categories",
                                   bin_1mm = TRUE)

calc_ind <- crabpack::calc_bioabund(crab_data = specimen_data,
                                     species = "SNOW",
                                     region = "EBS",
                                     years = c(styr:endyr),
                                     sex = "male",
                                     size_min=30,
                                     bin_1mm = FALSE,
                                     crab_category='all_categories')

# compute breaks from midpoints (half-step on each side)
breaks <- c(min(sizes) - 2.5, sizes + 2.5)  # e.g., 25–135 range

# assign bins and aggregate
n_at_size_new <- filter(n_at_size,SHELL_TEXT%in%c("new_hardshell","soft_molting")) %>%
  mutate(size_bin = cut(SIZE_1MM, breaks = breaks, labels = sizes, include.lowest = TRUE)) %>%
  group_by(YEAR, size_bin) %>%
  summarise(total_abundance = sum(ABUNDANCE, na.rm = TRUE), .groups = "drop")

n_at_size_old <- filter(n_at_size,SHELL_TEXT%in%c("oldshell","soft_movery_oldshelllting")) %>%
  mutate(size_bin = cut(SIZE_1MM, breaks = breaks, labels = sizes, include.lowest = TRUE)) %>%
  group_by(YEAR, size_bin) %>%
  summarise(total_abundance = sum(ABUNDANCE, na.rm = TRUE), .groups = "drop")

#==get probability of having undergone terminal molt
male_maturity_data <- crabpack::get_male_maturity(species = "SNOW",
                                                  region = "EBS",
                                                  district = "ALL",
                                                  channel = channel)

#==make a data.frame that fills in maturity data for missing years
#==and interpolates using GAMs to the size bins needed
out_p_mat<-matrix(ncol=length(sizes),nrow=length(years))

for(x in 1:length(years))
{
 tmp<-filter(male_maturity_data$male_mat_ratio,YEAR==years[x])
 if(nrow(tmp)>0) 
 {
  mod<-gam(data=tmp,PROP_MATURE~s(SIZE_BIN,k=15)) 
  out_p_mat[x,]<-predict(mod,newdata=data.frame(SIZE_BIN=sizes))
 }
}
no_dat<-which(is.na(out_p_mat[,1]))
for(y in no_dat)
 out_p_mat[y,]<-apply(out_p_mat,2,median,na.rm=T)
out_p_mat[out_p_mat<0]<-0
out_p_mat[out_p_mat>1]<-1

rownames(out_p_mat)<-years

#==split new shell crab at size into mature/immature
n_at_size_new_wide<-dcast(n_at_size_new,YEAR~size_bin)
n_at_size_imm<-n_at_size_new_wide[,-1]
rownames(n_at_size_imm)<-n_at_size_new_wide[,1]
n_at_size_new_mat<-n_at_size_new_wide[,-1]
rownames(n_at_size_new_mat)<-n_at_size_new_wide[,1]
for(x in 1:nrow(n_at_size_imm))
{
  n_at_size_imm[x,]<-n_at_size_imm[x,]*(1-out_p_mat[which(rownames(out_p_mat)==rownames(n_at_size_imm)[x]),])
  n_at_size_new_mat[x,]<-n_at_size_new_mat[x,]*(out_p_mat[which(rownames(out_p_mat)==rownames(n_at_size_new_mat)[x]),])  
}

n_at_size_mat<-dcast(n_at_size_old,YEAR~size_bin)[,-1]+n_at_size_new_mat

#==get index of mature and immature animals
use_imm_n<-n_at_size_imm[,-c(1,ncol(n_at_size_imm))]
use_mat_n<-n_at_size_mat[,-c(1,ncol(n_at_size_mat))]

imm_ind<-apply(use_imm_n,1,sum)
mat_ind<-apply(use_mat_n,1,sum)


#==need size transition matrix in here now too
# insert code here
# write.csv(size_trans,'models/ebs/data_files/size_trans.csv')
size_trans<-read.csv("models/ebs/data_files/size_trans.csv",header=F)



##################################################################
# write .DAT file
# .pin file below
#===============================================================

# Specify output file
outfile <- "models/ebs/data_files/snow_down.dat"
# Create an empty file (overwrite if exists)
file.create(outfile)

# Write header and values
cat("# eastern Bering Sea snow crab pop dy model\n", file = outfile)
cat(paste("# Generated on: ", Sys.Date()),"\n", file = outfile, append = TRUE)

# inputs
cat("# start year", file = outfile, append = TRUE,"\n")
cat(styr, file = outfile, append = TRUE,"\n")
cat("# end year", file = outfile, append = TRUE,"\n")
cat(endyr, file = outfile, append = TRUE,"\n")

cat("# number of years of survey data", file = outfile, append = TRUE,"\n")
cat(nrow(use_mat_n), file = outfile, append = TRUE,"\n")

cat("# years of survey data", file = outfile, append = TRUE,"\n")
cat(n_at_size_new_wide[,1], file = outfile, append = TRUE,"\n")

cat("# number of sizes", file = outfile, append = TRUE,"\n")
cat(ncol(use_mat_n), file = outfile, append = TRUE,"\n")

cat("# sizes (midpoints)", file = outfile, append = TRUE,"\n")
cat(colnames(use_mat_n), file = outfile, append = TRUE,"\n")

cat("# immature numbers", file = outfile, append = TRUE,"\n")
cat(imm_ind, file = outfile, append = TRUE,"\n")

cat("# immature numbers at size", file = outfile, append = TRUE,"\n")
write.table(
  round(use_imm_n),
  file = outfile,
  append = TRUE,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cat("# mature numbers", file = outfile, append = TRUE,"\n")
cat(mat_ind, file = outfile, append = TRUE,"\n")

cat("# mature numbers at size", file = outfile, append = TRUE,"\n")
write.table(
  round(use_mat_n),
  file = outfile,
  append = TRUE,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cat("# proportion undegoing terminal molt at size", file = outfile, append = TRUE,"\n")
#==need to lag by a year to make dynamics work right in the model
#==final year needs to be the average, which happens to be the first year
in_p_mat<-out_p_mat[-1,-1]
avg_p<-matrix(out_p_mat[1,-1],nrow=1)
in_p_mat<-rbind(in_p_mat,avg_p)

write.table(
  round(in_p_mat,3),
  file = outfile,
  append = TRUE,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


cat("# size transition matrix", file = outfile, append = TRUE,"\n")
cat(mat_ind, file = outfile, append = TRUE,"\n")

cat("# mature numbers at size", file = outfile, append = TRUE,"\n")
write.table(
  size_trans,
  file = outfile,
  append = TRUE,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

emp_sel<-unlist(read.csv("models/ebs/data_files/empirical_selectivity.csv",header=F))
cat("# empirical selectivity", file = outfile, append = TRUE,"\n")
cat(round(emp_sel,3), file = outfile, append = TRUE,"\n")


cat("# CV immature numbers", file = outfile, append = TRUE,"\n")
cat(round(filter(calc_ind,CATEGORY=='small_male')$ABUNDANCE_CV,3)
, file = outfile, append = TRUE,"\n")

cat("# CV mature numbers", file = outfile, append = TRUE,"\n")
cat(round(filter(calc_ind,CATEGORY=='legal_male')$ABUNDANCE_CV,3)
    , file = outfile, append = TRUE,"\n")

#===weights and model options below
cat("# mat_eff_samp", file = outfile, append = TRUE,"\n")
cat(25, file = outfile, append = TRUE,"\n")
cat("# imm_eff_samp", file = outfile, append = TRUE,"\n")
cat(25, file = outfile, append = TRUE,"\n")
cat("# log_mu_m", file = outfile, append = TRUE,"\n")
cat(c(-1.2,-1.2,-1.2), file = outfile, append = TRUE,"\n")

cat("# est_m_devs", file = outfile, append = TRUE,"\n")
cat(1, file = outfile, append = TRUE,"\n")
cat("# est_q_devs", file = outfile, append = TRUE,"\n")
cat(-1, file = outfile, append = TRUE,"\n")
cat("# est_m_mat_devs", file = outfile, append = TRUE,"\n")
cat(1, file = outfile, append = TRUE,"\n")
cat("# est_q_mat_devs", file = outfile, append = TRUE,"\n")
cat(-1, file = outfile, append = TRUE,"\n")
cat("# est_sigma_m", file = outfile, append = TRUE,"\n")
cat(-1, file = outfile, append = TRUE,"\n")

cat("# sigma_m_mu (prior on avg, notn the devs, devs are in the .pin file", file = outfile, append = TRUE,"\n")
cat(c(0.01,0.008,0.01), file = outfile, append = TRUE,"\n")

cat("# smooth_q_weight", file = outfile, append = TRUE,"\n")
cat(0.001, file = outfile, append = TRUE,"\n")
cat("# smooth_m_weight", file = outfile, append = TRUE,"\n")
cat(c(0.01,0.05), file = outfile, append = TRUE,"\n")
cat("# est_log_m_mu", file = outfile, append = TRUE,"\n")
cat(1, file = outfile, append = TRUE,"\n")
cat("# est_sigma_q", file = outfile, append = TRUE,"\n")
cat(-1, file = outfile, append = TRUE,"\n")
cat("# est_m_lg_devs", file = outfile, append = TRUE,"\n")
cat(-1, file = outfile, append = TRUE,"\n")
cat("# smooth_f_weight", file = outfile, append = TRUE,"\n")
cat(0.01, file = outfile, append = TRUE,"\n")

cat("# surve_sel_cv", file = outfile, append = TRUE,"\n")
cat(0.05, file = outfile, append = TRUE,"\n")
cat("# surve_sel_cv2", file = outfile, append = TRUE,"\n")
cat(0.005, file = outfile, append = TRUE,"\n")
cat("# smooth_surv_weight", file = outfile, append = TRUE,"\n")
cat(0, file = outfile, append = TRUE,"\n")
cat("# est_sel", file = outfile, append = TRUE,"\n")
cat(1, file = outfile, append = TRUE,"\n")
cat("# large_cutoff", file = outfile, append = TRUE,"\n")
cat(100, file = outfile, append = TRUE,"\n")
cat("# est_tv_fish_sel", file = outfile, append = TRUE,"\n")
cat(0, file = outfile, append = TRUE,"\n")

cat("\n# end of file\n", file = outfile, append = TRUE)

######################################################
# .PIN file creation
#++++++++++++++++++++++++++++++++++++++++++++++++++

# Specify output file
outfile <- "models/ebs/data_files/snow_down.pin"
file.create(outfile)

# Write header and values
cat("# eastern Bering Sea snow crab pop dy model .pin file\n", file = outfile)
cat(paste("# Generated on: ", Sys.Date()),"\n", file = outfile, append = TRUE)

cat("# log_n_imm", file = outfile, append = TRUE,"\n")
cat(0, file = outfile, append = TRUE,"\n")
cat("# log_n_mat", file = outfile, append = TRUE,"\n")
cat(0, file = outfile, append = TRUE,"\n")

cat("# nat_m_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")
cat("# nat_m_mat_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")
cat("# nat_m_lg_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")
cat("# nat_q_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")
cat("# nat_q_mat_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")
cat("# nat_q_lg_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")

cat("# log_avg_rec", file = outfile, append = TRUE,"\n")
cat(0, file = outfile, append = TRUE,"\n")
cat("# rec_devs", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")

cat("# sigma_m", file = outfile, append = TRUE,"\n")
cat(c(1,1,1), file = outfile, append = TRUE,"\n")
cat("# log_m_mu", file = outfile, append = TRUE,"\n")
cat(c(-1.35,-0.77,-1.1), file = outfile, append = TRUE,"\n")
cat("# prop_rec", file = outfile, append = TRUE,"\n")
cat(c(14.4,2.9), file = outfile, append = TRUE,"\n")
cat("# sigma_q", file = outfile, append = TRUE,"\n")
cat(c(.1,.1), file = outfile, append = TRUE,"\n")

cat("# log_f", file = outfile, append = TRUE,"\n")
cat(0.13, file = outfile, append = TRUE,"\n")
cat("# f_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")
cat("# fish_ret_sel_50", file = outfile, append = TRUE,"\n")
cat(97, file = outfile, append = TRUE,"\n")
cat("# fish_ret_sel_50_post", file = outfile, append = TRUE,"\n")
cat(102, file = outfile, append = TRUE,"\n")
cat("# fish_ret_sel_slope", file = outfile, append = TRUE,"\n")
cat(0.21, file = outfile, append = TRUE,"\n")
cat("# fish_ret_sel_slope_post", file = outfile, append = TRUE,"\n")
cat(0.29, file = outfile, append = TRUE,"\n")
cat("# fish_tot_sel_offset", file = outfile, append = TRUE,"\n")
cat(-9, file = outfile, append = TRUE,"\n")
cat("# fisH_tot_sel_offset_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")

cat("# fish_tot_sel_slope", file = outfile, append = TRUE,"\n")
cat(-1.6, file = outfile, append = TRUE,"\n")
cat("# surv_omega", file = outfile, append = TRUE,"\n")
cat(0.56, file = outfile, append = TRUE,"\n")
cat("# surv_alpha1", file = outfile, append = TRUE,"\n")
cat(0.1, file = outfile, append = TRUE,"\n")
cat("# surv_beta1", file = outfile, append = TRUE,"\n")
cat(110, file = outfile, append = TRUE,"\n")
cat("# surv_alpha2", file = outfile, append = TRUE,"\n")
cat(0.18, file = outfile, append = TRUE,"\n")
cat("# surv_beta2", file = outfile, append = TRUE,"\n")
cat(0.41, file = outfile, append = TRUE,"\n")