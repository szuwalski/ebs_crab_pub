#==script for automatically updating the .DAT file for EBS snow crab
library(crabpack)
library(mgcv)
library(dplyr)
library(reshape2)
library(ggplot2)

#==define size ranges of interest
bin_size<-5
sizes<-seq(47.5,167.5,bin_size)
styr<-1975
endyr<-2025
years<-seq(styr,endyr)

channel <- "API"
in_species<-"RKC"
in_region<-"EBS"
in_district<-"BB"
specimen_data <- crabpack::get_specimen_data(species = in_species,
                                             region = in_region,
                                             years = c(styr:endyr),
                                             district=in_district,
                                             channel = channel)

#==get indices (do this for comparison, but ultimately do not need this...except for CVs)
indices <- crabpack::calc_bioabund(crab_data = specimen_data,
                                       species = in_species,
                                       region = in_region,
                                       district=in_district,
                                       years =c(styr:endyr),
                                       sex = "male",
                                       size_min=min(sizes)-(bin_size/2),
                                       shell_condition = "all_categories")

#==get n at carapace width by shell condition
n_at_size <- crabpack::calc_bioabund(crab_data = specimen_data,
                                   species = in_species,
                                   region = in_region,
                                   years = c(styr:endyr),
                                   district=in_district,
                                   sex = "male",
                                   size_min=min(sizes)-(bin_size/2),
                                   shell_condition = "all_categories",
                                   bin_1mm = TRUE)

calc_ind <- crabpack::calc_bioabund(crab_data = specimen_data,
                                     species = in_species,
                                     region = in_region,
                                     district=in_district,
                                     years = c(styr:endyr),
                                     sex = "male",
                                     size_min=min(sizes)-(bin_size/2),
                                     bin_1mm = FALSE)

# compute breaks from midpoints (half-step on each side)
breaks <- c(min(sizes) - 2.5, sizes + 2.5)  # e.g., 25–135 range

# assign bins and aggregate
n_at_size_in <- n_at_size %>%
  mutate(size_bin = cut(SIZE_1MM, breaks = breaks, labels = sizes, include.lowest = TRUE)) %>%
  group_by(YEAR, size_bin) %>%
  summarise(total_abundance = sum(ABUNDANCE, na.rm = TRUE), .groups = "drop")

#==split new shell crab at size into mature/immature
n_at_size_dat<-dcast(n_at_size_in,YEAR~size_bin)
n_at_size_dat_in<-n_at_size_dat[,-1]
rownames(n_at_size_dat_in)<-n_at_size_dat[,1]

#==get index of mature and immature animals
#==add zeros
yrs<-n_at_size_dat[,1]
allyr<-seq(min(yrs),max(yrs))
use_n<-matrix(0,ncol=length(sizes),nrow=length(allyr))
rownames(use_n)<-allyr
use_n_sc<-use_n

for(x in 1:length(allyr))
{
  use_n[x,]  <- unlist(n_at_size_dat_in[(match(allyr[x],as.numeric(rownames(n_at_size_dat_in)))),-c(ncol(n_at_size_dat_in))])
  if(!is.na(sum(use_n[x,])))
   use_n_sc[x,]<-use_n[x,]/sum(use_n[x,])
}
use_n[is.na(use_n)]<-0
use_n_sc[is.na(use_n_sc)]<-0
surv_ind<-apply(use_n,1,sum)

#==put code for calculating size transition matrix in here 


##################################################################
# write .DAT file
# .pin file below
#===============================================================

# Specify output file
outfile <- "models/bbrkc/test/rkc.dat"
# Create an empty file (overwrite if exists)
file.create(outfile)

# Write header and values
cat("# bristol bay red king crab pop dy model\n", file = outfile)
cat(paste("# Generated on: ", Sys.Date()),"\n", file = outfile, append = TRUE)

# inputs
cat("# start year", file = outfile, append = TRUE,"\n")
cat(styr, file = outfile, append = TRUE,"\n")
cat("# end year", file = outfile, append = TRUE,"\n")
cat(endyr, file = outfile, append = TRUE,"\n")

cat("# number of years of survey data", file = outfile, append = TRUE,"\n")
cat(nrow(use_n), file = outfile, append = TRUE,"\n")

cat("# years of survey data", file = outfile, append = TRUE,"\n")
cat(rownames(use_n), file = outfile, append = TRUE,"\n")

cat("# number of sizes", file = outfile, append = TRUE,"\n")
cat(ncol(use_n), file = outfile, append = TRUE,"\n")

cat("# sizes (midpoints)", file = outfile, append = TRUE,"\n")
cat(sizes, file = outfile, append = TRUE,"\n")

cat("# survey numbers", file = outfile, append = TRUE,"\n")
cat(surv_ind, file = outfile, append = TRUE,"\n")

cat("# survey numbers at size", file = outfile, append = TRUE,"\n")
write.table(
  round(use_n_sc,4),
  file = outfile,
  append = TRUE,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cat("# probability of molting", file = outfile, append = TRUE,"\n")
input_p_molt<-1 - (1/(1+exp(-0.148*(sizes-138.79))))
cat(input_p_molt, file = outfile, append = TRUE,"\n")

in_stm<-read.csv("models/bbrkc/test/size_trans_bbrkc.csv",header=F)
cat("# size transition matrix", file = outfile, append = TRUE,"\n")
write.table(
  in_stm,
  file = outfile,
  append = TRUE,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cat("# CV survey numbers", file = outfile, append = TRUE,"\n")
tmp<-calc_ind%>%
  select(YEAR,ABUNDANCE_CV)

yrs<-tmp$YEAR
allyr<-seq(min(yrs),max(yrs))
surv_cv<-rep(0,length=length(allyr))
names(surv_cv)<-allyr

for(x in 1:length(allyr))
  surv_cv[x]  <- unlist(tmp[(match(allyr[x],yrs)),2])
surv_cv[is.na(surv_cv)]<-0

cat(round(surv_cv,2), file = outfile, append = TRUE,"\n")

#===weights and model options below
cat("# surv_eff_samp", file = outfile, append = TRUE,"\n")
cat(25, file = outfile, append = TRUE,"\n")
cat("# log_mu_m", file = outfile, append = TRUE,"\n")
cat(c(-1.7), file = outfile, append = TRUE,"\n")

cat("# est_m_devs", file = outfile, append = TRUE,"\n")
cat(1, file = outfile, append = TRUE,"\n")
cat("# est_q_devs", file = outfile, append = TRUE,"\n")
cat(-1, file = outfile, append = TRUE,"\n")
cat("# est_sigma_m", file = outfile, append = TRUE,"\n")
cat(-1, file = outfile, append = TRUE,"\n")

cat("# sigma_m_mu (prior on avg, notn the devs, devs are in the .pin file", file = outfile, append = TRUE,"\n")
cat(c(0.01), file = outfile, append = TRUE,"\n")

cat("# smooth_q_weight", file = outfile, append = TRUE,"\n")
cat(0.001, file = outfile, append = TRUE,"\n")
cat("# smooth_m_weight", file = outfile, append = TRUE,"\n")
cat(c(0.01), file = outfile, append = TRUE,"\n")
cat("# est_log_m_mu", file = outfile, append = TRUE,"\n")
cat(-1, file = outfile, append = TRUE,"\n")
cat("# est_sigma_q", file = outfile, append = TRUE,"\n")
cat(-1, file = outfile, append = TRUE,"\n")

cat("# smooth_f_weight", file = outfile, append = TRUE,"\n")
cat(0.01, file = outfile, append = TRUE,"\n")

cat("# surve_sel_cv", file = outfile, append = TRUE,"\n")
cat(0.05, file = outfile, append = TRUE,"\n")
cat("# smooth_surv_weight", file = outfile, append = TRUE,"\n")
cat(0, file = outfile, append = TRUE,"\n")
cat("# est_molt_prob", file = outfile, append = TRUE,"\n")
cat(1, file = outfile, append = TRUE,"\n")


cat("\n# end of file\n", file = outfile, append = TRUE)

######################################################
# .PIN file creation
#++++++++++++++++++++++++++++++++++++++++++++++++++

# Specify output file
outfile <-  "models/bbrkc/test/rkc.pin"
file.create(outfile)

# Write header and values
cat("# bristol bay red king crab pop dy model .pin file\n", file = outfile)
cat(paste("# Generated on: ", Sys.Date()),"\n", file = outfile, append = TRUE)

cat("# log_n_init", file = outfile, append = TRUE,"\n")
cat(log(rep(mean(use_n[1,]),length(sizes))), file = outfile, append = TRUE,"\n")
cat("# nat_m_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")
cat("# nat_q_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")

cat("# log_avg_rec", file = outfile, append = TRUE,"\n")
cat(20, file = outfile, append = TRUE,"\n")
cat("# rec_devs", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")

cat("# sigma_m", file = outfile, append = TRUE,"\n")
cat(c(1), file = outfile, append = TRUE,"\n")
cat("# log_m_mu", file = outfile, append = TRUE,"\n")
cat(c(-1.7), file = outfile, append = TRUE,"\n")
cat("# prop_rec", file = outfile, append = TRUE,"\n")
cat(c(14.4,2.9), file = outfile, append = TRUE,"\n")
cat("# sigma_q", file = outfile, append = TRUE,"\n")
cat(c(.1,.1), file = outfile, append = TRUE,"\n")

cat("# log_f", file = outfile, append = TRUE,"\n")
cat(-1, file = outfile, append = TRUE,"\n")
cat("# f_dev", file = outfile, append = TRUE,"\n")
cat(rep(0,length(years)), file = outfile, append = TRUE,"\n")
cat("# fish_ret_sel_50", file = outfile, append = TRUE,"\n")
cat(136, file = outfile, append = TRUE,"\n")
cat("# fish_ret_sel_slope", file = outfile, append = TRUE,"\n")
cat(0.54, file = outfile, append = TRUE,"\n")
cat("# fish_tot_sel_50", file = outfile, append = TRUE,"\n")
cat(114, file = outfile, append = TRUE,"\n")
cat("# fish_tot_sel_slope", file = outfile, append = TRUE,"\n")
cat(0.1, file = outfile, append = TRUE,"\n")

cat("# surv_sel_50", file = outfile, append = TRUE,"\n")
cat(78, file = outfile, append = TRUE,"\n")
cat("# surv_sel_slope", file = outfile, append = TRUE,"\n")
cat(0.1, file = outfile, append = TRUE,"\n")


cat("# molt_sel_50", file = outfile, append = TRUE,"\n")
cat(139, file = outfile, append = TRUE,"\n")
cat("# molt_sel_slope", file = outfile, append = TRUE,"\n")
cat(0.55, file = outfile, append = TRUE,"\n")



#####################################################
# catch stuff now
#####################################################

cat_dat<-read.csv("models/bbrkc/test/total_catch.csv")
tot_cat<-filter(cat_dat,group!="female"&target_stock=="BBRKC")%>%
  group_by(crab_year)%>%
  summarize(tot_n=sum(total_catch_n))
write.csv(tot_cat,"models/bbrkc/test/obs_tot_cat_in.csv")
