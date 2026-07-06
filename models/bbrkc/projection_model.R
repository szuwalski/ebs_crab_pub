
#==SET UP AN OPTION TO STOP THE FITTING AT A GIVEN YEAR
#==AND THEN PROJECT FROM THAT YEAR AND COMPARE TO THE ACTUAL MODEL ESTIMATES
#==can do this for all the red king crab species...publish the mortality paper first

#==BIG FUN
library(PBSmodelling)
library(patchwork)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggridges)
library(rsoi)
library(zoo)
library(mgcv)
library(DHARMa)
library(itsadug)
library(mgcViz)
library(directlabels)

#==read in parameters/model output/covariates
outs<-readList("models/bbrkc/test/rkc.rep")
#outs<-readList("models/smbkc/test/bkc.rep")
all_output<-read.csv("data/all_output.csv")
alt_metrics_calc<-read.csv("data/alt_metrics_calc.csv")
unc_mort<-read.csv("data/uncertainty_mort.csv")
unique(all_output$process)

#==fit relationships between process and covariates
use_species<-"BBRKC"
m_dat<-filter(unc_mort,stock==use_species)
alt_d<-filter(alt_metrics_calc,stock==use_species)
dens<-filter(all_output,species==use_species&process=="Abundance")
colnames(dens)[2]<-"Abund"
gam_dat<-merge(m_dat,alt_d,by='Year')
gam_dat<-merge(gam_dat,dens,by='Year')

gam_dat <- gam_dat %>%
  dplyr::select(Year, est_m, Temperature, Size, Abund, sd) %>%
  filter(
    is.finite(Year),
    is.finite(est_m),
    is.finite(Size),
    is.finite(Abund),
    is.finite(sd),
    sd > 0
  ) %>%
  mutate(m_resp = 1 - exp(-est_m))

m_mod<-gam(data=gam_dat,1-exp(-est_m)~s(Temperature,k=4)+s(Size,k=4)+s(Abund,k=4),family = betar(link = "logit"),weights=1/sd)
summary(m_mod)
plot(m_mod,pages=1)

chonk<-melt(gam_dat,"Year")
ggplot(chonk)+
  geom_line(aes(x=Year,y=value))+
  facet_wrap(~variable,scales='free')+
  theme_bw()

#==write projection function
detach(outs)
attach(outs)
proj_popdy<-function(mod_m=FALSE,
                     proj_yr=5,
                     proj_f=.5,
                     proj_temp=3,
                     discard_survival=0.85)
{
  hist_yrs<-seq(styr,endyr)
  all_years<-seq(styr,endyr+proj_yr)  
  size_n<-length(sizes)
  n_size_pred<-matrix(ncol=size_n,nrow=length(all_years))
  n_size_pred[1:length(hist_yrs),] <-sweep(matrix(ncol=size_n,as.numeric(`pred numbers at size`)),1,numbers_pred,"*")
  pred_retained_size<-matrix(ncol=size_n,nrow=length(all_years))
  gurp<-matrix(ncol=size_n,as.numeric(`pred_retained_size_comp`))
  gurp[gurp=="NaN"]<-0
  pred_retained_size[1:(length(hist_yrs)-1),] <-sweep(gurp,1,pred_retained_n,"*")
  retain_fish_sel = ret_fish_sel

  nat_m<-matrix(ncol=size_n,nrow=length(all_years))
  selectivity<-matrix(ncol=size_n,nrow=length(all_years))
  f_mort<-rep(proj_f,length(all_years))
  
  #==input historical par estimates
for(x in 1:length(hist_yrs))
  {
    nat_m[x,] <- `natural mortality`[x,]
    selectivity[x,] = `survey selectivity`[x,]
    f_mort[x]<- est_fishing_mort[x]
  }

for(y in (length(hist_yrs)):(length(hist_yrs)+proj_yr-1))
{
  
  nat_m[y,]<-median(`natural mortality`[,1])
  if(mod_m==TRUE)
  {
  #==implement model for M
  #==calculate density and average size
  tmp_abund<-sum(n_size_pred[y,])/sum(n_size_pred[1,])
  tmp_size<-weighted.mean(sizes,weights=n_size_pred[y,])
  tmp_temp<-proj_temp
  mod_pred<-predict(m_mod,newdata=list(Temperature=tmp_temp,Size=tmp_size,Abund=tmp_abund),type='response')
  nat_m[y,]<- -log(1-mod_pred)
  }
  
  temp_n <- n_size_pred[y,] * exp(-1*(0.17)*exp( nat_m[y,]))
  
  # fishery
  
    temp_catch_n = temp_n * (1 -exp(-(f_mort[y]*total_fish_sel)))
    pred_retained_size[y,] = temp_catch_n*retain_fish_sel 
    #pred_tot_size_comp[y,] = temp_catch_n
    temp_n = temp_n * exp(-((f_mort[y]*total_fish_sel)))
    temp_n = temp_n + temp_catch_n*(1-retain_fish_sel*(discard_survival))	

  # growth
  trans_n = size_trans %*% (temp_n*in_prob_molt)
  temp_n = trans_n + (temp_n*(1-in_prob_molt))
  
  # recruitment
  temp_n[1] <- temp_n[1] + median(recruits)*temp_prop_rec[1]
  temp_n[2] <- temp_n[2] + median(recruits)*temp_prop_rec[2]
  temp_n[3] <- temp_n[3] + median(recruits)*temp_prop_rec[3]
  
  # natural mortality		
  n_size_pred[y+1,] = temp_n * exp(-1*(0.83)*exp(nat_m[y,]))
  
}

list(n_size_pred=n_size_pred,
     f_mort=f_mort,
     nat_m=nat_m,
     retained_catch=pred_retained_size,
     years=all_years)

}


no_mod_warm<-proj_popdy(mod_m=FALSE,
           proj_yr=5,
           proj_f=median(est_fishing_mort),
           proj_temp=3,
           discard_survival=0.85)

no_mod_cold<-proj_popdy(mod_m=FALSE,
                        proj_yr=5,
                        proj_f=median(est_fishing_mort),
                        proj_temp=1,
                        discard_survival=0.85)

no_mod_hot<-proj_popdy(mod_m=FALSE,
                        proj_yr=5,
                        proj_f=median(est_fishing_mort),
                        proj_temp=5,
                        discard_survival=0.85)

mod_warm<-proj_popdy(mod_m=TRUE,
                        proj_yr=5,
                        proj_f=median(est_fishing_mort),
                        proj_temp=3,
                        discard_survival=0.85)

mod_cold<-proj_popdy(mod_m=TRUE,
                        proj_yr=5,
                        proj_f=median(est_fishing_mort),
                        proj_temp=1,
                        discard_survival=0.85)

mod_hot<-proj_popdy(mod_m=TRUE,
                       proj_yr=5,
                       proj_f=median(est_fishing_mort),
                       proj_temp=5,
                       discard_survival=0.85)


#==pull the numbers at size, plot three columns for temp, two diff polygons for 
proj_m<-data.frame(year=c(mod_hot$year,mod_warm$year,mod_cold$year,no_mod_hot$year,no_mod_warm$year,no_mod_cold$year),
           mort=c(mod_hot$nat_m[,1],mod_warm$nat_m[,1],mod_cold$nat_m[,1],no_mod_hot$nat_m[,1],no_mod_warm$nat_m[,1],no_mod_cold$nat_m[,1]),
           temp=rep(c(rep('hot',length(mod_hot$year)),rep('warm',length(mod_hot$year)),rep('cold',length(mod_hot$year))),2),
           model=c(rep('GAM',3*length(mod_hot$year)),rep('median',3*length(mod_hot$year))))

ggplot(proj_m)+
  geom_line(aes(x=year,y=mort,col=temp))+
  facet_wrap(~model)+
  theme_bw()

#===run the projection with and without environmental stuff, overplot the projections. is there a difference?
#===plot the 





library(dplyr)
library(mgcv)
library(ggplot2)
library(tidyr)
library(purrr)

# ---------------------------
# user settings
# ---------------------------
stock_name  <- "BBRKC"
n_peel      <- 10   # number of times to peel off the most recent year
k_smooth    <- 4

# ---------------------------
# build data
# ---------------------------
m_dat  <- filter(unc_mort, stock == stock_name)
alt_d  <- filter(alt_metrics_calc, stock == stock_name)
dens   <- filter(all_output, species == stock_name & process == "Abundance")

colnames(dens)[2] <- "Abund"

gam_dat <- merge(m_dat, alt_d, by = "Year")
gam_dat <- merge(gam_dat, dens, by = "Year")

# keep only needed columns and drop bad rows
gam_dat <- gam_dat %>%
  dplyr::select(Year, est_m, Temperature, Size, Abund, sd) %>%
  filter(
    is.finite(Year),
    is.finite(est_m),
    is.finite(Temperature),
    is.finite(Size),
    is.finite(Abund),
    is.finite(sd),
    sd > 0
  ) %>%
  mutate(resp = 1 - exp(-est_m))

# betar requires response strictly inside (0,1)
eps <- 1e-6
gam_dat <- gam_dat %>%
  mutate(resp = pmin(pmax(resp, eps), 1 - eps))

# adjust n_peel if too large
all_years <- sort(unique(gam_dat$Year))
n_peel <- min(n_peel, length(all_years) - 5)  # leave at least a few years to fit

# ---------------------------
# helper to fit one peeled model
# ---------------------------
fit_one_peel <- function(dat, peel_i, k_smooth = 4, n_grid = 200) {
  yrs <- sort(unique(dat$Year))
  max_year_included <- yrs[length(yrs) - peel_i]
  
  dsub <- dat %>% filter(Year <= max_year_included)
  
  mod <- gam(
    resp ~ s(Temperature, k = k_smooth) +
      s(Size, k = k_smooth) +
      s(Abund, k = k_smooth),
    data = dsub,
    family = betar(link = "logit"),
    weights = 1 / sd,
    method = "REML"
  )
  
  # common reference values for "other" covariates
  temp_med  <- median(dsub$Temperature, na.rm = TRUE)
  size_med  <- median(dsub$Size, na.rm = TRUE)
  abund_med <- median(dsub$Abund, na.rm = TRUE)
  
  # prediction grids
  grid_temp <- data.frame(
    Temperature = seq(min(dat$Temperature, na.rm = TRUE), max(dat$Temperature, na.rm = TRUE), length.out = n_grid),
    Size = size_med,
    Abund = abund_med
  )
  
  grid_size <- data.frame(
    Temperature = temp_med,
    Size = seq(min(dat$Size, na.rm = TRUE), max(dat$Size, na.rm = TRUE), length.out = n_grid),
    Abund = abund_med
  )
  
  grid_abund <- data.frame(
    Temperature = temp_med,
    Size = size_med,
    Abund = seq(min(dat$Abund, na.rm = TRUE), max(dat$Abund, na.rm = TRUE), length.out = n_grid)
  )
  
  # term-specific predictions on link scale
  p_temp <- predict(mod, newdata = grid_temp, type = "terms", se.fit = TRUE)
  p_size <- predict(mod, newdata = grid_size, type = "terms", se.fit = TRUE)
  p_abund <- predict(mod, newdata = grid_abund, type = "terms", se.fit = TRUE)
  
  out_temp <- data.frame(
    covariate = "Temperature",
    x = grid_temp$Temperature,
    fit = p_temp$fit[, "s(Temperature)"],
    se = p_temp$se.fit[, "s(Temperature)"],
    max_year_included = max_year_included
  )
  
  out_size <- data.frame(
    covariate = "Size",
    x = grid_size$Size,
    fit = p_size$fit[, "s(Size)"],
    se = p_size$se.fit[, "s(Size)"],
    max_year_included = max_year_included
  )
  
  out_abund <- data.frame(
    covariate = "Abund",
    x = grid_abund$Abund,
    fit = p_abund$fit[, "s(Abund)"],
    se = p_abund$se.fit[, "s(Abund)"],
    max_year_included = max_year_included
  )
  
  bind_rows(out_temp, out_size, out_abund)
}

# ---------------------------
# fit all peeled models
# peel_i = 0 means full data
# peel_i = 1 means remove most recent year
# etc.
# ---------------------------
smooth_dat <- map_dfr(0:n_peel, ~fit_one_peel(gam_dat, peel_i = .x, k_smooth = k_smooth))

# make year a factor so colors are discrete and ordered
smooth_dat <- smooth_dat %>%
  mutate(max_year_included = factor(max_year_included, levels = sort(unique(max_year_included))))

# ---------------------------
# plot
# ---------------------------
ggplot(smooth_dat, aes(x = x, y = fit, color = max_year_included, group = max_year_included)) +
  geom_line(linewidth = 1) +
  facet_wrap(~covariate, scales = "free_x", ncol = 3) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Estimated smooth (partial effect, link scale)",
    color = "Data through year",
    title = paste0(stock_name, ": GAM smooths under sequential retrospective peeling")
  )