#=======================
# RKC in RTMB
#=======================
library(RTMB)
#devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
#library(TMBhelper)

data<-NULL
# read in survey data
input<-"models/pirkc/rkc.DAT"
data$styr <- scan(input,skip=2,n=1,quiet=T)
data$endyr <- scan(input,skip=4,n=1,quiet=T)
data$year_n <- scan(input,skip=6,n=1,quiet=T)
data$years <- scan(input,skip=8,n=data$year_n,quiet=T)
data$size_n <- scan(input,skip=10,n=1,quiet=T)
data$sizes <- scan(input,skip=12,n=data$size_n,quiet=T)
data$n_obs<-scan(input,skip=14,n=data$year_n,quiet=T)
data$n_at_size_obs <- matrix(scan(input,skip=16,n=data$year_n*data$size_n,quiet=T),ncol=data$size_n,byrow=T)
data$molt_prob_in<-scan(input,skip=65,n=data$size_n,quiet=T)
data$size_trans <- (matrix(scan(input,skip=67,n=data$size_n*data$size_n,quiet=T),ncol=data$size_n,byrow=T))
data$surv_cv <- scan(input,skip=104,n=data$year_n,quiet=T)

# read in catch data
cat_input<-"models/pirkc/catch_dat.DAT"
data$cat_year_n <- scan(cat_input,skip=1,n=1,quiet=T)
data$cat_years <- scan(cat_input,skip=3,n=data$cat_year_n,quiet=T)
data$ret_cat <- scan(cat_input,skip=5,n=data$cat_year_n,quiet=T)
data$ret_cat_sc_year_n <- scan(cat_input,skip=7,n=1,quiet=T)
data$discard_survival <- 0.8
data$surv_samp_n<-30
data$ret_samp_n<-50
data$tot_samp_n<-50
data$sm_f_wt<-0
data$sm_r_wt<-0
data$sm_m_wt<-5

#=read in parameters
par_file<-"models/pirkc/rkc.par"
params<-NULL
params$log_n_init<-scan(par_file,skip=2,n=data$size_n,quiet=T)
params$nat_m_dev<-rep(0,data$year_n)
params$log_avg_rec<-scan(par_file,skip=8,n=1,quiet=T)
params$rec_dev<-rep(0,data$year_n)
params$log_m_mu<-scan(par_file,skip=14,n=1,quiet=T)
params$prop_rec_in<-scan(par_file,skip=16,n=2,quiet=T)
params$log_f<-scan(par_file,skip=20,n=1,quiet=T)
params$f_dev<-scan(par_file,skip=22,n=data$cat_year_n,quiet=T)
params$fish_ret_sel_50<-scan(par_file,skip=24,n=1,quiet=T)
params$fish_ret_sel_slope<-scan(par_file,skip=26,n=1,quiet=T)
params$fish_tot_sel_50<-scan(par_file,skip=28,n=1,quiet=T)
params$fish_tot_sel_slope<-scan(par_file,skip=30,n=1,quiet=T)
params$surv_sel_50<-scan(par_file,skip=32,n=1,quiet=T)
params$surv_sel_slope<-scan(par_file,skip=34,n=1,quiet=T)
params$molt_50<-scan(par_file,skip=36,n=1,quiet=T)
params$molt_slope<-scan(par_file,skip=38,n=1,quiet=T)
params$sigma_f<-.3
params$sigma_r<-.5
params$sigma_m<-.5


if (TRUE) {
  attach(data)       ## discouraged but
  attach(params) ## VERY handy while developing!
}

pop_dy<-function(log_n_init,nat_m_dev,log_avg_rec,rec_dev,sigma_r,
                 sigma_m,log_m_mu,prop_rec_in,
                 log_f,f_dev,sigma_f,fish_ret_sel_50,fish_ret_sel_slope,
                 fish_tot_sel_50,fish_tot_sel_slope,
                 surv_sel_50,surv_sel_slope,molt_50,molt_slope)
{
  tot_fish_sel <- 1/(1+exp(-fish_tot_sel_slope*(sizes-fish_tot_sel_50)))
  ret_fish_sel <- 1/(1+exp(-fish_ret_sel_slope*(sizes-fish_ret_sel_50)))
  surv_sel <- 1/(1+exp(-surv_sel_slope*(sizes-surv_sel_50)))
  molt_prob <- 1-(1/(1+exp(-molt_slope*(sizes-molt_50))))
  
  prop_rec<-c(20,prop_rec_in)
  tot_prop_rec<-sum(prop_rec)
  prop_rec<-prop_rec/tot_prop_rec
  
  f_mort<-rep(0,year_n)
  fish_ind <- which(!is.na(match(years,cat_years)))
  f_mort[fish_ind] <- exp(log_f+f_dev)
  tot_rec <- exp(log_avg_rec+rec_dev)
  
  in_m_dev<-nat_m_dev-mean(nat_m_dev)
  
  nat_m<-exp(log_m_mu+in_m_dev)

  #==storage
  n_pred<-matrix(0,ncol=size_n,nrow=year_n)
  n_pred[1,]<-exp(log_n_init)
  pred_tot_n_at_l<-matrix(0,ncol=size_n,nrow=year_n)
  pred_ret_n_at_l<-matrix(0,ncol=size_n,nrow=year_n)
  
  for(x in 1:(year_n-1))
  {
    temp_n <- n_pred[x,]*exp(-0.17*(nat_m[x])) 
    
    # fishery
    pred_tot_n_at_l[x,] <- temp_n * (1-exp(-f_mort[x]*tot_fish_sel))  
    pred_ret_n_at_l[x,] <- pred_tot_n_at_l[x,] * ret_fish_sel
    
    # take crab out
    temp_n <- temp_n*exp(-f_mort[x]*tot_fish_sel) 
    
    # toss some back, some survive
    temp_n <- temp_n + ((1-ret_fish_sel)*pred_tot_n_at_l[x,])*discard_survival
    
    # some crab molt
    trans_n <- size_trans%*%(temp_n*molt_prob)
    
    # others don't
    temp_n <-trans_n + (temp_n*(1-molt_prob))
    
    # recruitment
    
    for(y in 1:3)
      temp_n[y] <- temp_n[y] + tot_rec[x]*prop_rec[y]
    
    # remaining mortality
    n_pred[x+1,]<-temp_n*exp(-0.83*(nat_m[x])) 
  }
  
  #==negative log likelihood
  #===survey index
  num_like <- 0
  temp_n_pred<-n_pred
  for(p in 1:nrow(n_pred))
    temp_n_pred[p,]<-n_pred[p,]*surv_sel
  
  tot_n_pred<-apply(temp_n_pred,1,sum,na.rm=T)
  
  for(z in 1:year_n)
  {
    term <-((log(tot_n_pred[z])-log(n_obs[z]+0.01))^2)/(2*sqrt(log(1+surv_cv[z]^2)))
    num_like <- num_like + term*(n_obs[z]!=0)
  }
  
  #===retained catch 
  ret_like <- 0
  ret_n_pred<-apply(pred_ret_n_at_l,1,sum)
  fish_ind <-which(!is.na(match(years,cat_years)))
  for(z in 1:length(fish_ind))
    ret_like <- ret_like + ((log(ret_n_pred[fish_ind[z]])-log(ret_cat[z]))^2)/(2*sqrt(log(1+0.05^2)))
  
  #===total catch 
  tot_like <- 0
  # tot_c_pred<-apply(pred_tot_n_at_l,1,sum)
  # tot_ind <-which(!is.na(match(years,dcat_years)))
  # for(z in 1:length(tot_ind))
  #   tot_like <- tot_like + ((log(tot_c_pred[tot_ind[z]])-log(tot_cat[z]))^2)/(2*sqrt(log(1+0.1^2)))
  
  #===survey size composition
  num_sc_like<-0
  for(z in 1:year_n)
    for(w in 1:size_n)
    {
      term <- surv_samp_n*(n_at_size_obs[z,w]+0.0001) * log((surv_sel[w]*n_pred[z,w]/tot_n_pred[z]+0.0001)/(n_at_size_obs[z,w]+0.0001))
      num_sc_like <- num_sc_like + term*(n_obs[z]!=0)*(n_at_size_obs[z,w]>0.001)
    }
  
  num_sc_like <- -1*num_sc_like
  
  #===retained size composition
  ret_sc_like<-0
  # ret_sc_ind <-which(!is.na(match(years,ret_cat_sc_years)))
  # for(z in 1:(length(ret_sc_ind)-1))
  # {
  #   for(w in 1:size_n)
  #   {
  #     term <- ret_samp_n*(ret_at_size_obs[z,w]+0.0001) * log(0.0001+(pred_ret_n_at_l[ret_sc_ind[z],w]/ret_n_pred[ret_sc_ind[z]])/(ret_at_size_obs[z,w]+0.0001))
  #     ret_sc_like <- ret_sc_like +term*(ret_at_size_obs[z,w]>0.001)
  #   }
  # }
  # ret_sc_like <- -1*ret_sc_like
  
  #===total size composition
  tot_sc_like<-0
  # tot_sc_ind <-which(!is.na(match(years,tot_cat_sc_years)))
  # for(z in 1:(length(tot_sc_ind)-1))
  # {
  #   for(w in 1:size_n)
  #   {
  #     term <-  tot_samp_n*(tot_at_size_obs[z,w]+0.0001) * log(0.0001+(pred_tot_n_at_l[tot_sc_ind[z],w]/tot_c_pred[tot_sc_ind[z]])/(tot_at_size_obs[z,w]+0.0001))
  #     tot_sc_like <-tot_sc_like + term*(tot_at_size_obs[z,w]>0.001 )
  #   }
  # }
  # tot_sc_like <- -1*tot_sc_like  
  
  #==random effects
  nat_m_like<-0
  nat_m_like<- -sum(dnorm(nat_m_dev,  0, sigma_m, log = TRUE))
  
  sig_m_like<-0
  #sig_m_like<-sum(dnorm(sigma_m,0.5,0.01))
  
  f_like<-0
  f_like   <- -sum(dnorm(f_dev,  0, sigma_f, log = TRUE))
  
  rec_like<-0
  rec_like <- -sum(dnorm(rec_dev,0, sigma_r, log = TRUE))
  
  nat_m_prior<-0
  #nat_m_prior<-((0.21-exp(log_m_mu))^2)/(2*0.001^2)
  
  f_prior<-0
  #f_prior<-((0.5-exp(log_f))^2)/(2*1^2)
  
  smooth_f<-0
  smooth_r<-0
  smooth_m<-0
  # smooth_f<-sm_f_wt*sum(diff(f_dev,difference=2))
  # smooth_r<-sm_r_wt*sum(diff(rec_dev,difference=2))
  # smooth_m<- -sm_m_wt*sum(diff(nat_m_dev,difference=2))
  
  obj_fun <- num_like + ret_like + tot_like + num_sc_like + ret_sc_like + 
    tot_sc_like + nat_m_like + nat_m_prior + f_like + rec_like + sig_m_like + f_prior+
      smooth_f + smooth_r + smooth_m
  
  REPORT(num_like)
  REPORT(ret_like)
  REPORT(tot_like)
  REPORT(num_sc_like)
  REPORT(ret_sc_like)
  REPORT(tot_sc_like)
  REPORT(nat_m_like)
  REPORT(nat_m_prior)
  REPORT(f_like)  
  REPORT(rec_like)   
  
  REPORT(tot_n_pred)
  REPORT(ret_n_pred)
  REPORT(nat_m)
  REPORT(f_mort)
  REPORT(tot_rec)
  REPORT(surv_sel)
  REPORT(ret_fish_sel)
  REPORT(tot_fish_sel)  
  
  return(obj_fun)
  
  #==reports
  
}


#==============================
# rec devs + nat m as random effect
map<-list(sigma_f=factor(NA),
          fish_ret_sel_50=factor(NA),
          fish_ret_sel_slope=factor(NA),
          fish_tot_sel_50=factor(NA),
          fish_tot_sel_slope=factor(NA),
          log_m_mu=factor(NA),
          sigma_m=factor(NA))

obj <- MakeADFun(function(p)do.call(pop_dy,p), params, DLL="bbrkc", map=map,random=c("rec_dev","nat_m_dev"))
lower <- obj$par*0-Inf
upper <- obj$par*0+Inf
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper,control=list(eval.max=10000,iter.max=10000)))
system.time(opt <- nlminb(opt$par, obj$fn, obj$gr, lower=lower, upper=upper,control=list(eval.max=10000,iter.max=10000)))
rep <- sdreport(obj)

#==plot fits
plot(data$n_obs~data$years,
     pch=20,ylab="Male abundance",
     xlab="Year",las=1)
lines(obj$report()$tot_n_pred ~ data$years,lwd=2)

plot(data$ret_cat~data$cat_years,
     pch=20,ylab="Retained catch",
     xlab="Year",las=1)
lines(obj$report()$ret_n_pred ~ data$years,lwd=2)

plot(obj$report()$'tot_rec'~data$years,type='l')
plot(obj$report()$'nat_m'~data$years,type='l')
plot(obj$report()$'f_mort'~data$years,type='l')

plot(obj$report()$'surv_sel',type='l',ylim=c(0,1))
lines(obj$report()$'ret_fish_sel',type='l',col=2)
lines(obj$report()$'tot_fish_sel',type='l',col=3)


# Get the standard report
sdr <- sdreport(obj)
exp(-1.789)
# Extract the random effects estimates and standard errors
# The 'summary()' function with select="random" returns a matrix
sdr_summary  <- summary(sdr, select = "random")
# Convert to data frame and capture row names as 'Parameter'
re_df <- as.data.frame(sdr_summary)
re_df$Parameter <- rownames(re_df)
colnames(re_df) <- c("Estimate", "StdError", "Parameter")

# Create a 'Type' column for faceting and an 'Index' for the x-axis
# Note: This assumes they appear in blocks (e.g., all rec_dev then all nat_m_dev)
re_df$Type <- ifelse(grepl("rec_dev", re_df$Parameter), "Recruitment Deviations", "Natural Mortality Deviations")

# Add a time index within each group (Type)
re_df <- do.call(rbind, lapply(split(re_df, re_df$Type), function(d) {
  d$Year <- 1:nrow(d) # Replace with your actual year vector if available
  return(d)
}))

# Calculate confidence intervals for ribbons
re_df$lwr <- exp(-1.789 + (re_df$Estimate - 1.96 * re_df$StdError))
re_df$upr <- exp(-1.789+ (re_df$Estimate + 1.96 * re_df$StdError))

sdr['log_avg_rec',]

# 1. Get the summary for all fixed effects
fixed_effects <- summary(sdreport(obj), select = "fixed")
my_param <- fixed_effects["log_avg_rec", ]
param_value <- fixed_effects["log_avg_rec", "Estimate"]

library(ggplot2)

ggplot(re_df) +
  geom_ribbon(aes(x = Year,ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~Type, scales = "free_y", ncol = 1) +
  labs(title = "Random Effect Deviations with 95% Uncertainty",
       y = "Estimate (log-scale)", x = "Year/Index") +
  theme_minimal()

