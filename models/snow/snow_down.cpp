#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <snow_down.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  styr.allocate("styr");
  endyr.allocate("endyr");
  dat_yr.allocate("dat_yr");
  years.allocate(1,dat_yr,"years");
  size_n.allocate("size_n");
  sizes.allocate(1,size_n,"sizes");
  imm_n_obs.allocate(styr,endyr,"imm_n_obs");
  imm_n_size_obs.allocate(styr,endyr,1,size_n,"imm_n_size_obs");
  mat_n_obs.allocate(styr,endyr,"mat_n_obs");
  mat_n_size_obs.allocate(styr,endyr,1,size_n,"mat_n_size_obs");
  prop_term_molt.allocate(styr,endyr,1,size_n,"prop_term_molt");
  size_trans.allocate(1,size_n,1,size_n,"size_trans");
  survey_sel.allocate(1,size_n,"survey_sel");
  sigma_numbers_imm.allocate(styr,endyr,"sigma_numbers_imm");
  sigma_numbers_mat.allocate(styr,endyr,"sigma_numbers_mat");
  mat_eff_samp.allocate("mat_eff_samp");
  imm_eff_samp.allocate("imm_eff_samp");
  log_mu_m_prior.allocate(1,3,"log_mu_m_prior");
  est_m_devs.allocate("est_m_devs");
  est_q_devs.allocate("est_q_devs");
  est_m_mat_devs.allocate("est_m_mat_devs");
  est_q_mat_devs.allocate("est_q_mat_devs");
  est_sigma_m.allocate("est_sigma_m");
  sigma_m_mu.allocate(1,3,"sigma_m_mu");
  smooth_q_weight.allocate("smooth_q_weight");
  smooth_m_weight.allocate(1,2,"smooth_m_weight");
  est_log_m_mu.allocate("est_log_m_mu");
  est_sigma_q.allocate("est_sigma_q");
  est_m_lg_devs.allocate("est_m_lg_devs");
  smooth_f_weight.allocate("smooth_f_weight");
  surv_sel_cv.allocate("surv_sel_cv");
  surv_sel_cv_2.allocate("surv_sel_cv_2");
  smooth_surv_weight.allocate("smooth_surv_weight");
  est_sel.allocate("est_sel");
  large_cutoff.allocate("large_cutoff");
  est_tv_fish_sel.allocate("est_tv_fish_sel");
cout<<"imm_n_obs"<<imm_n_obs<<endl; 
cout<<"mat_n_obs"<<mat_n_obs<<endl; 
cout<<"survey_sel"<<survey_sel<<endl; 
cout<<"est_m_devs"<<est_m_devs<<endl;
cout<<"sigma_m_mu"<<sigma_m_mu<<endl;
cout<<"smooth_m_weight"<<smooth_m_weight<<endl;
cout<<"est_m_lg_devs"<<est_m_lg_devs<<endl;
 ad_comm::change_datafile_name("catch_dat.DAT");
  ret_cat_yr_n.allocate("ret_cat_yr_n");
  ret_cat_yrs.allocate(1,ret_cat_yr_n,"ret_cat_yrs");
  ret_cat_numbers.allocate(1,ret_cat_yr_n,"ret_cat_numbers");
  disc_cat_numbers.allocate(1,ret_cat_yr_n,"disc_cat_numbers");
  sigma_numbers_ret.allocate("sigma_numbers_ret");
  sigma_numbers_disc.allocate("sigma_numbers_disc");
  discard_survival.allocate("discard_survival");
  ret_eff_samp.allocate("ret_eff_samp");
  disc_eff_samp.allocate("disc_eff_samp");
  ret_cat_size.allocate(1,ret_cat_yr_n,1,size_n,"ret_cat_size");
  disc_cat_size.allocate(1,ret_cat_yr_n,1,size_n,"disc_cat_size");
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_n_imm.allocate(1,size_n,1.01,30,1,"log_n_imm");
  log_n_mat.allocate(1,size_n,1.01,30,1,"log_n_mat");
  nat_m_dev.allocate(styr,endyr,-4,4,est_m_devs,"nat_m_dev");
  nat_m_mat_dev.allocate(styr,endyr,-4,4,est_m_mat_devs,"nat_m_mat_dev");
  nat_m_lg_dev.allocate(styr,endyr,-4,4,est_m_mat_devs,"nat_m_lg_dev");
  q_dev.allocate(styr,endyr,-0.2,0.2,est_q_devs,"q_dev");
  q_mat_dev.allocate(styr,endyr,-0.2,0.2,est_q_mat_devs,"q_mat_dev");
  q_lg_dev.allocate(styr,endyr,-0.2,0.2,est_q_mat_devs,"q_lg_dev");
  log_avg_rec.allocate(1,40,"log_avg_rec");
  rec_devs.allocate(styr,endyr,-10,10,1,"rec_devs");
  sigma_m.allocate(1,3,0.01,4,est_sigma_m,"sigma_m");
  log_m_mu.allocate(1,3,-5,3,est_log_m_mu,"log_m_mu");
  prop_rec.allocate(1,2,0.00001,200,"prop_rec");
  sigma_q.allocate(1,2,0.01,4,est_sigma_q,"sigma_q");
  log_f.allocate(-5,5,"log_f");
  f_dev.allocate(1,ret_cat_yr_n,-5,5,"f_dev");
  fish_ret_sel_50.allocate(25,150,"fish_ret_sel_50");
  fish_ret_sel_50_post.allocate(25,150,"fish_ret_sel_50_post");
  fish_ret_sel_slope.allocate(0.0001,20,"fish_ret_sel_slope");
  fish_ret_sel_slope_post.allocate(0.0001,20,"fish_ret_sel_slope_post");
  fish_tot_sel_offset.allocate(-10,20,"fish_tot_sel_offset");
  fish_tot_sel_offset_dev.allocate(1,ret_cat_yr_n,-40,40,"fish_tot_sel_offset_dev");
  fish_tot_sel_slope.allocate(-10,0,"fish_tot_sel_slope");
  surv_omega.allocate(0,1,est_sel,"surv_omega");
  surv_alpha1.allocate(0.0001,20,est_sel,"surv_alpha1");
  surv_beta1.allocate(25,150,est_sel,"surv_beta1");
  surv_alpha2.allocate(0.0001,20,est_sel,"surv_alpha2");
  surv_beta2.allocate(25,150,est_sel,"surv_beta2");
  imm_n_size_pred.allocate(styr,endyr,1,size_n,"imm_n_size_pred");
  #ifndef NO_AD_INITIALIZE
    imm_n_size_pred.initialize();
  #endif
  mat_n_size_pred.allocate(styr,endyr,1,size_n,"mat_n_size_pred");
  #ifndef NO_AD_INITIALIZE
    mat_n_size_pred.initialize();
  #endif
  nat_m.allocate(styr,endyr,1,size_n,"nat_m");
  #ifndef NO_AD_INITIALIZE
    nat_m.initialize();
  #endif
  nat_m_mat.allocate(styr,endyr,1,size_n,"nat_m_mat");
  #ifndef NO_AD_INITIALIZE
    nat_m_mat.initialize();
  #endif
  selectivity.allocate(styr,endyr,1,size_n,"selectivity");
  #ifndef NO_AD_INITIALIZE
    selectivity.initialize();
  #endif
  selectivity_mat.allocate(styr,endyr,1,size_n,"selectivity_mat");
  #ifndef NO_AD_INITIALIZE
    selectivity_mat.initialize();
  #endif
  total_fish_sel.allocate(styr,endyr,1,size_n,"total_fish_sel");
  #ifndef NO_AD_INITIALIZE
    total_fish_sel.initialize();
  #endif
  retain_fish_sel.allocate(styr,endyr,1,size_n,"retain_fish_sel");
  #ifndef NO_AD_INITIALIZE
    retain_fish_sel.initialize();
  #endif
  pred_retained_n.allocate(styr,endyr,"pred_retained_n");
  #ifndef NO_AD_INITIALIZE
    pred_retained_n.initialize();
  #endif
  pred_discard_n.allocate(styr,endyr,"pred_discard_n");
  #ifndef NO_AD_INITIALIZE
    pred_discard_n.initialize();
  #endif
  pred_retained_size_comp.allocate(styr,endyr,1,size_n,"pred_retained_size_comp");
  #ifndef NO_AD_INITIALIZE
    pred_retained_size_comp.initialize();
  #endif
  pred_discard_size_comp.allocate(styr,endyr,1,size_n,"pred_discard_size_comp");
  #ifndef NO_AD_INITIALIZE
    pred_discard_size_comp.initialize();
  #endif
  f_mort.allocate(styr,endyr,"f_mort");
  #ifndef NO_AD_INITIALIZE
    f_mort.initialize();
  #endif
  temp_imm.allocate(1,size_n,"temp_imm");
  #ifndef NO_AD_INITIALIZE
    temp_imm.initialize();
  #endif
  temp_mat.allocate(1,size_n,"temp_mat");
  #ifndef NO_AD_INITIALIZE
    temp_mat.initialize();
  #endif
  trans_imm.allocate(1,size_n,"trans_imm");
  #ifndef NO_AD_INITIALIZE
    trans_imm.initialize();
  #endif
  temp_catch_imm.allocate(1,size_n,"temp_catch_imm");
  #ifndef NO_AD_INITIALIZE
    temp_catch_imm.initialize();
  #endif
  temp_catch_mat.allocate(1,size_n,"temp_catch_mat");
  #ifndef NO_AD_INITIALIZE
    temp_catch_mat.initialize();
  #endif
  surv_sel.allocate(1,size_n,"surv_sel");
  #ifndef NO_AD_INITIALIZE
    surv_sel.initialize();
  #endif
  sum_imm_numbers_obs.allocate(styr,endyr,"sum_imm_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_imm_numbers_obs.initialize();
  #endif
  sum_mat_numbers_obs.allocate(styr,endyr,"sum_mat_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_mat_numbers_obs.initialize();
  #endif
  imm_numbers_pred.allocate(styr,endyr,"imm_numbers_pred");
  #ifndef NO_AD_INITIALIZE
    imm_numbers_pred.initialize();
  #endif
  mat_numbers_pred.allocate(styr,endyr,"mat_numbers_pred");
  #ifndef NO_AD_INITIALIZE
    mat_numbers_pred.initialize();
  #endif
  sum_ret_numbers_obs.allocate(styr,endyr,"sum_ret_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_ret_numbers_obs.initialize();
  #endif
  sum_disc_numbers_obs.allocate(styr,endyr,"sum_disc_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_disc_numbers_obs.initialize();
  #endif
  total_population_n.allocate(styr,endyr,"total_population_n");
  fished_population_n.allocate(styr,endyr,"fished_population_n");
  imm_num_like.allocate("imm_num_like");
  #ifndef NO_AD_INITIALIZE
  imm_num_like.initialize();
  #endif
  mat_num_like.allocate("mat_num_like");
  #ifndef NO_AD_INITIALIZE
  mat_num_like.initialize();
  #endif
  ret_cat_like.allocate("ret_cat_like");
  #ifndef NO_AD_INITIALIZE
  ret_cat_like.initialize();
  #endif
  disc_cat_like.allocate("disc_cat_like");
  #ifndef NO_AD_INITIALIZE
  disc_cat_like.initialize();
  #endif
  imm_like.allocate("imm_like");
  #ifndef NO_AD_INITIALIZE
  imm_like.initialize();
  #endif
  mat_like.allocate("mat_like");
  #ifndef NO_AD_INITIALIZE
  mat_like.initialize();
  #endif
  ret_comp_like.allocate("ret_comp_like");
  #ifndef NO_AD_INITIALIZE
  ret_comp_like.initialize();
  #endif
  disc_comp_like.allocate("disc_comp_like");
  #ifndef NO_AD_INITIALIZE
  disc_comp_like.initialize();
  #endif
  use_term_molt.allocate(styr,endyr,1,size_n,"use_term_molt");
  #ifndef NO_AD_INITIALIZE
    use_term_molt.initialize();
  #endif
  nat_m_like.allocate("nat_m_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_like.initialize();
  #endif
  nat_m_mat_like.allocate("nat_m_mat_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_mat_like.initialize();
  #endif
  nat_m_mu_like.allocate("nat_m_mu_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_mu_like.initialize();
  #endif
  nat_m_lg_like.allocate("nat_m_lg_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_lg_like.initialize();
  #endif
  nat_m_mat_mu_like.allocate("nat_m_mat_mu_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_mat_mu_like.initialize();
  #endif
  nat_m_lg_mu_like.allocate("nat_m_lg_mu_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_lg_mu_like.initialize();
  #endif
  smooth_q_like.allocate("smooth_q_like");
  #ifndef NO_AD_INITIALIZE
  smooth_q_like.initialize();
  #endif
  smooth_m_like.allocate("smooth_m_like");
  #ifndef NO_AD_INITIALIZE
  smooth_m_like.initialize();
  #endif
  q_like.allocate("q_like");
  #ifndef NO_AD_INITIALIZE
  q_like.initialize();
  #endif
  q_mat_like.allocate("q_mat_like");
  #ifndef NO_AD_INITIALIZE
  q_mat_like.initialize();
  #endif
  smooth_f_like.allocate("smooth_f_like");
  #ifndef NO_AD_INITIALIZE
  smooth_f_like.initialize();
  #endif
  surv_sel_prior.allocate("surv_sel_prior");
  #ifndef NO_AD_INITIALIZE
  surv_sel_prior.initialize();
  #endif
  smooth_surv_like.allocate("smooth_surv_like");
  #ifndef NO_AD_INITIALIZE
  smooth_surv_like.initialize();
  #endif
  f_prior.allocate("f_prior");
  #ifndef NO_AD_INITIALIZE
  f_prior.initialize();
  #endif
  temp_prop_rec.allocate(1,3,"temp_prop_rec");
  #ifndef NO_AD_INITIALIZE
    temp_prop_rec.initialize();
  #endif
  tot_prop_rec.allocate("tot_prop_rec");
  #ifndef NO_AD_INITIALIZE
  tot_prop_rec.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  f =0.0;
 dvariable fpen=0.0;
 for(int size=1;size<=size_n;size++)
  {
   imm_n_size_pred(styr,size) = exp(log_n_imm(size));
   mat_n_size_pred(styr,size) = exp(log_n_mat(size));
   surv_sel(size) = (surv_omega / (1 + exp(-surv_alpha1*(sizes(size)-surv_beta1))) ) + ( (1-surv_omega)/(1+exp(-surv_alpha2*(sizes(size)-surv_beta2)))  );
   }
  // total selectivity  
  for(int year=1;year<=ret_cat_yr_n;year++)
  for(int size=1;size<=size_n;size++)
   {
	total_fish_sel(ret_cat_yrs(year),size) = 1 / (1+exp(-exp(fish_tot_sel_slope)*(sizes(size)-(fish_ret_sel_50-fish_tot_sel_offset)))) ; 
    if(est_tv_fish_sel==1)
	{	total_fish_sel(ret_cat_yrs(year),size) = 1 / (1+exp(-exp(fish_tot_sel_slope)*(sizes(size)-(fish_ret_sel_50-(fish_tot_sel_offset +fish_tot_sel_offset_dev(year)))))) ; }
 	//{	total_fish_sel(ret_cat_yrs(year),size) = 1 / (1+exp(-exp(fish_tot_sel_slope+fish_tot_sel_slope_dev(year))*(sizes(size)-(fish_ret_sel_50-(fish_tot_sel_offset +fish_tot_sel_offset_dev(year)))))) ; }
	}
 // retained selectivity; changes after rationalization
  for(int year=styr;year<=endyr;year++)
  for(int size=1;size<=size_n;size++)
   {
	if(year>2004)
	 retain_fish_sel(year,size) = 1 / (1+exp(-fish_ret_sel_slope_post*(sizes(size)-fish_ret_sel_50_post)));
	if(year<=2004)
     retain_fish_sel(year,size) = 1 / (1+exp(-fish_ret_sel_slope*(sizes(size)-fish_ret_sel_50)));
   }	  
   nat_m_dev(2020) = 0;
   nat_m_mat_dev(2020) = 0;
 for(int year=styr;year<=endyr;year++)
  for(int size=1;size<=size_n;size++)
  {
  nat_m(year,size) = log_m_mu(1) + nat_m_dev(year);
  nat_m_mat(year,size) = log_m_mu(2) + nat_m_dev(year);
  selectivity(year,size) = surv_sel(size) + q_dev(year);
  selectivity_mat(year,size) = surv_sel(size) + q_dev(year);
  //==option for estimating separate devs by maturity state
  if(est_q_mat_devs>0)
   selectivity_mat(year,size) = surv_sel(size) + q_mat_dev(year);  
  if(est_m_mat_devs>0)
   nat_m_mat(year,size) = log_m_mu(2) + nat_m_mat_dev(year);
  if(est_m_lg_devs>0 & sizes(size) > large_cutoff)
  {
   nat_m(year,size) = log_m_mu(3) + nat_m_lg_dev(year);
   nat_m_mat(year,size) = log_m_mu(3) + nat_m_lg_dev(year);
  }
  }
  temp_prop_rec.initialize();
  tot_prop_rec = 20;
  for(int i=1;i<=2;i++)
   tot_prop_rec += prop_rec(i);
   temp_prop_rec(1) = 20/ tot_prop_rec;
  for(int i = 1;i<=2;i++)
   temp_prop_rec(i+1) = prop_rec(i)/tot_prop_rec;
 // use_term_molt = prop_term_molt;
 // for(int size = 1; size<=size_n;size++)
 //  use_term_molt(2020,size) = prob_mat_2020(size);
  use_term_molt = prop_term_molt;
  //for(int size = 1; size<=16;size++)
  // use_term_molt(2020,size) = 1 / (1+exp(-prob_mat_2020(size)));
 //==make last year rec dev and f dev equal to the average for now
 rec_devs(endyr) =0;
 f_mort.initialize();
  for(int year=1;year<=ret_cat_yr_n;year++)
   f_mort(ret_cat_yrs(year)) = exp(log_f+f_dev(year));
 for(int year=styr;year<endyr;year++)
  {
  for (int size=1;size<=size_n;size++) 
      {
		temp_imm(size) = imm_n_size_pred(year,size) * exp(-1*(0.59)*exp( nat_m(year,size)));
		temp_mat(size) = mat_n_size_pred(year,size) * exp(-1*(0.59)*exp( nat_m_mat(year,size)));
	   }	
  // fishery
	   for(int size=1;size<=size_n;size++)
	   {	   
	   temp_catch_imm(size) = temp_imm(size) * (1 -exp(-(f_mort(year))*total_fish_sel(year,size)));
	   temp_catch_mat(size) = temp_mat(size) * (1 -exp(-(f_mort(year))*total_fish_sel(year,size)));
	   pred_retained_size_comp(year,size) = temp_catch_imm(size)*retain_fish_sel(year,size) + temp_catch_mat(size)*retain_fish_sel(year,size); 
	   pred_discard_size_comp(year,size) = temp_catch_imm(size)*(1-retain_fish_sel(year,size)) + temp_catch_mat(size)*(1-retain_fish_sel(year,size));
	   temp_imm(size) = temp_imm(size) * exp(-(f_mort(year))*total_fish_sel(year,size));
	   temp_mat(size) = temp_mat(size) * exp(-(f_mort(year))*total_fish_sel(year,size));
	   temp_imm(size) += temp_catch_imm(size)*(1-retain_fish_sel(year,size))*(discard_survival);
	   temp_mat(size) += temp_catch_mat(size)*(1-retain_fish_sel(year,size))*(discard_survival);
	   }	
	  // growth
	   trans_imm = size_trans * temp_imm;
	   // recruitment
       trans_imm(1) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(1);
       trans_imm(2) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(2);
	   trans_imm(3) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(3);
	  // maturity
	   for (int size=1;size<=size_n;size++) 
	   {
	  	temp_imm(size) = (trans_imm(size) * (1-use_term_molt(year,size))); 	
	  	temp_mat(size) = (trans_imm(size) * use_term_molt(year,size) + temp_mat(size));
       }
      // natural mortality		
	  for (int size=1;size<=size_n;size++) 
	   {
	   imm_n_size_pred(year+1,size) = temp_imm(size) * exp(-1*(0.41)*exp(nat_m(year,size)));
	   mat_n_size_pred(year+1,size) = temp_mat(size) * exp(-1*(0.41)*exp(nat_m_mat(year,size)));
       }
    }
   evaluate_the_objective_function();
}

void model_parameters::evaluate_the_objective_function(void)
{
  // make total numbers by maturity state from obs and preds
  imm_numbers_pred.initialize();
  mat_numbers_pred.initialize();
  sum_imm_numbers_obs.initialize();
  sum_mat_numbers_obs.initialize(); 
  pred_retained_n.initialize();
  pred_discard_n.initialize();
  total_population_n.initialize();
  fished_population_n.initialize();
  for (int year=styr;year<=endyr;year++)
   for (int size=1;size<=size_n;size++)
   {
    imm_numbers_pred(year)    += selectivity(year,size)*imm_n_size_pred(year,size);
    mat_numbers_pred(year)    += selectivity_mat(year,size)*mat_n_size_pred(year,size);
	sum_imm_numbers_obs(year) += imm_n_size_obs(year,size);
	sum_mat_numbers_obs(year) += mat_n_size_obs(year,size);
	pred_retained_n(year)     += pred_retained_size_comp(year,size);
	pred_discard_n(year)      += pred_discard_size_comp(year,size);
	total_population_n(year)  += imm_n_size_pred(year,size)+mat_n_size_pred(year,size);
	fished_population_n(year)  += retain_fish_sel(year,size)*imm_n_size_pred(year,size)+retain_fish_sel(year,size)*mat_n_size_pred(year,size);
	   }
  // likelihoods
  imm_num_like = 0;
  for (int year=styr;year<=endyr;year++)
   if (year!=2020)
    imm_num_like += square( log(imm_numbers_pred(year)) - log(imm_n_obs(year))) / (2.0 * sqrt(log(1+square(sigma_numbers_imm(year)))));
  mat_num_like = 0;
  for (int year=styr;year<=endyr;year++)
   if (year!=2020)
    mat_num_like += square( log(mat_numbers_pred(year)) - log(mat_n_obs(year))) / (2.0 * sqrt(log(1+square(sigma_numbers_mat(year)))));
   // likelihoods
  ret_cat_like = 0;
  for (int year=1;year<=ret_cat_yr_n;year++)
    ret_cat_like += square( log(pred_retained_n(ret_cat_yrs(year))) - log(ret_cat_numbers(year))) / (2.0 * sqrt(log(1+square(sigma_numbers_ret))));
  disc_cat_like = 0;
  for (int year=1;year<=ret_cat_yr_n;year++)
    disc_cat_like += square( log(pred_discard_n(ret_cat_yrs(year))) - log(disc_cat_numbers(year))) / (2.0 * sqrt(log(1+square(sigma_numbers_disc))));
  // immature numbers at size data
  imm_like = 0;
  for (int year=styr;year<=endyr;year++)
   for (int size=1;size<=size_n;size++)
    if (imm_n_size_obs(year,size) >0 & imm_n_size_pred(year,size) >0)
     imm_like += imm_eff_samp*(imm_n_size_obs(year,size)/sum_imm_numbers_obs(year)) * log( (selectivity(year,size)*imm_n_size_pred(year,size)/imm_numbers_pred(year)) / (imm_n_size_obs(year,size)/sum_imm_numbers_obs(year)));
  imm_like = -1*imm_like;
  // mature numbers at size data
  mat_like = 0;
  for (int year=styr;year<=endyr;year++)
   for (int size=1;size<=size_n;size++)
    if (mat_n_size_obs(year,size) >0)
     mat_like += mat_eff_samp*(mat_n_size_obs(year,size)/sum_mat_numbers_obs(year)) * log( (selectivity_mat(year,size)*mat_n_size_pred(year,size)/mat_numbers_pred(year)) / (mat_n_size_obs(year,size)/sum_mat_numbers_obs(year)));
  mat_like = -1*mat_like;
  // retained catch at size data
  ret_comp_like = 0;
  for (int year=1;year<=ret_cat_yr_n;year++)
   for (int size=1;size<=size_n;size++)
    if (ret_cat_size(year,size) >0.001 & pred_retained_size_comp(ret_cat_yrs(year),size) >0.001)
     ret_comp_like += ret_eff_samp*(ret_cat_size(year,size)) * log( (pred_retained_size_comp(ret_cat_yrs(year),size)/pred_retained_n(ret_cat_yrs(year))) / (ret_cat_size(year,size)));
  ret_comp_like = -1*ret_comp_like;
  // discard catch at size data
  disc_comp_like = 0;
  for (int year=1;year<=ret_cat_yr_n;year++)
   for (int size=1;size<=size_n;size++)
    if (disc_cat_size(year,size) >0.001 & pred_discard_size_comp(ret_cat_yrs(year),size) >0.001)
     disc_comp_like += disc_eff_samp*(disc_cat_size(year,size)) * log( (pred_discard_size_comp(ret_cat_yrs(year),size)/pred_discard_n(ret_cat_yrs(year))) / (disc_cat_size(year,size)));
  disc_comp_like = -1*disc_comp_like;
 //penalties on m 
  nat_m_mu_like =0;
  nat_m_mu_like += pow(((log_m_mu(1))-(log_mu_m_prior(1)))/ (sqrt(2)*sqrt(sigma_m_mu(1))),2.0);
  nat_m_mat_mu_like =0;
  nat_m_mat_mu_like += pow(((log_m_mu(2))-(log_mu_m_prior(2)))/ (sqrt(2)*sqrt(sigma_m_mu(2))),2.0); 
  surv_sel_prior.initialize();
  if(current_phase()==est_sel)
  {
  for (int size=1;size<14;size++)
   surv_sel_prior +=  pow(((survey_sel(size))-(surv_sel(size)))/ (sqrt(2)*sqrt(surv_sel_cv)),2.0); 
  for (int size=14;size<=size_n;size++)
   surv_sel_prior +=  pow(((survey_sel(size))-(surv_sel(size)))/ (sqrt(2)*sqrt(surv_sel_cv_2)),2.0); 
  }
  f_prior.initialize();  
  if(est_tv_fish_sel==1)
  {
   for (int year=styr;year<=endyr;year++)
    f_prior += pow((f_mort(year)-1)/ (sqrt(2)*sqrt(50)),2.0);
   }
  if(est_m_devs>0)
  {
  nat_m_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_like += pow((nat_m(year,1)-log_m_mu(1))/ (sqrt(2)*sqrt(sigma_m(1))),2.0);
  nat_m_mat_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_mat_like += pow((nat_m_mat(year,1)-log_m_mu(2))/ (sqrt(2)*sqrt(sigma_m(2))),2.0);
   nat_m_lg_like =0;
   if(est_m_lg_devs>0)
  {
  for (int year=styr;year<=endyr;year++)
   nat_m_lg_like += pow((nat_m_mat(year,20)-log_m_mu(3))/ (sqrt(2)*sqrt(sigma_m(3))),2.0);
  }
  }
  if(est_q_devs>0)
  {
  q_like =0;
  for (int year=styr;year<=endyr;year++)
   q_like += pow((selectivity(year,4)-surv_sel(4))/ (sqrt(2)*sqrt(sigma_q(1))),2.0);
  q_mat_like =0;
  for (int year=styr;year<=endyr;year++)
   q_mat_like += pow((selectivity_mat(year,4)-surv_sel(4))/ (sqrt(2)*sqrt(sigma_q(2))),2.0);
  }
  smooth_q_like = 0;
  smooth_q_like = smooth_q_weight* (norm2(first_difference(first_difference(q_dev))) +norm2(first_difference(first_difference(q_mat_dev)))) ;
  smooth_m_like = 0;
  smooth_m_like = smooth_m_weight(1)* norm2(first_difference(first_difference(nat_m_dev))) + smooth_m_weight(2)*norm2(first_difference(first_difference(nat_m_mat_dev))) ;
 // this is a smoothness penalty between immature and mature 
  //for (int year=styr;year<=endyr;year++)
   //smooth_m_like += pow((nat_m_mat(year,1)-nat_m(year,1))/ (sqrt(2)*sqrt(sigma_m(3))),2.0);
  smooth_f_like = 0;
  smooth_f_like = smooth_f_weight* (norm2(first_difference(first_difference(f_dev))));
  smooth_surv_like = 0;
  //smooth_surv_like = smooth_surv_weight* (norm2(first_difference(first_difference(surv_sel))));
  f = imm_num_like + mat_num_like + ret_cat_like + disc_cat_like + imm_like + mat_like + ret_comp_like + disc_comp_like + 
  nat_m_like + nat_m_mat_like + nat_m_lg_like + nat_m_mu_like + nat_m_mat_mu_like + smooth_q_like + smooth_m_like + q_like + q_mat_like + smooth_f_like +
  surv_sel_prior + smooth_surv_like + f_prior;
  cout<<imm_num_like<< " " << mat_num_like << " " << imm_like << " " << mat_like << " " <<endl;
  cout<<ret_cat_like<< " " << ret_cat_like << " " << ret_comp_like << " " << disc_comp_like << " " <<endl;
  cout<<nat_m_like<< " " << disc_cat_like << " " << nat_m_lg_like << " " << nat_m_mu_like << " " << nat_m_mat_mu_like << " " <<endl;
  cout<<smooth_q_like<< " " << smooth_m_like << " " << q_like << " " << q_mat_like << " " <<smooth_f_like<<" " << surv_sel_prior<<" "<<f_prior<<" "<<endl;
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report<<"$likelihoods"<<endl;
  report<<imm_num_like<<" "<<mat_num_like<<" "<<imm_like<<" "<<mat_like<<" "<<endl;
  report<<ret_cat_like<< " " << disc_cat_like << " " << ret_comp_like << " " << disc_comp_like << " " <<endl;
  report<<"$penalties"<<endl;
  report<<nat_m_like<< " " << nat_m_mat_like << " " << nat_m_lg_like << " " << nat_m_mu_like << " " << nat_m_mat_mu_like << " " <<endl;
  report<<smooth_q_like<< " " << q_like << " " << q_mat_like << " " <<smooth_f_like<<" " << " " <<smooth_m_like<<" " <<endl;
  report<<surv_sel_prior<<" "<<"place2"<<" "<<surv_sel_prior<<" "<<smooth_surv_like<<" "<<f_prior<<" "<<endl;  
  report <<"$natural mortality" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << mfexp(nat_m(i))<<endl;
  }
  report <<"$size_trans" << endl;
  for(int i=1; i<=size_n; i++)
  {
    report << size_trans(i)<<endl;
  }
  report <<"$recruits" << endl;
  for(int i=styr; i<=endyr; i++)
   report << mfexp(log_avg_rec + rec_devs(i))<<endl;
  report <<"$immature numbers at size" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (elem_prod(selectivity(i),imm_n_size_pred(i)))/imm_numbers_pred(i)<<endl;
  }
  report <<"$mature numbers at size" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (elem_prod(selectivity_mat(i),mat_n_size_pred(i)))/mat_numbers_pred(i)<<endl;
  }
  report <<"$obs_imm_n_size" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << imm_n_size_obs(i)<<endl;
  }
  report <<"$obs_mat_n_size" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << mat_n_size_obs(i)<<endl;
  }
  report<<"$styr"<<endl;
  report<<styr<<endl;
  report<<"$endyr"<<endl;
  report<<endyr<<endl;
  report<<"$imm_numbers_pred"<<endl;
  report<<imm_numbers_pred<<endl;
  report<<"$mat_numbers_pred"<<endl;
  report<<mat_numbers_pred<<endl;
  report<<"$imm_n_obs"<<endl;
  report<<imm_n_obs<<endl;
  report<<"$mat_n_obs"<<endl;
  report<<mat_n_obs<<endl;
  report<<"$survey_sel_input"<<endl;
  report<<surv_sel<<endl;
  report<<"$ret_fish_sel"<<endl;
  for(int i=styr; i<=endyr; i++)
  report<<retain_fish_sel(i)<<endl;
  report<<"$total_fish_sel"<<endl;
  for(int i=styr; i<=endyr; i++)
  report<<total_fish_sel(i)<<endl;  
  report<<"$surv_sel"<<endl;
  report<<surv_sel<<endl;  
  report <<"$est_fishing_mort" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << f_mort(i) <<endl;
  }  
  report <<"$mature natural mortality" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << mfexp(nat_m_mat(i))<<endl;
  }
  report <<"$survey selectivity" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << selectivity(i)<<endl;
  }
  report <<"$mature survey selectivity" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << selectivity_mat(i)<<endl;
  }
  report <<"$pred_imm_pop_num" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (imm_n_size_pred(i))<<endl;
  }
  report <<"$pred_mat_pop_num" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (mat_n_size_pred(i))<<endl;
  }
  cout<<3<<endl;
  report <<"$pred_retained_n" << endl;
  for(int i=styr; i<endyr; i++)
  {
    report << (pred_retained_n(i))<<endl;
  }
  report <<"$pred_discard_n" << endl;
  for(int i=styr; i<endyr; i++)
  {
    report << (pred_discard_n(i))<<endl;
  }
  report <<"$ret_cat_numbers" << endl;
  for(int i=1; i<=ret_cat_yr_n; i++)
  {
    report << (ret_cat_numbers(i))<<endl;
  }
  report <<"$disc_cat_numbers" << endl;
  for(int i=1; i<=ret_cat_yr_n; i++)
  {
    report << (disc_cat_numbers(i))<<endl;
  }
  cout<<4<<endl;
  report <<"$pred_retained_size_comp" << endl;
  for(int i=styr; i<endyr; i++)
  {
    report << ((pred_retained_size_comp(i)/pred_retained_n(i)))<<endl;
  }
  report <<"$pred_discard_size_comp" << endl;
  for(int i=styr; i<endyr; i++)
  {
    report << ((pred_discard_size_comp(i)/pred_discard_n(i)))<<endl;
  }
  report <<"$obs_retained_size_comp" << endl;
  for(int i=1; i<=ret_cat_yr_n; i++)
  {
    report << (ret_cat_size(i))<<endl;
  }
  report <<"$obs_discard_size_comp" << endl;
  for(int i=1; i<=ret_cat_yr_n; i++)
  {
    report << (disc_cat_size(i))<<endl;
  }
 	  cout<<5<<endl; 
  report <<"$temp_prop_rec" << endl;
  report << temp_prop_rec << endl;
  report <<"$prob_term_molt" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (use_term_molt(i))<<endl;
  }
  report <<"$sizes" << endl;
  report << sizes << endl;	
  report <<"$imm_cv" << endl;
  report << sigma_numbers_imm << endl;	
  report <<"$mat_cv" << endl;
  report << sigma_numbers_mat << endl;	
  report <<"$ret_cat_yrs" << endl;
  report << ret_cat_yrs << endl;	 
    save_gradients(gradients);
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{10000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-3}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
