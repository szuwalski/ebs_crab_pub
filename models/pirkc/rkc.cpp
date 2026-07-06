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
#include <rkc.htp>

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
  n_obs.allocate(styr,endyr,"n_obs");
  n_size_obs.allocate(styr,endyr,1,size_n,"n_size_obs");
  prob_molt.allocate(1,size_n,"prob_molt");
  size_trans.allocate(1,size_n,1,size_n,"size_trans");
  sigma_numbers.allocate(styr,endyr,"sigma_numbers");
  sc_eff_samp.allocate("sc_eff_samp");
  log_mu_m_prior.allocate("log_mu_m_prior");
  est_m_devs.allocate("est_m_devs");
  est_q_devs.allocate("est_q_devs");
  est_sigma_m.allocate("est_sigma_m");
  sigma_m_mu.allocate("sigma_m_mu");
  smooth_q_weight.allocate("smooth_q_weight");
  smooth_m_weight.allocate("smooth_m_weight");
  est_log_m_mu.allocate("est_log_m_mu");
  est_sigma_q.allocate("est_sigma_q");
  smooth_f_weight.allocate("smooth_f_weight");
  surv_sel_cv.allocate("surv_sel_cv");
  smooth_surv_weight.allocate("smooth_surv_weight");
  smooth_init_weight.allocate("smooth_init_weight");
  est_molt.allocate("est_molt");
cout<<"n_obs"<<n_obs<<endl; 
cout<<"est_m_devs"<<est_m_devs<<endl;
cout<<"sigma_m_mu"<<sigma_m_mu<<endl;
cout<<"smooth_m_weight"<<smooth_m_weight<<endl;
cout<<"sigma_numbers"<<sigma_numbers<<endl;
 ad_comm::change_datafile_name("catch_dat.DAT");
  ret_cat_yr_n.allocate("ret_cat_yr_n");
  ret_cat_yrs.allocate(1,ret_cat_yr_n,"ret_cat_yrs");
  ret_cat_numbers.allocate(1,ret_cat_yr_n,"ret_cat_numbers");
  ret_sc_yr_n.allocate("ret_sc_yr_n");
  ret_sc_yrs.allocate(1,ret_sc_yr_n,"ret_sc_yrs");
  ret_cat_size.allocate(1,ret_sc_yr_n,1,size_n,"ret_cat_size");
  tot_cat_yr_n.allocate("tot_cat_yr_n");
  tot_cat_yrs.allocate(1,tot_cat_yr_n,"tot_cat_yrs");
  tot_cat_numbers.allocate(1,tot_cat_yr_n,"tot_cat_numbers");
  tot_sc_yr_n.allocate("tot_sc_yr_n");
  tot_sc_yrs.allocate(1,tot_sc_yr_n,"tot_sc_yrs");
  tot_cat_size.allocate(1,tot_sc_yr_n,1,size_n,"tot_cat_size");
  sigma_numbers_ret.allocate("sigma_numbers_ret");
  sigma_numbers_tot.allocate("sigma_numbers_tot");
  discard_survival.allocate("discard_survival");
  ret_eff_samp.allocate("ret_eff_samp");
  tot_eff_samp.allocate("tot_eff_samp");
cout<<"ret_cat_numbers"<<ret_cat_numbers<<endl; 
cout<<"tot_cat_numbers"<<tot_cat_numbers<<endl;
cout<<"tot_sc_yrs"<<tot_sc_yrs<<endl;
cout<<"sigma_numbers_ret"<<sigma_numbers_ret<<endl;
cout<<"tot_eff_samp"<<tot_eff_samp<<endl;
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
  log_n_init.allocate(1,size_n,1.01,30,1,"log_n_init");
  nat_m_dev.allocate(styr,endyr,-4,4,est_m_devs,"nat_m_dev");
  q_dev.allocate(styr,endyr,-0.2,0.2,est_q_devs,"q_dev");
  log_avg_rec.allocate(1,40,"log_avg_rec");
  rec_devs.allocate(styr,endyr,-12,12,1,"rec_devs");
  sigma_m.allocate(0.01,4,est_sigma_m,"sigma_m");
  log_m_mu.allocate(-5,3,est_log_m_mu,"log_m_mu");
  prop_rec.allocate(1,2,0.00001,200,"prop_rec");
  sigma_q.allocate(1,2,0.01,4,est_sigma_q,"sigma_q");
  log_f.allocate(-10,5,"log_f");
  f_dev.allocate(1,ret_cat_yr_n,-5,5,"f_dev");
  fish_ret_sel_50.allocate(25,150,-1,"fish_ret_sel_50");
  fish_ret_sel_slope.allocate(0.0001,20,-1,"fish_ret_sel_slope");
  fish_tot_sel_50.allocate(25,150,"fish_tot_sel_50");
  fish_tot_sel_slope.allocate(0.0001,20,"fish_tot_sel_slope");
  surv_sel_50.allocate(25,150,"surv_sel_50");
  surv_sel_slope.allocate(0.0001,20,"surv_sel_slope");
  molt_sel_50.allocate(25,180,"molt_sel_50");
  molt_sel_slope.allocate(0.0001,20,"molt_sel_slope");
  n_size_pred.allocate(styr,endyr,1,size_n,"n_size_pred");
  #ifndef NO_AD_INITIALIZE
    n_size_pred.initialize();
  #endif
  nat_m.allocate(styr,endyr,1,size_n,"nat_m");
  #ifndef NO_AD_INITIALIZE
    nat_m.initialize();
  #endif
  selectivity.allocate(styr,endyr,1,size_n,"selectivity");
  #ifndef NO_AD_INITIALIZE
    selectivity.initialize();
  #endif
  total_fish_sel.allocate(1,size_n,"total_fish_sel");
  #ifndef NO_AD_INITIALIZE
    total_fish_sel.initialize();
  #endif
  retain_fish_sel.allocate(1,size_n,"retain_fish_sel");
  #ifndef NO_AD_INITIALIZE
    retain_fish_sel.initialize();
  #endif
  surv_sel.allocate(1,size_n,"surv_sel");
  #ifndef NO_AD_INITIALIZE
    surv_sel.initialize();
  #endif
  pred_retained_n.allocate(styr,endyr,"pred_retained_n");
  #ifndef NO_AD_INITIALIZE
    pred_retained_n.initialize();
  #endif
  pred_tot_n.allocate(styr,endyr,"pred_tot_n");
  pred_retained_size_comp.allocate(styr,endyr,1,size_n,"pred_retained_size_comp");
  #ifndef NO_AD_INITIALIZE
    pred_retained_size_comp.initialize();
  #endif
  pred_tot_size_comp.allocate(styr,endyr,1,size_n,"pred_tot_size_comp");
  #ifndef NO_AD_INITIALIZE
    pred_tot_size_comp.initialize();
  #endif
  f_mort.allocate(styr,endyr,"f_mort");
  #ifndef NO_AD_INITIALIZE
    f_mort.initialize();
  #endif
  in_molt_prob.allocate(1,size_n,"in_molt_prob");
  #ifndef NO_AD_INITIALIZE
    in_molt_prob.initialize();
  #endif
  est_molt_prob.allocate(1,size_n,"est_molt_prob");
  #ifndef NO_AD_INITIALIZE
    est_molt_prob.initialize();
  #endif
  temp_n.allocate(1,size_n,"temp_n");
  #ifndef NO_AD_INITIALIZE
    temp_n.initialize();
  #endif
  trans_n.allocate(1,size_n,"trans_n");
  #ifndef NO_AD_INITIALIZE
    trans_n.initialize();
  #endif
  temp_catch_n.allocate(1,size_n,"temp_catch_n");
  #ifndef NO_AD_INITIALIZE
    temp_catch_n.initialize();
  #endif
  sum_numbers_obs.allocate(styr,endyr,"sum_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_numbers_obs.initialize();
  #endif
  numbers_pred.allocate(styr,endyr,"numbers_pred");
  sum_ret_numbers_obs.allocate(styr,endyr,"sum_ret_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_ret_numbers_obs.initialize();
  #endif
  sum_tot_numbers_obs.allocate(styr,endyr,"sum_tot_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_tot_numbers_obs.initialize();
  #endif
  num_like.allocate("num_like");
  #ifndef NO_AD_INITIALIZE
  num_like.initialize();
  #endif
  ret_cat_like.allocate("ret_cat_like");
  #ifndef NO_AD_INITIALIZE
  ret_cat_like.initialize();
  #endif
  tot_cat_like.allocate("tot_cat_like");
  #ifndef NO_AD_INITIALIZE
  tot_cat_like.initialize();
  #endif
  surv_sc_like.allocate("surv_sc_like");
  #ifndef NO_AD_INITIALIZE
  surv_sc_like.initialize();
  #endif
  ret_comp_like.allocate("ret_comp_like");
  #ifndef NO_AD_INITIALIZE
  ret_comp_like.initialize();
  #endif
  tot_comp_like.allocate("tot_comp_like");
  #ifndef NO_AD_INITIALIZE
  tot_comp_like.initialize();
  #endif
  nat_m_like.allocate("nat_m_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_like.initialize();
  #endif
  nat_m_mu_like.allocate("nat_m_mu_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_mu_like.initialize();
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
  smooth_init_like.allocate("smooth_init_like");
  #ifndef NO_AD_INITIALIZE
  smooth_init_like.initialize();
  #endif
  f_prior.allocate("f_prior");
  #ifndef NO_AD_INITIALIZE
  f_prior.initialize();
  #endif
  total_population_n.allocate(styr,endyr,"total_population_n");
  fished_population_n.allocate(styr,endyr,"fished_population_n");
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
   n_size_pred(styr,size) = exp(log_n_init(size));
   total_fish_sel(size) = 1 / (1+exp(-fish_tot_sel_slope*(sizes(size)-fish_tot_sel_50))) ; 
   retain_fish_sel(size) = 1 / (1+exp(-fish_ret_sel_slope*(sizes(size)-fish_ret_sel_50)));
   surv_sel(size) = 1 / (1+exp(-surv_sel_slope*(sizes(size)-surv_sel_50)));
   est_molt_prob(size) = 1 - (1 / (1+exp(-molt_sel_slope*(sizes(size)-molt_sel_50))));
   }
 // specify molting probability
    in_molt_prob = prob_molt;
 // estimate molting probability
 if(est_molt==1)
	in_molt_prob = est_molt_prob;
 for(int year=styr;year<=endyr;year++)
  for(int size=1;size<=size_n;size++)
  {
  nat_m(year,size) = log_m_mu + nat_m_dev(year);
  selectivity(year,size) = surv_sel(size) + q_dev(year);
  }
  temp_prop_rec.initialize();
  tot_prop_rec = 20;
  for(int i=1;i<=2;i++)
   tot_prop_rec += prop_rec(i);
   temp_prop_rec(1) = 20/ tot_prop_rec;
  for(int i = 1;i<=2;i++)
   temp_prop_rec(i+1) = prop_rec(i)/tot_prop_rec;
 //==make last year rec dev and f dev equal to the average for now
 rec_devs(endyr) =0;
 nat_m_dev(endyr)=0;
 f_mort.initialize();
  for(int year=1;year<=ret_cat_yr_n;year++)
    f_mort(ret_cat_yrs(year)) = exp(log_f+f_dev(year));
 for(int year=styr;year<endyr;year++)
  {
  for (int size=1;size<=size_n;size++) 
	temp_n(size) = n_size_pred(year,size) * exp(-1*(0.17)*exp( nat_m(year,size)));
	  // fishery
	   for(int size=1;size<=size_n;size++)
	   {	   
	   temp_catch_n(size) = temp_n(size) * (1 -exp(-(f_mort(year)*total_fish_sel(size))));
	   pred_retained_size_comp(year,size) = temp_catch_n(size)*retain_fish_sel(size) ; 
	   pred_tot_size_comp(year,size) = temp_catch_n(size);
	   temp_n(size) = temp_n(size) * exp(-((f_mort(year)*total_fish_sel(size))));
       temp_n(size) += temp_catch_n(size)*(1-retain_fish_sel(size))*(discard_survival)	;
	   }	
	  // growth
	   trans_n = size_trans * elem_prod(temp_n,in_molt_prob);
	   temp_n = trans_n + elem_prod(temp_n,1-in_molt_prob);
	   // recruitment
       temp_n(1) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(1);
       temp_n(2) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(2);
	   temp_n(3) += exp(log_avg_rec + rec_devs(year))*temp_prop_rec(3);
      // natural mortality		
	  for (int size=1;size<=size_n;size++) 
	   n_size_pred(year+1,size) = temp_n(size) * exp(-1*(0.83)*exp(nat_m(year,size)));
    }
   evaluate_the_objective_function();
}

void model_parameters::evaluate_the_objective_function(void)
{
  // make total numbers by maturity state from obs and preds
  numbers_pred.initialize();
  sum_numbers_obs.initialize();
  pred_retained_n.initialize();
  pred_tot_n.initialize();
  total_population_n.initialize();
  fished_population_n.initialize(); 
  for (int year=styr;year<=endyr;year++)
   for (int size=1;size<=size_n;size++)
   {
    total_population_n(year)+= n_size_pred(year,size);
    numbers_pred(year)    += selectivity(year,size)*n_size_pred(year,size);
	sum_numbers_obs(year) += n_size_obs(year,size);
	pred_retained_n(year) += pred_retained_size_comp(year,size);
	pred_tot_n(year)      += pred_tot_size_comp(year,size);
	fished_population_n(year)+=n_size_pred(year,size)*retain_fish_sel(size);
   }
  // likelihoods
  num_like = 0;
  for (int year=styr;year<=endyr;year++)
   if (year!=2020)
    num_like += square( log(numbers_pred(year)) - log(n_obs(year))) / (2.0 * square(sigma_numbers(year)));
   // likelihoods
  ret_cat_like = 0;
  for (int year=1;year<=ret_cat_yr_n;year++)
    ret_cat_like += square( log(pred_retained_n(ret_cat_yrs(year))) - log(ret_cat_numbers(year))) / (2.0 * square(sigma_numbers_ret));
  tot_cat_like = 0;
  for (int year=1;year<=tot_cat_yr_n;year++)
    tot_cat_like += square( log(pred_tot_n(tot_cat_yrs(year))) - log(tot_cat_numbers(year))) / (2.0 * square(sigma_numbers_tot));
  // immature numbers at size data
  surv_sc_like = 0;
  for (int year=styr;year<=endyr;year++)
   for (int size=1;size<=size_n;size++)
    if (n_size_obs(year,size) >0.01 & n_size_pred(year,size) >0.01)
     surv_sc_like += sc_eff_samp*(n_size_obs(year,size)/sum_numbers_obs(year)) * log( (selectivity(year,size)*n_size_pred(year,size)/numbers_pred(year)) / (n_size_obs(year,size)/sum_numbers_obs(year)));
  surv_sc_like = -1*surv_sc_like;
  // retained catch at size data
  ret_comp_like = 0;
  for (int year=1;year<=ret_sc_yr_n;year++)
   for (int size=1;size<=size_n;size++)
    if (ret_cat_size(year,size) >0.01 & pred_retained_size_comp(ret_sc_yrs(year),size) >0.01)
     ret_comp_like += ret_eff_samp*(ret_cat_size(year,size)) * log( (pred_retained_size_comp(ret_sc_yrs(year),size)/pred_retained_n(ret_sc_yrs(year))) / (ret_cat_size(year,size)));
  ret_comp_like = -1*ret_comp_like;
  // total catch at size data
  tot_comp_like = 0;
  for (int year=1;year<=tot_sc_yr_n;year++)
   for (int size=1;size<=size_n;size++)
    if (tot_cat_size(year,size) >0.001 & pred_tot_size_comp(tot_sc_yrs(year),size) >0.001)
     tot_comp_like += tot_eff_samp*(tot_cat_size(year,size)) * log( (pred_tot_size_comp(tot_sc_yrs(year),size)/pred_tot_n(tot_sc_yrs(year))) / (tot_cat_size(year,size)));
  tot_comp_like = -1*tot_comp_like;
 //penalties on m 
  nat_m_mu_like =0;
  nat_m_mu_like += pow(((log_m_mu)-(log_mu_m_prior))/ (sqrt(2)*sqrt(sigma_m_mu)),2.0);
  // f_prior.initialize();;
  // for (int year=styr;year<=endyr;year++)
   // f_prior += pow((exp(log_f+f_dev(year))-1)/ (sqrt(2)*sqrt(1)),2.0);
       // cout<<4<<endl;
  if(est_m_devs>0)
  {
  nat_m_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_like += pow((nat_m(year,1)-log_m_mu)/ (sqrt(2)*sqrt(sigma_m)),2.0);
  }
  if(est_q_devs>0)
  {
  q_like =0;
  for (int year=styr;year<=endyr;year++)
   q_like += pow((selectivity(year,4)-surv_sel(4))/ (sqrt(2)*sqrt(sigma_q(1))),2.0);
   }
  smooth_q_like = 0;
  smooth_q_like = smooth_q_weight* (norm2(first_difference(first_difference(q_dev)))) ;
  smooth_m_like = 0;
  smooth_m_like = smooth_m_weight* norm2(first_difference(first_difference(nat_m_dev))) ;
  smooth_f_like = 0;
  smooth_f_like = smooth_f_weight* (norm2(first_difference(first_difference(f_dev))));
  smooth_surv_like = 0;
  smooth_surv_like = smooth_surv_weight* (norm2(first_difference(first_difference(surv_sel))));
  smooth_init_like = 0;
  smooth_init_like = smooth_init_weight* (norm2(first_difference(first_difference(log_n_init))));
  f = num_like +  ret_cat_like + tot_cat_like + surv_sc_like + ret_comp_like + tot_comp_like + 
  nat_m_like +  nat_m_mu_like +  smooth_q_like + smooth_m_like + q_like +  smooth_f_like +
  surv_sel_prior + smooth_surv_like + f_prior + smooth_init_like;
  cout<<num_like<< " "  << surv_sc_like << " " <<endl;
  cout<<ret_cat_like<< " " << ret_cat_like << " " << ret_comp_like << " " << tot_comp_like << " " <<endl;
  cout<<nat_m_like<< " " << tot_cat_like << " " << nat_m_mu_like << " "  <<endl;
  cout<<smooth_q_like<< " " << smooth_m_like << " " << q_like << " "  <<smooth_f_like<<" " << surv_sel_prior<<" "<<f_prior<<" "<<smooth_init_like<<" "<<endl;
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
  report<<f<<" "<<num_like<<" "<<surv_sc_like<<" "<<nat_m_like<<" "<< nat_m_mu_like<<" "<<smooth_q_like<<" "<<smooth_m_like<<" "<<q_like<<" "<<endl;
  report <<"$natural mortality" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << mfexp(nat_m(i))<<endl;
  }
  report <<"$surv_n_cv" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << sigma_numbers(i)<<endl;
  }
  report <<"$size_trans" << endl;
  for(int i=1; i<=size_n; i++)
  {
    report << size_trans(i)<<endl;
  }
  report <<"$recruits" << endl;
  for(int i=styr; i<=endyr; i++)
   report << mfexp(log_avg_rec + rec_devs(i))<<endl;
  report <<"$pred numbers at size" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (elem_prod(selectivity(i),n_size_pred(i)))/numbers_pred(i)<<endl;
  }
  report<<"$styr"<<endl;
  report<<styr<<endl;
  report<<"$endyr"<<endl;
  report<<endyr<<endl;
  report<<"$numbers_pred"<<endl;
  report<<numbers_pred<<endl;
  report<<"$n_obs"<<endl;
  report<<n_obs<<endl;
  report <<"$n_obs_len" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (n_size_obs(i))<<endl;
  }
  report<<"$ret_fish_sel"<<endl;
  report<<retain_fish_sel<<endl;
  report<<"$total_fish_sel"<<endl;
  report<<total_fish_sel<<endl;  
  report <<"$est_fishing_mort" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << f_mort (i) <<endl;
  }  
  report <<"$survey selectivity" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << selectivity(i)<<endl;
  }
  report <<"$pred_pop_num" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (n_size_pred(i))<<endl;
  }
  report <<"$pred_retained_n" << endl;
  for(int i=styr; i<endyr; i++)
  {
    report << (pred_retained_n(i))<<endl;
  }
  report <<"$pred_tot_n" << endl;
  for(int i=styr; i<endyr; i++)
  {
    report << (pred_tot_n(i))<<endl;
  }
  report <<"$ret_cat_numbers" << endl;
  for(int i=1; i<=ret_cat_yr_n; i++)
  {
    report << (ret_cat_numbers(i))<<endl;
  }
  report <<"$tot_cat_numbers" << endl;
  for(int i=1; i<=tot_cat_yr_n; i++)
  {
    report << (tot_cat_numbers(i))<<endl;
  }
  report <<"$pred_retained_size_comp" << endl;
  for(int i=styr; i<endyr; i++)
  {
    report << ((pred_retained_size_comp(i)/pred_retained_n(i)))<<endl;
  }
  report <<"$pred_tot_size_comp" << endl;
  for(int i=styr; i<endyr; i++)
  {
    report << ((pred_tot_size_comp(i)/pred_tot_n(i)))<<endl;
  }
  report <<"$obs_retained_size_comp" << endl;
  for(int i=1; i<=ret_sc_yr_n; i++)
  {
    report << (ret_cat_size(i))<<endl;
  }
  report <<"$obs_tot_size_comp" << endl;
  for(int i=1; i<=tot_sc_yr_n; i++)
  {
    report << (tot_cat_size(i))<<endl;
  }
  report <<"$temp_prop_rec" << endl;
  report << temp_prop_rec << endl;
  report <<"$prob_molt" << endl;
  report << (prob_molt)<<endl;
  report <<"$in_prob_molt" << endl;
  report << (in_molt_prob)<<endl;
  report <<"$sigma_numbers" << endl;
  report << (sigma_numbers)<<endl;
  report <<"$sizes" << endl;
  report << (sizes)<<endl;
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
