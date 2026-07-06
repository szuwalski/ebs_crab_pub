DATA_SECTION
 //Bering Sea snow crab model
  
  init_int styr  																							
  init_int endyr 
  init_int dat_yr
  init_ivector years(1,dat_yr)
  init_int size_n
  init_vector sizes(1,size_n)
  init_vector n_obs(styr,endyr)
  init_matrix n_size_obs(styr,endyr,1,size_n)
  init_vector prob_molt(1,size_n)
  init_matrix size_trans(1,size_n,1,size_n) 
   
  init_vector sigma_numbers(styr,endyr)
  init_number sc_eff_samp

  init_number log_mu_m_prior
  init_number est_m_devs
  init_number est_q_devs
  init_number est_sigma_m
  init_number sigma_m_mu
  init_number smooth_q_weight
  init_number smooth_m_weight
  init_number est_log_m_mu
  init_number est_sigma_q
  init_number smooth_f_weight
  init_number surv_sel_cv
  init_number smooth_surv_weight
  init_number smooth_init_weight
  init_number est_molt
 
   !!cout<<"n_obs"<<n_obs<<endl; 
   !!cout<<"est_m_devs"<<est_m_devs<<endl;
   !!cout<<"sigma_m_mu"<<sigma_m_mu<<endl;
   !!cout<<"smooth_m_weight"<<smooth_m_weight<<endl;
   !!cout<<"sigma_numbers"<<sigma_numbers<<endl;

 //==read in catch data
 !! ad_comm::change_datafile_name("catch_dat.DAT");
  init_int ret_cat_yr_n
  init_vector ret_cat_yrs(1,ret_cat_yr_n)
  init_vector ret_cat_numbers(1,ret_cat_yr_n)

  init_int ret_sc_yr_n
  init_vector ret_sc_yrs(1,ret_sc_yr_n)  
  init_matrix ret_cat_size(1,ret_sc_yr_n,1,size_n)
   
  init_int tot_cat_yr_n
  init_vector tot_cat_yrs(1,tot_cat_yr_n)
  init_vector tot_cat_numbers(1,tot_cat_yr_n)

  init_int tot_sc_yr_n
  init_vector tot_sc_yrs(1,tot_sc_yr_n)  
  init_matrix tot_cat_size(1,tot_sc_yr_n,1,size_n)
  
  init_number sigma_numbers_ret
  init_number sigma_numbers_tot
  init_number discard_survival
  init_number ret_eff_samp
  init_number tot_eff_samp
  
   !!cout<<"ret_cat_numbers"<<ret_cat_numbers<<endl; 
   !!cout<<"tot_cat_numbers"<<tot_cat_numbers<<endl;
   !!cout<<"tot_sc_yrs"<<tot_sc_yrs<<endl;
   !!cout<<"sigma_numbers_ret"<<sigma_numbers_ret<<endl;
   !!cout<<"tot_eff_samp"<<tot_eff_samp<<endl;
  
PARAMETER_SECTION
  init_bounded_vector log_n_init(1,size_n,1.01,30,1)
  init_bounded_dev_vector nat_m_dev(styr,endyr,-4,4,est_m_devs)
  init_bounded_vector q_dev(styr,endyr,-0.2,0.2,est_q_devs)
  init_bounded_number log_avg_rec(1,40)
  init_bounded_dev_vector rec_devs(styr,endyr,-12,12,1)
  init_bounded_number sigma_m(0.01,4,est_sigma_m)
  init_bounded_number log_m_mu(-5,3,est_log_m_mu)
  init_bounded_vector prop_rec(1,2,0.00001,200)
  init_bounded_vector sigma_q(1,2,0.01,4,est_sigma_q)
  
  init_bounded_number log_f(-10,5)
  init_bounded_dev_vector f_dev(1,ret_cat_yr_n,-5,5)
  init_bounded_number fish_ret_sel_50(25,150,-1)
  init_bounded_number fish_ret_sel_slope(0.0001,20,-1)
  init_bounded_number fish_tot_sel_50(25,150)
  init_bounded_number fish_tot_sel_slope(0.0001,20) 
  init_bounded_number surv_sel_50(25,150)
  init_bounded_number surv_sel_slope(0.0001,20) 
  
  init_bounded_number molt_sel_50(25,180)
  init_bounded_number molt_sel_slope(0.0001,20) 
  
  matrix n_size_pred(styr,endyr,1,size_n)
  matrix nat_m(styr,endyr,1,size_n)
  matrix selectivity(styr,endyr,1,size_n)
  
  vector total_fish_sel(1,size_n)
  vector retain_fish_sel(1,size_n)
  vector surv_sel(1,size_n)
  vector pred_retained_n(styr,endyr)
  sdreport_vector pred_tot_n(styr,endyr)
  matrix pred_retained_size_comp(styr,endyr,1,size_n)
  matrix pred_tot_size_comp(styr,endyr,1,size_n)
  
  vector f_mort(styr,endyr)
  vector in_molt_prob(1,size_n)
  vector est_molt_prob(1,size_n) 
   
  vector temp_n(1,size_n)
  vector trans_n(1,size_n)
  vector temp_catch_n(1,size_n)

  vector sum_numbers_obs(styr,endyr)
  sdreport_vector numbers_pred(styr,endyr)
  vector sum_ret_numbers_obs(styr,endyr)
  vector sum_tot_numbers_obs(styr,endyr) 
  
  number num_like
  number ret_cat_like
  number tot_cat_like
  number surv_sc_like
  number ret_comp_like
  number tot_comp_like
   
  number nat_m_like
  number nat_m_mu_like
  number smooth_q_like
  number smooth_m_like
  number q_like
  number q_mat_like
  number smooth_f_like
  number surv_sel_prior
  number smooth_surv_like
    number smooth_init_like
  number f_prior
  
    sdreport_vector total_population_n(styr,endyr)
    sdreport_vector fished_population_n(styr,endyr)	
  vector temp_prop_rec(1,3)
  number tot_prop_rec
  
  objective_function_value f
  
 
 
//==============================================================================
PROCEDURE_SECTION
 dvariable fpen=0.0;
// initial year numbers at size and selectivity
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
   
//==============================================================================
FUNCTION evaluate_the_objective_function

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
  
// ========================y==================================================   
REPORT_SECTION
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
  
RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 10000
  convergence_criteria 1e-3

