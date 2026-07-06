DATA_SECTION
 //Bering Sea snow crab model
  
  init_int styr  																							
  init_int endyr 
  init_int dat_yr
  init_ivector years(1,dat_yr)
  init_int size_n
  init_vector sizes(1,size_n)
  init_vector imm_n_obs(styr,endyr)
  init_matrix imm_n_size_obs(styr,endyr,1,size_n)
  init_vector mat_n_obs(styr,endyr)
  init_matrix mat_n_size_obs(styr,endyr,1,size_n)
  init_matrix prop_term_molt(styr,endyr,1,size_n)
  init_matrix size_trans(1,size_n,1,size_n) 
  init_vector survey_sel(1,size_n)
   
  init_vector sigma_numbers_imm(styr,endyr)
  init_vector sigma_numbers_mat(styr,endyr)

  init_number mat_eff_samp
  init_number imm_eff_samp

  init_vector log_mu_m_prior(1,3)
  init_number est_m_devs
  init_number est_q_devs
  init_number est_m_mat_devs
  init_number est_q_mat_devs
  init_number est_sigma_m
  init_vector sigma_m_mu(1,3)
  init_number smooth_q_weight
  init_vector smooth_m_weight(1,2)
  init_number est_log_m_mu
  init_number est_sigma_q
  init_number est_m_lg_devs
  init_number smooth_f_weight
  init_number surv_sel_cv
  init_number surv_sel_cv_2
  init_number smooth_surv_weight
  init_number est_sel
  init_number large_cutoff
  init_number est_tv_fish_sel
 
  !!cout<<"imm_n_obs"<<imm_n_obs<<endl; 
  !!cout<<"mat_n_obs"<<mat_n_obs<<endl; 
  !!cout<<"survey_sel"<<survey_sel<<endl; 
  !!cout<<"est_m_devs"<<est_m_devs<<endl;
  !!cout<<"sigma_m_mu"<<sigma_m_mu<<endl;
  !!cout<<"smooth_m_weight"<<smooth_m_weight<<endl;
  !!cout<<"est_m_lg_devs"<<est_m_lg_devs<<endl;

 //==read in catch data
 !! ad_comm::change_datafile_name("catch_dat.DAT");
  init_int ret_cat_yr_n
  init_vector ret_cat_yrs(1,ret_cat_yr_n)
  
  init_vector ret_cat_numbers(1,ret_cat_yr_n)
  init_vector disc_cat_numbers(1,ret_cat_yr_n) 
  
  init_number sigma_numbers_ret
  init_number sigma_numbers_disc
  init_number discard_survival
  init_number ret_eff_samp
  init_number disc_eff_samp
  
  init_matrix ret_cat_size(1,ret_cat_yr_n,1,size_n)
  init_matrix disc_cat_size(1,ret_cat_yr_n,1,size_n)
  
PARAMETER_SECTION
  init_bounded_vector log_n_imm(1,size_n,1.01,30,1)
  init_bounded_vector log_n_mat(1,size_n,1.01,30,1)
  init_bounded_dev_vector nat_m_dev(styr,endyr,-4,4,est_m_devs)
  init_bounded_dev_vector nat_m_mat_dev(styr,endyr,-4,4,est_m_mat_devs)
  init_bounded_dev_vector nat_m_lg_dev(styr,endyr,-4,4,est_m_mat_devs)
  init_bounded_vector q_dev(styr,endyr,-0.2,0.2,est_q_devs)
  init_bounded_vector q_mat_dev(styr,endyr,-0.2,0.2,est_q_mat_devs)
  init_bounded_vector q_lg_dev(styr,endyr,-0.2,0.2,est_q_mat_devs)
  init_bounded_number log_avg_rec(1,40)
  init_bounded_dev_vector rec_devs(styr,endyr,-10,10,1)
  init_bounded_vector sigma_m(1,3,0.01,4,est_sigma_m)
  init_bounded_vector log_m_mu(1,3,-5,3,est_log_m_mu)
  init_bounded_vector prop_rec(1,2,0.00001,200)
  init_bounded_vector sigma_q(1,2,0.01,4,est_sigma_q)
  
  init_bounded_number log_f(-5,5)
  init_bounded_dev_vector f_dev(1,ret_cat_yr_n,-5,5)
  init_bounded_number fish_ret_sel_50(25,150)
  init_bounded_number fish_ret_sel_50_post(25,150)
  init_bounded_number fish_ret_sel_slope(0.0001,20)
  init_bounded_number fish_ret_sel_slope_post(0.0001,20)
  init_bounded_number fish_tot_sel_offset(-10,20)
  init_bounded_dev_vector fish_tot_sel_offset_dev(1,ret_cat_yr_n,-40,40)
  
  init_bounded_number fish_tot_sel_slope(-10,0)
  //init_bounded_dev_vector fish_tot_sel_slope_dev(1,ret_cat_yr_n,-10,10)
  
  init_bounded_number surv_omega(0,1,est_sel)
  init_bounded_number surv_alpha1(0.0001,20,est_sel)
  init_bounded_number surv_beta1(25,150,est_sel)
  init_bounded_number surv_alpha2(0.0001,20,est_sel)
  init_bounded_number surv_beta2(25,150,est_sel)
  
  matrix imm_n_size_pred(styr,endyr,1,size_n)
  matrix mat_n_size_pred(styr,endyr,1,size_n)
  matrix nat_m(styr,endyr,1,size_n)
  matrix nat_m_mat(styr,endyr,1,size_n)
  matrix selectivity(styr,endyr,1,size_n)
  matrix selectivity_mat(styr,endyr,1,size_n)
  
  matrix total_fish_sel(styr,endyr,1,size_n)
  matrix retain_fish_sel(styr,endyr,1,size_n)
  vector pred_retained_n(styr,endyr)
  vector pred_discard_n(styr,endyr)
  matrix pred_retained_size_comp(styr,endyr,1,size_n)
  matrix pred_discard_size_comp(styr,endyr,1,size_n)
  
  vector f_mort(styr,endyr)

  vector temp_imm(1,size_n)
  vector temp_mat(1,size_n)
  vector trans_imm(1,size_n)
  vector temp_catch_imm(1,size_n)
  vector temp_catch_mat(1,size_n)

  vector surv_sel(1,size_n)

  vector sum_imm_numbers_obs(styr,endyr)
  vector sum_mat_numbers_obs(styr,endyr)
  vector imm_numbers_pred(styr,endyr)
  vector mat_numbers_pred(styr,endyr)
  vector sum_ret_numbers_obs(styr,endyr)
  vector sum_disc_numbers_obs(styr,endyr) 
  
  sdreport_vector total_population_n(styr,endyr)
  sdreport_vector fished_population_n(styr,endyr) 
   
  number imm_num_like
  number mat_num_like
  number ret_cat_like
  number disc_cat_like
  number imm_like
  number mat_like
  number ret_comp_like
  number disc_comp_like
  
  matrix use_term_molt(styr,endyr,1,size_n)
  
  number nat_m_like
  number nat_m_mat_like
  number nat_m_mu_like
  number nat_m_lg_like
  number nat_m_mat_mu_like
  number nat_m_lg_mu_like
  number smooth_q_like
  number smooth_m_like
  number q_like
  number q_mat_like
  number smooth_f_like
  number surv_sel_prior
  number smooth_surv_like
  number f_prior
  
  vector temp_prop_rec(1,3)
  number tot_prop_rec
  
  objective_function_value f
  
 
 
//==============================================================================
PROCEDURE_SECTION
 dvariable fpen=0.0;
// initial year numbers at size and selectivity


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
   
//==============================================================================
FUNCTION evaluate_the_objective_function

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
  
// ========================y==================================================   
REPORT_SECTION
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
  
RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 10000
  convergence_criteria 1e-3

