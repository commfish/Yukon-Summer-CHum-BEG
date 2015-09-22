//==============================================================================
//  Reconst_SS_SR.TPL 
//  This model is Yukon Summer Run reconstruciton Satate-Space model 
//  Written by: Toshihide "Hamachan" Hamazaki
//  Date: 06/18/2015
//==============================================================================

//=============================================================================
// 1.0  Data Entry 
//==============================================================================
DATA_SECTION
  init_int fyear;  	// First year
  init_int lyear;  	// Last year
  init_int tyear;   // Transition year
  init_int fage;   	// First age 
  init_int lage;   	// Last age 
  init_int nweir;  	// number of Weir projects 
  init_int naerial; // number of Aerial projects
  init_matrix obs_H(1,4,fyear,lyear);    // Harvests from 1 to 4 
  
// Read Pilot Sonar data: This is the anchora data   
  init_vector obs_plt(fyear,lyear);       // Pilot Sonar Passage estimate 
  init_vector obs_plt_sd(fyear,lyear);    // SD of Pilot Sonar Passage estimates
  init_vector obs_w_andr(fyear,lyear);    // Observations Andreafsky Weir 
  init_vector obs_a_andr_e(fyear,lyear);  // Observations Andreafsky Aerial East
  init_vector obs_a_andr_w(fyear,lyear);  // Observations Andreafsky Aerial West
  init_vector obs_anvk(fyear,lyear);      // Observations Anvik Sonar 1973-2002
  
// Read weir data   
  init_matrix obs_wS(1,nweir,fyear,lyear);

// Read Aerial data  
  init_matrix obs_aS(1,naerial,fyear,lyear);
  
// Read Age Compostion data  
  init_vector efN(fyear,lyear);		// effective sample size 
  init_matrix obs_age_p(fage,lage,fyear,lyear);	//Annual age composition
  init_number SDRec;			//Weight Recruitment
  init_number SDma;  			//Weight age composition 

// Read Control data    
  !! ad_comm::change_datafile_name("proj.ctl");
  init_vector lnalpha_lup(1,3);   // Ricker ln alpha
  init_vector s_beta_lup(1,3);    // Ricker beta X 100000
  init_vector phi_lup(1,3);       // AR1 phi
  init_vector ln_mu_R0_lup(1,3);  // Mean Recruitment 
  init_vector ln_Rdevs_re_lup(1,3); // Recruitment variation
  init_vector resid0_lup(1,3);      // First year residual
  init_vector ln_mu_ma_lup(1,3);			// Mean maturity years 
  init_vector ln_madev_re_lup(1,3); // deviation of maturity years   
  init_vector ln_wS_lup(1,3);		// ln transformed slope for upriver weir/tower model
  init_vector ln_aS_lup(1,3);		// ln transformed slope for upriver aerial model
  init_vector ln_rwS_lup(1,3);	// ln transformed sd for weir model
  init_vector ln_raS_lup(1,3);	// ln transformed sd for aerial model  
  init_vector p_east_lup(1,3);	// proportion of Andreafsky east fork
  init_vector ln_wandr_lup(1,3);        // ln transformed slope for Andreafsky River: 
  init_vector ln_aandr_lup(1,3);       // ln transformed slope for Andreafsky Aerial
  init_vector ln_anvk_lup(1,3);       // ln transformed slope for Anvik Sonar: 1973-2002, 2003-2013
  init_vector ln_ranvk_lup(1,3);       // ln transformed sd for Anvik Sonar
 
  !! cout << "Data Section Completed" << endl;
//!!cout<<SDma<<endl;
//!!exit(9);

//==========================================================================
INITIALIZATION_SECTION
//==========================================================================
  ln_wS   3.0;        // ln transformed slope for upriver Weir/Tower model 
  ln_aS  5.0;         // ln transformed slope for upriver aerial model 
  ln_rwS   0.0;        //  for Kwethluk weir model
  ln_raS   0.0;       // log transformed slope for Kwethluk aerial model
  ln_wandr  4.0;       // slope for Anndreafsky
  ln_aandr  5.0;       // slope for Anndreafsky East Aerial
  ln_anvk  1.0;       // slope for Anvik Sonar 1973-2002
  ln_alpha   1.2;     // Ricker ln Alpha 
  s_beta    0.7;       // Ricker beta
  phi        0.0;
  ln_mu_ma   1.2; 
  ln_mu_R0   15.0;
  p_east     0.5;
  q          1.0;
//  resid0    -0.15;
	
//==============================================================================
// 2.0  Define Parameters 
//==============================================================================
PARAMETER_SECTION
// State-Space Model parameters
  init_bounded_number ln_alpha(lnalpha_lup);			// Ricker ln alpha
  init_bounded_number s_beta(s_beta_lup);				// Ricker beta X 10000
  init_bounded_number phi(phi_lup);				// AR1 phi 
  init_bounded_number ln_mu_R0(ln_mu_R0_lup);			// Mean Recruitment   
  init_bounded_dev_vector ln_Rdevs_re(fyear-lage,lyear-fage,-20.0,20.0,2); // Recruitment variation
  init_bounded_number resid0(-10.0,10.0,1);			// First year residual
  init_bounded_number ln_mu_ma(0.6,2.0,1);			// Mean maturity years 
  init_bounded_dev_vector ln_madev_re(fyear-lage,lyear-fage,-1.0,1.0,2); // deviation of maturity years   

// Run reconstruction parameters 
  init_bounded_vector ln_wS(1,nweir,2.0,7.0,1);		// ln transformed slope for upriver weir/tower model
  init_bounded_vector ln_aS(1,naerial,1.0,8.0,1);		// ln transformed slope for upriver aerial model
  init_bounded_number ln_rwS(-3.0,10.0,3);	// ln transformed sd for weir model
  init_bounded_number ln_raS(-3.0,5.0,3);	// ln transformed sd for aerial model  
  init_bounded_number p_east(0.0,0.6,1);	// proportion of Andreafsky east fork
  init_bounded_vector ln_wandr(1,2,0.0,5.0,1);        // ln transformed slope for Andreafsky River: 1983-2002 
  init_bounded_number ln_aandr(1.0,6.0,1);       // ln transformed slope for Andreafsky Aerial
  init_bounded_vector ln_anvk(1,2,0.0,3.0,1);       // ln transformed slope for Anvik Sonar: 1973-2002, 2003-2013
  init_bounded_number ln_ranvk(-3.0,3.0,1);       // ln transformed sd for Anvik Sonar
  init_bounded_number q(0.1,2.0,1);       // Pilot Station during previous. 
 
//  init_bounded_vector ln_rwS(1,nweir,-3.0,10.0,2);	// ln transformed sd for weir model
//  init_bounded_vector ln_raS(1,naerial,-3.0,5.0,2);	// ln transformed sd for aerial model  



//  Working parameters: Run reconstructin 
  vector wS(1,nweir);			// slope for weir model
  vector aS(1,naerial);       // slope for aerial model
  number rwS;				// sd for weir model
  number raS;				// sd for aerial model
//  vector rwS(1,nweir);		// sd for weir model
//  vector raS(1,naerial);		
  vector wandr(1,2);        // log transformed slope for Andreafsky River: 1983-2002 
  number aandr;       // log transformed slope for Andreafsky Aerial
  vector anvk(1,2);       // log transformed slope for Anvik Sonar: 1973-2002, 2003-2013
  number ranvk; 
  

//  Working parameters: State-Space Model    
// Ricker Spawner-Recruit Pamaeters 	
  number alpha;
  number beta;   
  // Expected ln Recruit
  vector ln_R(fyear-lage,lyear-fage);
  vector ln_predR1(fyear,lyear-fage); 
  vector ma(fyear-lage,lyear-fage);           // Annual mean maturity age 
  matrix g(fage,lage,fyear-lage,lyear-fage);  //maturity schedule logistic 
  matrix p(fage,lage,fyear-lage,lyear-fage);  //proporion of mature for eage 
  matrix N_ta(fage,lage,fyear,lyear);			// Total Run by age 
  matrix est_age_p(fage,lage,fyear,lyear);  // Expected run age proportion 
   
// Estimated paramters size wih SD 
  sdreport_vector N(fyear,lyear);			//Total Run with SD
  sdreport_vector S(fyear,lyear);			//Total Escapement with SD	
  sdreport_vector R(fyear-lage,lyear-fage); //Recuritment with SD 
  sdreport_vector resid(fyear-1,lyear-fage); //Expected residuals
  vector obs_HT(fyear,lyear);			//Total Harvest
  vector obs_H24(fyear,lyear);			//Sum of Harvests 2 to 4 
  vector obs_H34(fyear,lyear);			//Sum of Harvests 3 to 4
  
//  vector S(fyear,lyear);			//Total Escapement with SD	
//  vector R(fyear-lage,lyear-fage); //Recuritment with SD 
//  vector resid(fyear-1,lyear-fage); //Expected residuals
  vector S_up(fyear,lyear);	    //Esc above Anvik 
  vector S_adr(fyear,lyear);	// Andreafsky escapement
  vector S_anv(fyear,lyear);    // Anvik escapement	
  vector N_plt(fyear,lyear);	// Pilot Station Run

// Likelihoood 	
  number fpen;
  vector tfw(1,nweir);   // Likelihood for weir model
  vector tfa(1,naerial);  // likelihood for aerial model
  vector tfr(1,5);          // likelihood for inriver model  
  vector tf(1,3);
  objective_function_value f;
	
 !!cout<<"parameter_section Done"<<endl;

//==============================================================================
// 5.0  Calculate sum Harvest data        
//==============================================================================
PRELIMINARY_CALCS_SECTION
  int i;
  obs_HT = colsum(obs_H);                    // Harvest Sum  
  for (i=fyear;i<=lyear;i++)
   {
    obs_H24(i) = obs_H(2,i)+obs_H(3,i)+obs_H(4,i);
	obs_H34(i) = obs_H(3,i)+obs_H(4,i);   
   }
 
  cout<<"preliminary_calcs_section Done"<<endl;
  
//==============================================================================
// 5.0  Procedure       
//==============================================================================
PROCEDURE_SECTION
  fpen = 0.0;
  convert_parameters_into_rates();
//  cout <<"OK for convert_parameters..."<<endl;

  get_initial_years_recruit();
//  cout <<"OK for get_initial_year..."<<endl;
 
  get_population_dynamics();
//  cout <<"OK for get_populaiton_dynamics..."<<endl;


  evaluate_the_objective_function();
 
   if (mceval_phase())
    {
     ofstream MCMCreport("post_samp.csv", ios::app);
	   MCMCreport <<ln_alpha << "," << s_beta << "," << phi << "," <<tf(2) << ","<< N << S<< R <<endl;
	   MCMCreport.close();
//	 cout << ln_alpha << "," << s_beta << endl;
	}

  	 
//==============================================================================
// 6.0  Function:  Convert parameters into rates      
//==============================================================================
FUNCTION convert_parameters_into_rates
  int i,j;
  wS=exp(ln_wS);       // slope for weir model
  aS=exp(ln_aS);       // slope for aerial model
  rwS=exp(ln_rwS);       // SD for weir model
  raS = exp(ln_raS);       // SD for aerial model 
  wandr = exp(ln_wandr);        // slope for Andreafsky River
  aandr = exp(ln_aandr);       // slope for Andreafsky Aerial
  anvk = exp(ln_anvk);       //  slope for Anvik
  ranvk = exp(ln_ranvk);       //  slope for Anvik  
   
  alpha = exp(ln_alpha);  				// Ricker alpha 
  beta =  s_beta/1000000;    				// Ricker beta
  resid(fyear-1) = resid0;				// First year residual   
  ln_R = ln_mu_R0;                      //   
  ma = mfexp(ln_mu_ma+ln_madev_re);      // maturity index  
// calculate maturity functions   
  for (i=fyear-lage;i<=lyear-fage;i++){
    for (j=fage;j<=lage;j++){
// matrutiry probability logistic function	
     g(j,i) = 1.0/(1.0+exp(ma(i)*(lage-j)+log(1.0/0.999-1.0)));  
		}
    }
// calculate maturity functions   
  p(fage) = g(fage);
  for (i=fyear-lage;i<=lyear-fage;i++){
    for (j=fage+1;j<=lage;j++){
     p(j,i) = g(j,i) - g(j-1,i);  
		}
    }
	
 
//==============================================================================
// 7.0  Function:  get_initial_years_abundance
// 		This function estimates previous lyears of recruitment 
//      Both abundance and length composition is estimated     
//==============================================================================
FUNCTION get_initial_years_recruit
  int i;
// GET RECRUITMENTS FOR YEARS WITH NO SR LINK
  for (i=fyear-lage;i<fyear;i++){
      ln_R(i)=ln_Rdevs_re(i)+ln_mu_R0;
      R(i)=exp(ln_R(i));
  } 
  
//==============================================================================
// 6.0  Function:  generate_maturity_schedule     
//==============================================================================
FUNCTION get_population_dynamics
  int i,j;
// CALCULATE NUMBERS AT AGE
  for (i=fyear;i<=lyear;i++) 
  {
    for (j=fage;j<=lage;j++) 
		{
// Assign each year's run by age 
           N_ta(j,i) = p(j,i-j)*R(i-j);
		}
// N is expected run size 
    N(i) = colsum(N_ta)(i);
// S is expected Escapement: Expected run size - harvest 

// est_age_p is expected age proportion 
 for (j=fage;j<=lage;j++) 
		{
// Assign each year's run by age 
           est_age_p(j,i) = N_ta(j,i)/N(i);
		}		
    S(i) = N(i) - obs_HT(i);
	S(i) = posfun(S(i), 350000,fpen);
	
// Calculate Escapement of Andreafsky, Anvik, Pilot Run  =======================
	if(i <= tyear) 
		{
		S_adr(i) = (S(i)+obs_H24(i))/wandr(1);   // Andreafsky River Escapement 
		S_anv(i) = (S(i)+obs_H(4,i)-S_adr(i))/anvk(1);  // Anvik River Escapement 
		}
	else
		{
		S_adr(i) = (S(i)+obs_H24(i))/wandr(2);
		S_anv(i) = (S(i)+obs_H(4,i)-S_adr(i))/anvk(2);
		}	
	S_up(i) = S(i) - S_anv(i) - S_adr(i);  // Upriver Escapement 
	N_plt(i) = S(i) + obs_H34(i) - S_adr(i);  // Run size at Pilot
	
// Calculte 
  if (i<=lyear-fage)
		{
// Expected log Recruit without AR1 process         
		ln_predR1(i) = ln_alpha - beta*S(i) + log(S(i));
		resid(i) = phi*resid(i-1) + ln_Rdevs_re(i);		
// Expected log Recruit with AR1 process  
		ln_R(i) = ln_predR1(i) + resid(i);
// Expected Recruit with AR1 process  
        R(i) = exp(ln_R(i));
		}		
  }
		
		
//==============================================================================
// 10.0  Likelihood Calculation  
//==============================================================================
FUNCTION evaluate_the_objective_function
  int i,j,k;
//observation model
  f=0.0;
  tfw = 0.0;              // initialilze to 0
  tfa = 0.0;              // initialilze to 0
  tfr = 0.0;              // initialilze to 0
  tf = 0.0;
  for (i=fyear;i<=lyear;i++)
  {
//=============================================================================
// 6.6 Andreafsky Escapement likelihood
//==============================================================================
//Andreafsky East Weir
  if(obs_w_andr(i)>0) 
    {
        tfr(1) += log(rwS)+0.5*square(log(obs_w_andr(i))-log(S_adr(i)*p_east))/square(rwS); 
    } 
//Andreafsky Areal East 
  if(obs_a_andr_e(i)>0) 
	{
        tfr(2) += log(raS)+0.5*square(log(obs_a_andr_e(i))-log(S_adr(i)*p_east/aandr))/square(raS); 
	}	 
//Andreafsky Aerial west  
   if(obs_a_andr_w(i)>0) 
    {
		tfr(3) += log(raS)+0.5*square(log(obs_a_andr_w(i))-log(S_adr(i)*(1-p_east)/aandr))/square(raS); 
    }
	
//==============================================================================
// 6.4  Pilot Station Likelihood
//==============================================================================
  if(obs_plt(i)>0)
    {
	if(i < 1995)
	{
	tfr(4) += square(log(obs_plt(i))-log(q*N_plt(i)))/log(square(obs_plt_sd(i)/obs_plt(i))+1);	
    }
	else
	{	
	tfr(4) += square(log(obs_plt(i))-log(N_plt(i)))/log(square(obs_plt_sd(i)/obs_plt(i))+1);
    }	
    }      
  
//==============================================================================
// 6.3 Middle (Anvik) Escapement likelihood
//==============================================================================
//Anvik Sonar 
  if(obs_anvk(i)>0) 
   {
        tfr(5) += log(ranvk)+0.5*square(log(obs_anvk(i))-log(S_anv(i)))/square(ranvk);
   }
   
//============= Weir likelihood Calculation ====================================
  for(j=1;j<=nweir;j++)
   {
  if(obs_wS(j,i)>0) 
    {
       tfw(j) += log(rwS)+0.5*square(log(obs_wS(j,i))-log(S_up(i)/wS(j)))/square(rwS);  
	} 
   }
   
//===  Aerial survey based likelihood calculation ==============================
  for(k=1;k<=naerial;k++)
   {
  if(obs_aS(k,i)>0) 
     {
       tfa(k) += (log(raS)+0.5*square(log(obs_aS(k,i))-log(S(i)/aS(k)))/square(raS));  
	 } 
    }
 }  
  
  
//Log likelhihood for Run age proportion   
  tf(1) = -(sum(elem_prod(efN,colsum(elem_prod(obs_age_p,log(est_age_p+1.e-3))))) - sum(elem_prod(efN,colsum(elem_prod(obs_age_p,log(obs_age_p+1.e-3))))));   
//deviation in recruits  
  tf(2) = norm2(ln_Rdevs_re)/(2*SDRec*SDRec);  
//deviation in maturity scheduel  
  tf(3) = norm2(ln_madev_re)/(2*SDma*SDma);    

// Sum all likelihood ==========================================================

  f = sum(tf)+sum(tfw)+sum(tfa)+sum(tfr)+fpen;    
 

REPORT_SECTION

  report << "f" << endl << f << endl;  // Total likelihood
// Individual likelihood
  report << "tfw" << endl << tfw << endl;	// Weir 
  report << "tfa" << endl << tfa << endl;	// Aerial 
  report << "tfr" << endl << tfr << endl;   // 
  report << "tf" << endl << tf <<endl;		//State-Space 
  report << "N_ta" << endl << N_ta << endl;

  
//==============================================================================

GLOBALS_SECTION
  #include <math.h>
  #include <admodel.h>
  #include <time.h>

  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
  int header;

// ===========================================================================

TOP_OF_MAIN_SECTION
  arrmblsize = 10000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  time(&start);

// ===========================================================================
FINAL_SECTION

 // Output summary stuff
 time(&finish);
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;

 
