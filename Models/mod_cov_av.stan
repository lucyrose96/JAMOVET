// Stan model for joint analysis of asymptomatic and symptomatic cases
// in COVID-19 vaccine trial data

// Final covariates model (with bias-adjustment)

data {
  int<lower=0> N; // Sample size
  real<lower=0> N_real; // Sample size (real for division)
  // Number in each small category 
  //ref = 0 for all variables, cov1 = number only for cov1 (reference for all other variable), cov1_cov2= number only for cov1 and 2
  real<lower=0> N_ref;
  real<lower=0> N_cov1;
  real<lower=0> N_cov2;// Number in each categorical variable for calculating average VE's. Age (cov3) is already averaged because it's centred
  real<lower=0> N_cov4;
  real<lower=0> N_cov5;
  real<lower=0> N_cov4_cov5;
  real<lower=0> N_cov2_cov4;
  real<lower=0> N_cov2_cov5;
  real<lower=0> N_cov2_cov4_cov5;
  real<lower=0> N_cov1_cov5;
  real<lower=0> N_cov1_cov4;
  real<lower=0> N_cov1_cov4_cov5;
  //Number obese or not obese
  real<lower=0> N_cov4n; // N not obese
  real<lower=0> N_cov4y; // N obese
  int<lower=0> C[N]; // N symptomatic cases
  int<lower=0> A[N]; // N asymptomatic cases
  int<lower=0, upper=1> vaccinated[N]; // vaccinated indicator variable
  real <lower=0> pers_yrs_at_risk[N]; // person-years at risk
  real <lower=0, upper=1> rel_sens_asym[N]; //relative sensitivity of detecting an asymptomatic to symptomatic infection
  real <lower=0, upper=1> fp[N];  //probability of receiving a false positive over follow-up
  real cov1[N]; // HCW 0 covid patients
  real cov2[N]; // HCW 1+ covid patients
  real cov3[N]; // Age (years)
  real cov4[N]; // Obese (y/n)
  real cov5[N]; // White (y/n)

   }

parameters {
  //PS
  real <lower=-50, upper=50>alpha_ref;
  real <lower=-50, upper=50>alpha1;
  real <lower=-50, upper=50>alpha2;
  real <lower=-50, upper=50>alpha3; 
  real <lower=-50, upper=50>alpha4;
  real <lower=-50, upper=50>alpha5;
  
  //VEpr
  real <lower=-50, upper=50> beta_ref;
  real <lower=-50, upper=50> beta3;
  real <lower=-50, upper=50> beta5;
  
  //FOI
  real <lower=-50, upper=50>gam_ref;
  real <lower=-50, upper=50>gam1;
  real <lower=-50, upper=50>gam2;
  real <lower=-50, upper=50>gam3; 
  
  //VEin
  real <lower=-50, upper=50> delta_ref;
  real <lower=-50, upper=50> delta3;
 
}

transformed parameters {
  real <lower=0> mu_c[N]; // rate symptomatic infection
  real <lower=0> mu_a[N]; // rate asymptomatic infection
  real <lower=0, upper =1> ps[N]; // probability of symptoms in unvaccinated
  real <lower=0> lam[N]; // rate of infection in unvaccinated

  for (i in 1:N){
    
    ps[i] = inv_logit(alpha_ref + alpha1*cov1[i] + alpha2*cov2[i] + alpha3*cov3[i] + alpha4*cov4[i] + alpha5*cov5[i] + 
                      beta_ref*vaccinated[i] + beta3*vaccinated[i]*cov3[i] + beta5*vaccinated[i]*cov5[i] );
                  
    lam[i] = exp(gam_ref + gam1*cov1[i] + gam2*cov2[i] + gam3*cov3[i]  + 
                 delta_ref*vaccinated[i] + delta3*vaccinated[i]*cov3[i])*pers_yrs_at_risk[i];
    
    mu_c[i] = ps[i]*lam[i];
    mu_a[i] = (1-ps[i])*lam[i]*rel_sens_asym[i]+fp[i];
    
  }
}

model {
  target+= poisson_lpmf(C | mu_c);
  target+= poisson_lpmf(A | mu_a);
}

generated quantities{
  
//Absolute values
 // FOI
    //Placebo
    real FOI_p_ref=exp(gam_ref); 
    real FOI_p_cov1=exp(gam_ref+gam1);
    real FOI_p_cov2=exp(gam_ref+gam2);
    real FOI_p_cov3=exp(gam_ref+gam3);//10 years older than reference
    real FOI_p_cov4=exp(gam_ref);
    real FOI_p_cov5=exp(gam_ref);
    real FOI_p_cov3_30yrs=exp(gam_ref+gam3*3);//estimate for participant 30 years older than the reference age (30 years)
    
    real FOI_p_cov2_cov4=exp(gam_ref+gam2);
    real FOI_p_cov1_cov4=exp(gam_ref+gam1);
    real FOI_p_cov4_cov5=exp(gam_ref);
    real FOI_p_cov2_cov5=exp(gam_ref+gam2);
    real FOI_p_cov2_cov4_cov5=exp(gam_ref+gam2);
    real FOI_p_cov1_cov5=exp(gam_ref+gam1);
    real FOI_p_cov1_cov4_cov5=exp(gam_ref+gam1);
    
    //Vaccine
    real FOI_v_ref=exp(gam_ref+delta_ref); 
    real FOI_v_cov1=exp(gam_ref+gam1 + delta_ref);
    real FOI_v_cov2=exp(gam_ref+gam2 + delta_ref);
    real FOI_v_cov3=exp(gam_ref+gam3 + delta_ref+delta3);
    real FOI_v_cov4=exp(gam_ref+delta_ref);
    real FOI_v_cov5=exp(gam_ref+delta_ref);
    real FOI_v_cov3_30yrs=exp(gam_ref+gam3*3 + delta_ref+delta3*3);
    
  // PS
    //Placebo
    real PS_p_ref=inv_logit(alpha_ref);
    real PS_p_cov1=inv_logit(alpha_ref+alpha1);
    real PS_p_cov2=inv_logit(alpha_ref+alpha2);
    real PS_p_cov3=inv_logit(alpha_ref+alpha3);
    real PS_p_cov4=inv_logit(alpha_ref+alpha4);
    real PS_p_cov5=inv_logit(alpha_ref+alpha5);
    real PS_p_cov3_30yrs=inv_logit(alpha_ref+alpha3*3);
    
    real PS_p_cov2_cov4=inv_logit(alpha_ref+alpha2+alpha4);
    real PS_p_cov1_cov4=inv_logit(alpha_ref+alpha1+alpha4);
    real PS_p_cov4_cov5=inv_logit(alpha_ref+alpha4+alpha5);
    real PS_p_cov2_cov5=inv_logit(alpha_ref+alpha2+alpha5);
    real PS_p_cov2_cov4_cov5=inv_logit(alpha_ref+alpha2+alpha4+alpha5);
    real PS_p_cov1_cov5=inv_logit(alpha_ref+alpha1+alpha5);
    real PS_p_cov1_cov4_cov5=inv_logit(alpha_ref+alpha1+alpha4+alpha5);
    
    
    //Vaccine
    real PS_v_ref=inv_logit(alpha_ref + beta_ref);
    real PS_v_cov1=inv_logit(alpha_ref+alpha1 + beta_ref);
    real PS_v_cov2=inv_logit(alpha_ref+alpha2 + beta_ref);
    real PS_v_cov3=inv_logit(alpha_ref+alpha3 + beta_ref+beta3);
    real PS_v_cov4=inv_logit(alpha_ref+alpha4 + beta_ref);
    real PS_v_cov5=inv_logit(alpha_ref+alpha5 + beta_ref+beta5);
    real PS_v_cov3_30yrs=inv_logit(alpha_ref+alpha3*3 + beta_ref+beta3*3);
   
   //VEin
    real RR_in_ref=exp(delta_ref);
    real RR_in_cov1=exp(delta_ref);
    real RR_in_cov2=exp(delta_ref);
    real RR_in_cov3=exp(delta_ref+delta3);
    real RR_in_cov4=exp(delta_ref);
    real RR_in_cov5=exp(delta_ref);
    real RR_in_cov3_30yrs=exp(delta_ref+delta3*3);//don't need to do any more for averages as this is already averaged (the ref)
    
    real VE_in_ref=1-RR_in_ref;
    real VE_in_cov1=1-RR_in_cov1;
    real VE_in_cov2=1-RR_in_cov2;
    real VE_in_cov3=1-RR_in_cov3;
    real VE_in_cov4=1-RR_in_cov4;
    real VE_in_cov5=1-RR_in_cov5;
    real VE_in_cov3_30yrs=1-RR_in_cov3_30yrs;
    
  //VEpr
    real OR_pr_ref=exp(beta_ref); 
    real RR_pr_ref= OR_pr_ref/(1-PS_p_ref*(1-OR_pr_ref));
    real VE_pr_ref= 1-RR_pr_ref;
    
    real OR_pr_cov1=exp(beta_ref); 
    real RR_pr_cov1=  OR_pr_cov1/(1-PS_p_cov1*(1-OR_pr_cov1));
    real VE_pr_cov1= 1-RR_pr_cov1;
    
    real OR_pr_cov2=exp(beta_ref); 
    real RR_pr_cov2=  OR_pr_cov2/(1-PS_p_cov2*(1-OR_pr_cov2));
    real VE_pr_cov2= 1-RR_pr_cov2;
    
    real OR_pr_cov3=exp(beta_ref+beta3); 
    real RR_pr_cov3=  OR_pr_cov3/(1-PS_p_cov3*(1-OR_pr_cov3));
    real VE_pr_cov3= 1-RR_pr_cov3;
    
    real OR_pr_cov4=exp(beta_ref); 
    real RR_pr_cov4=  OR_pr_cov4/(1-PS_p_cov4*(1-OR_pr_cov4));
    real VE_pr_cov4= 1-RR_pr_cov4;
    
    real OR_pr_cov5=exp(beta_ref+beta5); 
    real RR_pr_cov5=  OR_pr_cov5/(1-PS_p_cov5*(1-OR_pr_cov5));
    real RR_pr_cov5_2=  OR_pr_cov5/(1-PS_p_cov5*(1-OR_pr_cov5));
    real VE_pr_cov5= 1-RR_pr_cov5;
    
    real OR_pr_cov3_30yrs=exp(beta_ref+beta3*3); 
    real RR_pr_cov3_30yrs=  OR_pr_cov3_30yrs/(1-PS_p_cov3_30yrs*(1-OR_pr_cov3_30yrs));
    real VE_pr_cov3_30yrs= 1-RR_pr_cov3_30yrs;
    
    real OR_pr_cov2_cov4= exp(beta_ref); 
    real RR_pr_cov2_cov4=  OR_pr_cov2_cov4/(1-PS_p_cov2_cov4*(1-OR_pr_cov2_cov4));
    real VE_pr_cov2_cov4= 1-RR_pr_cov2_cov4;
    
    real OR_pr_cov1_cov4= exp(beta_ref);
    real RR_pr_cov1_cov4=  OR_pr_cov1_cov4/(1-PS_p_cov1_cov4*(1-OR_pr_cov1_cov4));
    real VE_pr_cov1_cov4= 1-RR_pr_cov1_cov4;
    
    real OR_pr_cov4_cov5= exp(beta_ref+beta5);
    real RR_pr_cov4_cov5=  OR_pr_cov4_cov5/(1-PS_p_cov4_cov5*(1-OR_pr_cov4_cov5));
    real VE_pr_cov4_cov5= 1-RR_pr_cov4_cov5;
    
    real OR_pr_cov2_cov5= exp(beta_ref+beta5);
    real RR_pr_cov2_cov5=  OR_pr_cov2_cov5/(1-PS_p_cov2_cov5*(1-OR_pr_cov2_cov5));
    real VE_pr_cov2_cov5= 1-RR_pr_cov2_cov5;
    
    real OR_pr_cov2_cov4_cov5= exp(beta_ref+beta5);
    real RR_pr_cov2_cov4_cov5=  OR_pr_cov2_cov4_cov5/(1-PS_p_cov2_cov4_cov5*(1-OR_pr_cov2_cov4_cov5));
    real VE_pr_cov2_cov4_cov5= 1-RR_pr_cov2_cov4_cov5;
   
    real OR_pr_cov1_cov5= exp(beta_ref+beta5);
    real RR_pr_cov1_cov5=  OR_pr_cov1_cov5/(1-PS_p_cov1_cov5*(1-OR_pr_cov1_cov5));
    real VE_pr_cov1_cov5= 1-RR_pr_cov1_cov5;
   
    real OR_pr_cov1_cov4_cov5= exp(beta_ref+beta5);
    real RR_pr_cov1_cov4_cov5=  OR_pr_cov1_cov4_cov5/(1-PS_p_cov1_cov4_cov5*(1-OR_pr_cov1_cov4_cov5));
    real VE_pr_cov1_cov4_cov5= 1-RR_pr_cov1_cov4_cov5;
   
    
  //VEsym
  real VE_sym_ref=1-(1-VE_in_ref)*(1-VE_pr_ref);
  real VE_sym_cov1=1-(1-VE_in_cov1)*(1-VE_pr_cov1);
  real VE_sym_cov2=1-(1-VE_in_cov2)*(1-VE_pr_cov2);
  real VE_sym_cov3=1-(1-VE_in_cov3)*(1-VE_pr_cov3);
  real VE_sym_cov4=1-(1-VE_in_cov4)*(1-VE_pr_cov4);
  real VE_sym_cov5=1-(1-VE_in_cov5)*(1-VE_pr_cov5);
  real VE_sym_cov3_30yrs=1-(1-VE_in_cov3_30yrs)*(1-VE_pr_cov3_30yrs);
  
  real VE_sym_cov2_cov4=1-(1-VE_in_ref)*(1-VE_pr_cov2_cov4); //note VEin is averaged anyway so here just use the reference 
  real VE_sym_cov1_cov4=1-(1-VE_in_ref)*(1-VE_pr_cov1_cov4); 
  real VE_sym_cov4_cov5=1-(1-VE_in_ref)*(1-VE_pr_cov4_cov5); 
  real VE_sym_cov2_cov5=1-(1-VE_in_ref)*(1-VE_pr_cov2_cov5); 
  real VE_sym_cov2_cov4_cov5=1-(1-VE_in_ref)*(1-VE_pr_cov2_cov4_cov5); 
  real VE_sym_cov1_cov5=1-(1-VE_in_ref)*(1-VE_pr_cov1_cov5); 
  real VE_sym_cov1_cov4_cov5=1-(1-VE_in_ref)*(1-VE_pr_cov1_cov4_cov5); 
  
  //VEasym
  real VE_asym_ref=1-((1-PS_p_ref)*(1-VE_in_ref)+PS_p_ref*(1-VE_in_ref)*VE_pr_ref)/(1-PS_p_ref);
  real VE_asym_cov1=1-((1-PS_p_cov1)*(1-VE_in_cov1)+PS_p_cov1*(1-VE_in_cov1)*VE_pr_cov1)/(1-PS_p_cov1); 
  real VE_asym_cov2=1-((1-PS_p_cov2)*(1-VE_in_cov2)+PS_p_cov2*(1-VE_in_cov2)*VE_pr_cov2)/(1-PS_p_cov2); 
  real VE_asym_cov3=1-((1-PS_p_cov3)*(1-VE_in_cov3)+PS_p_cov3*(1-VE_in_cov3)*VE_pr_cov3)/(1-PS_p_cov3); 
  real VE_asym_cov4=1-((1-PS_p_cov4)*(1-VE_in_cov4)+PS_p_cov4*(1-VE_in_cov4)*VE_pr_cov4)/(1-PS_p_cov4); 
  real VE_asym_cov5=1-((1-PS_p_cov5)*(1-VE_in_cov5)+PS_p_cov5*(1-VE_in_cov5)*VE_pr_cov5)/(1-PS_p_cov5); 
  real VE_asym_cov3_30yrs=1-((1-PS_p_cov3_30yrs)*(1-VE_in_cov3_30yrs)+PS_p_cov3_30yrs*(1-VE_in_cov3_30yrs)*VE_pr_cov3_30yrs)/(1-PS_p_cov3_30yrs); 
  
  real VE_asym_cov2_cov4= 1-((1-PS_p_cov2_cov4)*(1-VE_in_ref)+PS_p_cov2_cov4*(1-VE_in_ref)*VE_pr_cov2_cov4)/(1-PS_p_cov2_cov4);  
  real VE_asym_cov1_cov4= 1-((1-PS_p_cov1_cov4)*(1-VE_in_ref)+PS_p_cov1_cov4*(1-VE_in_ref)*VE_pr_cov1_cov4)/(1-PS_p_cov1_cov4);  
  real VE_asym_cov4_cov5= 1-((1-PS_p_cov4_cov5)*(1-VE_in_ref)+PS_p_cov4_cov5*(1-VE_in_ref)*VE_pr_cov4_cov5)/(1-PS_p_cov4_cov5);  
  real VE_asym_cov2_cov5= 1-((1-PS_p_cov2_cov5)*(1-VE_in_ref)+PS_p_cov2_cov5*(1-VE_in_ref)*VE_pr_cov2_cov5)/(1-PS_p_cov2_cov5);  
  real VE_asym_cov2_cov4_cov5= 1-((1-PS_p_cov2_cov4_cov5)*(1-VE_in_ref)+PS_p_cov2_cov4_cov5*(1-VE_in_ref)*VE_pr_cov2_cov4_cov5)/(1-PS_p_cov2_cov4_cov5);  
  real VE_asym_cov1_cov5= 1-((1-PS_p_cov1_cov5)*(1-VE_in_ref)+PS_p_cov1_cov5*(1-VE_in_ref)*VE_pr_cov1_cov5)/(1-PS_p_cov1_cov5);  
  real VE_asym_cov1_cov4_cov5= 1-((1-PS_p_cov1_cov4_cov5)*(1-VE_in_ref)+PS_p_cov1_cov4_cov5*(1-VE_in_ref)*VE_pr_cov1_cov4_cov5)/(1-PS_p_cov1_cov4_cov5); 

  //Variable effects
    //FOI
  real FOI_diff_cov1=exp(gam1);
  real FOI_diff_cov2=exp(gam2);
  real FOI_diff_cov3=exp(gam3);
   //PS
  real OR_ps_diff_cov1=exp(alpha1);
  real PS_diff_cov1= OR_ps_diff_cov1/(1-PS_p_ref*(1-OR_ps_diff_cov1));
  real OR_ps_diff_cov2=exp(alpha2);
  real PS_diff_cov2= OR_ps_diff_cov2/(1-PS_p_ref*(1-OR_ps_diff_cov2));
  real OR_ps_diff_cov3=exp(alpha3);
  real PS_diff_cov3= OR_ps_diff_cov3/(1-PS_p_ref*(1-OR_ps_diff_cov3));
  real OR_ps_diff_cov4=exp(alpha4);
  real PS_diff_cov4= OR_ps_diff_cov4/(1-PS_p_ref*(1-OR_ps_diff_cov4));
  real OR_ps_diff_cov5=exp(alpha5);
  real PS_diff_cov5= OR_ps_diff_cov5/(1-PS_p_ref*(1-OR_ps_diff_cov5));
   //VEin
  real RR_in_diff_cov3=exp(delta3); //difference in relative risk of infection as a result of vaccination between age 30 and 40
  real RR_in_diff_cov3_30yrs=exp(delta3*3); //60 vs 30
    //VEpr
   
  real RR_pr_diff_cov3=(PS_v_cov3/PS_p_cov3)/(PS_v_ref/PS_p_ref);
  real RR_pr_diff_cov3_30yrs=(PS_v_cov3/PS_p_cov3)/(PS_v_ref/PS_p_ref);
  real RR_pr_diff_cov5=(PS_v_cov5/PS_p_cov5)/(PS_v_ref/PS_p_ref);
   

  
  // Averages
  
  real FOI_p= FOI_p_ref*(N_ref/N_real) + FOI_p_cov1*(N_cov1/N_real) + FOI_p_cov2*(N_cov2/N_real) + FOI_p_cov4*(N_cov4/N_real) + FOI_p_cov5*(N_cov5/N_real)+
          FOI_p_cov2_cov4*(N_cov2_cov4/N_real) + FOI_p_cov1_cov4*(N_cov1_cov4/N_real) + FOI_p_cov4_cov5*(N_cov4_cov5/N_real) + FOI_p_cov2_cov5*(N_cov2_cov5/N_real)
          + FOI_p_cov2_cov4_cov5*(N_cov2_cov4_cov5/N_real) + FOI_p_cov1_cov5*(N_cov1_cov5/N_real) + FOI_p_cov1_cov4_cov5*(N_cov1_cov4_cov5/N_real);
          
  real PS_p= PS_p_ref*(N_ref/N_real) + PS_p_cov1*(N_cov1/N_real) + PS_p_cov2*(N_cov2/N_real) + PS_p_cov4*(N_cov4/N_real) + PS_p_cov5*(N_cov5/N_real)+
          PS_p_cov2_cov4*(N_cov2_cov4/N_real) + PS_p_cov1_cov4*(N_cov1_cov4/N_real) + PS_p_cov4_cov5*(N_cov4_cov5/N_real) + PS_p_cov2_cov5*(N_cov2_cov5/N_real)
          + PS_p_cov2_cov4_cov5*(N_cov2_cov4_cov5/N_real) + PS_p_cov1_cov5*(N_cov1_cov5/N_real) + PS_p_cov1_cov4_cov5*(N_cov1_cov4_cov5/N_real);
   
  real VE_in=VE_in_ref*(N_cov4n/N_real) + VE_in_cov4*(N_cov4y/N_real); //weighted average of VEin in non-obese (ref) and obese (cov4) participants
    //age is already averaged as it's centred

  real VE_pr= VE_pr_ref*(N_ref/N_real) + VE_pr_cov1*(N_cov1/N_real) + VE_pr_cov2*(N_cov2/N_real) + VE_pr_cov4*(N_cov4/N_real) + VE_pr_cov5*(N_cov5/N_real)+
          VE_pr_cov2_cov4*(N_cov2_cov4/N_real) + VE_pr_cov1_cov4*(N_cov1_cov4/N_real) + VE_pr_cov4_cov5*(N_cov4_cov5/N_real) + VE_pr_cov2_cov5*(N_cov2_cov5/N_real)
          + VE_pr_cov2_cov4_cov5*(N_cov2_cov4_cov5/N_real) + VE_pr_cov1_cov5*(N_cov1_cov5/N_real) + VE_pr_cov1_cov4_cov5*(N_cov1_cov4_cov5/N_real);
   
  real VE_sym= VE_sym_ref*(N_ref/N_real) + VE_sym_cov1*(N_cov1/N_real) + VE_sym_cov2*(N_cov2/N_real) + VE_sym_cov4*(N_cov4/N_real) + VE_sym_cov5*(N_cov5/N_real)+
          VE_sym_cov2_cov4*(N_cov2_cov4/N_real) + VE_sym_cov1_cov4*(N_cov1_cov4/N_real) + VE_sym_cov4_cov5*(N_cov4_cov5/N_real) + VE_sym_cov2_cov5*(N_cov2_cov5/N_real)
          + VE_sym_cov2_cov4_cov5*(N_cov2_cov4_cov5/N_real) + VE_sym_cov1_cov5*(N_cov1_cov5/N_real) + VE_sym_cov1_cov4_cov5*(N_cov1_cov4_cov5/N_real);
      
  real VE_asym= VE_asym_ref*(N_ref/N_real) + VE_asym_cov1*(N_cov1/N_real) + VE_asym_cov2*(N_cov2/N_real) + VE_asym_cov4*(N_cov4/N_real) + VE_asym_cov5*(N_cov5/N_real)+
          VE_asym_cov2_cov4*(N_cov2_cov4/N_real) + VE_asym_cov1_cov4*(N_cov1_cov4/N_real) + VE_asym_cov4_cov5*(N_cov4_cov5/N_real) + VE_asym_cov2_cov5*(N_cov2_cov5/N_real)
          + VE_asym_cov2_cov4_cov5*(N_cov2_cov4_cov5/N_real) + VE_asym_cov1_cov5*(N_cov1_cov5/N_real) + VE_asym_cov1_cov4_cov5*(N_cov1_cov4_cov5/N_real);
  
 
  // log likelihood
vector[N] log_lik;
  for (i in 1:N){
  log_lik[i] = poisson_lpmf(C[i] | mu_c[i]) + poisson_lpmf(A[i] | mu_a[i]);
  }
  
}



