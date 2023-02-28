// Stan model for joint analysis of asymptomatic and symptomatic cases
// in COVID-19 vaccine trial data

// Final bias-adjusted model. 

data {
  int<lower=0> N; // Sample size
  int<lower=0> C[N]; // N symptomatic cases
  int<lower=0> A[N]; // N asymptomatic cases
  int<lower=0, upper=1> vaccinated[N]; // vaccinated indicator variable
  real <lower=0> pers_yrs_at_risk[N]; // person-years at risk
  real <lower=0, upper=1> rel_sens_asym[N]; //relative sensitivity of detecting an asymptomatic to symptomatic infection
  real <lower=0, upper=1> fp[N]; //probability of receiving a false positive over follow-up
   }

parameters {
  real alpha;
  real beta;
  real gam;
  real delta;
}

transformed parameters {
  real <lower=0> mu_c[N]; // rate symptomatic infection
  real <lower=0> mu_a[N]; // rate asymptomatic infection
  real <lower=0, upper =1> ps[N]; // probability of symptoms in unvaccinated
  real <lower=0> lam[N]; // rate of infection in unvaccinated

  for (i in 1:N){
    
    ps[i] = inv_logit(alpha + beta*vaccinated[i]);
    lam[i] = exp(gam + delta*vaccinated[i])*pers_yrs_at_risk[i];
    
    mu_c[i] = ps[i]*lam[i];
    mu_a[i] = (1-ps[i])*lam[i]*rel_sens_asym[i]+fp[i];
    
  }
}

model {
  target+= poisson_lpmf(C | mu_c);
  target+= poisson_lpmf(A | mu_a);
}

generated quantities{
  
  real FOI_p=exp(gam); //rate of infection - placebo
  real FOI_v=exp(gam+delta); // rate of infection - vaccine
  
  real RR_in=exp(delta); //risk ratio infection V vs P
  real VE_in=1-RR_in; // VEin
  
  real PS_p=inv_logit(alpha); //prob. symtoms - placebo
  real PS_v=inv_logit(alpha+beta); //prob. symtoms - vaccine
  real OR_pr=exp(beta); //odds ratio symptoms|infection vax vs unvax
  
  real RR_pr= OR_pr/(1-PS_p*(1-OR_pr)); //risk ratio symptoms|infection vax vs unvax
  real VE_pr=1-RR_pr; //VEpr
  
  real VE_sym=1-(1-VE_in)*(1-VE_pr); //VEsym
  real VE_asym=1-((1-PS_p)*(1-VE_in)+PS_p*(1-VE_in)*VE_pr)/(1-PS_p); //VEasym
  
  //log likelihood
  vector[N] log_lik;
  for (i in 1:N){
    log_lik[i] = poisson_lpmf(C[i] | mu_c[i]) + poisson_lpmf(A[i] | mu_a[i]);
  }
}


