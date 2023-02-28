// Stan model for  analysis of any infection
// in COVID-19 vaccine trial data
//Nicholas Grassly 25 Aug 2021

data {
  int<lower=0> N;
  int<lower=0> C[N];
  int <lower=0, upper=1> vaccinated[N];
  real <lower=0> pers_yrs_at_risk[N];
}

parameters {
  real alpha;
  real <lower=-30, upper=1> VE_sym;
}

transformed parameters {
  real <lower=0> mu_c[N];
  
 for (i in 1:N){
    mu_c[i] =(1-VE_sym*vaccinated[i])*exp(alpha)*pers_yrs_at_risk[i];
    
  }
}


model {
  target+= poisson_lpmf(C | mu_c);
}

