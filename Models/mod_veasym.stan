// Stan model for  analysis of any infection
// in COVID-19 vaccine trial data
//Nicholas Grassly 25 Aug 2021

data {
  int<lower=0> N;
  int<lower=0> A[N];
  int <lower=0, upper=1> vaccinated[N];
  real <lower=0> pers_yrs_at_risk[N];
}

parameters {
  real alpha;
  real <lower=-50, upper=1> VE_asym;
}

transformed parameters {
  real <lower=0> mu_a[N];
  
 for (i in 1:N){
    mu_a[i] =(1-VE_asym*vaccinated[i])*exp(alpha)*pers_yrs_at_risk[i];
    
  }
}


model {
  target+= poisson_lpmf(A | mu_a);
}

