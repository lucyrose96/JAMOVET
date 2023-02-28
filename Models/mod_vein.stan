// Stan model for  analysis of any infection
// in COVID-19 vaccine trial data
//Nicholas Grassly 25 Aug 2021

data {
  int<lower=0> N;
  int<lower=0> Total[N];
  int <lower=0, upper=1> vaccinated[N];
  real <lower=0> pers_yrs_at_risk[N];
}

parameters {
  real alpha;
  real <lower=-30, upper=1> VE_any;
}

transformed parameters {
  real <lower=0> mu[N];
  
 for (i in 1:N){
      mu[i] =(1-VE_any*vaccinated[i])*exp(alpha)*pers_yrs_at_risk[i];
  }
}


model {
  target+= poisson_lpmf(Total | mu);
}


