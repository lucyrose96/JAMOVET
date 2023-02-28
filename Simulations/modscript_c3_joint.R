#setup 
library(rstan)
library(foreach)
library(doParallel)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
numCores <- 16
registerDoParallel(numCores)

model <- stan_model(file = 'mod_bias.stan')
c3<-read.csv("c3.csv")
num_sims<-length(unique(c3$sim))
sims_vec=unique(c3$sim)


timetaken<-system.time(
  results <- foreach(i = sims_vec, .combine=rbind, .packages = c("rstan", "dplyr"))%dopar%{
    
    d                     <-data.frame(Total=c3[which(c3$sim == i),]$Total,
                                       C=c3[which(c3$sim == i),]$C,
                                       A=c3[which(c3$sim == i),]$A,
                                       fp=c3[which(c3$sim == i),]$fp_rate*c3$time_at_risk_years,
                                       rel_sens_asym=c3[which(c3$sim == i),]$rel_sens_asym,
                                       time_at_risk_years=c3[which(c3$sim == i),]$time_at_risk_years,
                                       vaccinated=c3[which(c3$sim == i),]$arm_01,
                                       intercept=c3[which(c3$sim == i),]$intercept
    )
    trial_dat_c3 <- list(N=nrow(d), 
                         Total=d$Total,
                         C=d$C,
                         A=d$A,
                         vaccinated=d$vaccinated,
                         pers_yrs_at_risk=d$time_at_risk_years,
                         fp=d$fp,
                         rel_sens_asym=d$rel_sens_asym
                         
    )
    as.data.frame(sampling(object = model, data = trial_dat_c3, chains = 4, cores = 4))%>%
      select(c("VE_in", "VE_pr", "VE_sym", "VE_asym", "FOI_p", "FOI_v", "PS_p", "PS_v"))%>%
      summarise(
        mean_VEin=mean(VE_in),
        median_VEin=median(VE_in),
        mean_VEpr=mean(VE_pr),
        median_VEpr=median(VE_pr),
        mean_VEsym=mean(VE_sym),
        median_VEsym=median(VE_sym),
        mean_VEasym=mean(VE_asym),
        median_VEasym=median(VE_asym),
        mean_FOI_p=mean(FOI_p),
        median_FOI_p=median(FOI_p),
        mean_FOI_v=mean(FOI_v),
        median_FOI_v=median(FOI_v),
        mean_PS_p=mean(PS_p),
        median_PS_p=median(PS_p),
        mean_PS_v=mean(PS_v),
        median_PS_v=median(PS_v),
      )
 }
)
write.csv(results, "results_c3_joint.csv", row.names = FALSE)

