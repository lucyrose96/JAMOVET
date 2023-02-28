### Setup data for Rstan simulation model - biased specificity ###


#1. Setup
library(tidyverse)
set.seed(1)
lci<-function(x){quantile(x, 0.025)}
uci<-function(x){quantile(x, 0.975)}

#2. Set common parameters
  num_subjects<-10000 #number of subjects
  num_sims=5 #number of simulations (5 for demonstration, 1000 were run for the paper)
  lam<-c(0.1) #force of infection
  s<-c(0.6) #probability of symptoms
  VE_in<-0.5 #VE against infection
  VE_pr<-0.5 #VE against progression to symptoms

#2. Simulate unbiased data
  # a. set test parameters
  
  sens_sym<-1 #sensitivity to symptomatic infections
  sens_asym<-1 #sensitivity to asymptomatic infections
  rel_sens_asym<-sens_asym/sens_sym #relative sensitivity
  spec_asym<-0.999 #specificity to asymptomatic infections
  follow_up=1 # follow up time in years
  num_tests=52 #number of tests over all of follow up (non-censored)
  test_freq= num_tests/follow_up # number of tests taken each year times by duration of follow-up
  fp_rate<-test_freq*(1-spec_asym) #rate of false positives
  
  # b. create data frame for vaccine and placebo arms
  c2_v<-data.frame(matrix(nrow = (num_subjects/2)*num_sims, ncol=4, data = NA ))
  colnames(c2_v)<-c("arm","A","C", "sim")
  c2_v$arm<-"V"
  c2_v$arm_01<-1
  c2_v$sim<-rep(1:num_sims, each =num_subjects/2)

  c2_p<-data.frame(matrix(nrow = (num_subjects/2)*num_sims, ncol=4, data = NA ))
  colnames(c2_p)<-c("arm","A","C", "sim")
  c2_p$arm<-"P"
  c2_p$arm_01<-0
  c2_p$sim<-rep(1:num_sims, each =num_subjects/2)

# c. simulate for vaccine and placebo arms
  c2_p<-c2_p%>%
  mutate(time_to_asym_years= rexp(nrow(c2_p), lam*(1-s)*sens_asym),
         time_to_fp_years=pmin(rexp(nrow(c2_p), test_freq*(1-spec_asym)), Inf, na.rm = TRUE),
         time_to_sym_years= rexp(nrow(c2_p), lam*s*sens_sym),
         time_at_risk_years= pmin(time_to_asym_years, time_to_fp_years, time_to_sym_years, follow_up),
         A=ifelse(time_to_asym_years == time_at_risk_years | time_to_fp_years == time_at_risk_years, 1, 0),
         C=ifelse(time_to_sym_years == time_at_risk_years, 1, 0),
         FP_a=ifelse(time_to_fp_years == time_at_risk_years, 1, 0),
         TP_a=ifelse(time_to_asym_years == time_at_risk_years, 1, 0),
         TP=ifelse(time_to_asym_years == time_at_risk_years | time_to_sym_years == time_at_risk_years, 1, 0),
         Total = C+A,
         None=ifelse(Total == 0, 1, 0)
  )

  c2_v<-c2_v%>%
  mutate(time_to_asym_years= rexp(nrow(c2_p), (1-VE_in)*lam*((1-s)+s*VE_pr)*sens_asym ),
         time_to_fp_years=pmin(rexp(nrow(c2_p), test_freq*(1-spec_asym)), Inf, na.rm = TRUE),
         time_to_sym_years= rexp(nrow(c2_p), (1-VE_in)*(1-VE_pr)*lam*s*sens_sym),
         time_at_risk_years= pmin(time_to_asym_years, time_to_fp_years, time_to_sym_years, follow_up),
         A=ifelse(time_to_asym_years == time_at_risk_years | time_to_fp_years == time_at_risk_years, 1, 0),
         C=ifelse(time_to_sym_years == time_at_risk_years, 1, 0),
         FP_a=ifelse(time_to_fp_years == time_at_risk_years, 1, 0),
         TP_a=ifelse(time_to_asym_years == time_at_risk_years, 1, 0),
         TP=ifelse(time_to_asym_years == time_at_risk_years | time_to_sym_years == time_at_risk_years, 1, 0),
         Total = C+A,
         None=ifelse(Total == 0, 1, 0)
  )

# e. join datasets
  c2<-rbind(c2_p, c2_v)
  c2$fp_rate<-fp_rate
  c2$rel_sens_asym<-rel_sens_asym
  c2$intercept=1
  c2$outcome="None"
  c2$outcome[which(c2$FP_a == 1)]<-"False asymptomatic"
  c2$outcome[which(c2$TP_a == 1)]<-"True asymptomatic"
  c2$outcome[which(c2$C == 1)]<-"True symptomatic"
  
# f. summarise
  
summary_c2<-  c2 %>%
    group_by(sim) %>%
    summarise(FP_a = sum(FP_a),
              TP_a = sum(TP_a),
              C = sum(C),
              TP = sum(TP_a + C),
              None=sum(None)
              ) %>%
    summarise(
      FP_a_mean = mean(FP_a),
      FP_a_lci = lci(FP_a),
      FP_a_uci = uci(FP_a),
      TP_a_mean = mean(TP_a),
      TP_a_lci = lci(TP_a),
      TP_a_uci = uci(TP_a),
      C_mean = mean(C),
      C_lci = lci(C),
      C_uci = uci(C),
      None_mean=mean(None),
      None_lci=lci(None),
      None_uci=uci(None),
      TP_mean=mean(TP),
      TP_lci=lci(TP),
      TP_uci=uci(TP),
      PPV_mean= mean( (TP)/(TP + FP_a) ),
      PPV_lci= lci( (TP)/(TP + FP_a) ),
      PPV_uci= uci( (TP)/(TP + FP_a) ))
summary_c2

summary_c2_v<-  c2_v %>%
  group_by(sim) %>%
  summarise(FP_a = sum(FP_a),
            TP_a = sum(TP_a),
            C = sum(C),
            TP = sum(TP_a + C),
            None=sum(None)
  ) %>%
  summarise(
    FP_a_mean = mean(FP_a),
    FP_a_lci = lci(FP_a),
    FP_a_uci = uci(FP_a),
    TP_a_mean = mean(TP_a),
    TP_a_lci = lci(TP_a),
    TP_a_uci = uci(TP_a),
    C_mean = mean(C),
    C_lci = lci(C),
    C_uci = uci(C),
    None_mean=mean(None),
    None_lci=lci(None),
    None_uci=uci(None),
    TP_mean=mean(TP),
    TP_lci=lci(TP),
    TP_uci=uci(TP),
    PPV_mean= mean( (TP)/(TP + FP_a) ),
    PPV_lci= lci( (TP)/(TP + FP_a) ),
    PPV_uci= uci( (TP)/(TP + FP_a) ))
summary_c2_v

summary_c2_p<-  c2_p %>%
  group_by(sim) %>%
  summarise(FP_a = sum(FP_a),
            TP_a = sum(TP_a),
            C = sum(C),
            TP = sum(TP_a + C),
            None=sum(None)
  ) %>%
  summarise(
    FP_a_mean = mean(FP_a),
    FP_a_lci = lci(FP_a),
    FP_a_uci = uci(FP_a),
    TP_a_mean = mean(TP_a),
    TP_a_lci = lci(TP_a),
    TP_a_uci = uci(TP_a),
    C_mean = mean(C),
    C_lci = lci(C),
    C_uci = uci(C),
    None_mean=mean(None),
    None_lci=lci(None),
    None_uci=uci(None),
    TP_mean=mean(TP),
    TP_lci=lci(TP),
    TP_uci=uci(TP),
    PPV_mean= mean( (TP)/(TP + FP_a) ),
    PPV_lci= lci( (TP)/(TP + FP_a) ),
    PPV_uci= uci( (TP)/(TP + FP_a) ))
summary_c2_p


#4. save simulated dataset in both formats
write.csv(c2, "c2.csv", row.names = FALSE)
