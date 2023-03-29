####### RECALL OF TRUE TIE RATE #######
PP <- 20
RTT_test <- seq(0.4, 0.999, length.out=PP)
G_list <- dat_list <- vector("list",PP)

for(i in 1:PP){
   set.seed(420)

####### RTT test pars
N_id_rtt = 100
N_groups_rtt = 3
gprobs_rtt = c(0.32, 0.36, 0.32)
in_block_rtt = 0.008
out_block_rtt = 0.0006
sr_mu_rtt = c(0,0)  
sr_sigma_rtt = c(1.8, 1.2) 
sr_rho_rtt = 0.5
dr_mu_rtt = c(0,0) 
dr_sigma_rtt = 2.2
dr_rho_rtt = 0.65
predictor_1_rtt = rbeta(100, 0.4*1.2, (1-0.4)*1.2) # Status is continuous, but most are low in status
sr_effects_1_rtt = c(2.1, -1.0)                       # High status are more likely to send food, and less likely to recieve
fpr_effects_1_rtt = c(-1.1, -0.9, 0.0)                # high status are less likely to exagerate ties in either direction, and no effect in layer 3
rtt_effects_1_rtt = c(-2.8, -2.5, 0.0)                # high status are less likely to recall who they gave to, but else 0
theta_effects_1_rtt = 0.0                             # No effect of status of question order
false_positive_rate_rtt = c(0.02, 0.02, 0.005)
recall_of_true_ties_rtt = c(RTT_test[i], RTT_test[i], 0.995)
fpr_sigma_rtt = c(0.8, 0.8, 0.1) 
rtt_sigma_rtt = c(1.78, 1.8, 0.1)
theta_mean_rtt = 0.2 
theta_sigma_rtt = 1.8 
N_responses_rtt = 2
N_periods_rtt = 12
effect_max = 2.2
effect_decay = 0.15
decay_curve_rtt = effect_max*exp(-effect_decay *c(1:N_periods_rtt))
flow_rate_rtt = rbeta(N_periods_rtt, 4, 30)


 G = simulate_selfreport_network( N_id = N_id_rtt, 
                                 N_groups = N_groups_rtt, 
                                 group_probs = gprobs_rtt,
                                 in_block = in_block_rtt, 
                                 out_block = out_block_rtt,
                                 sr_mu = sr_mu_rtt,  
                                 sr_sigma = sr_sigma_rtt, 
                                 sr_rho = sr_rho_rtt,
                                 dr_mu = dr_mu_rtt,  
                                 dr_sigma = dr_sigma_rtt, 
                                 dr_rho = dr_rho_rtt,
                                 individual_predictors = data.frame(Status=predictor_1_rtt),
                                 dyadic_predictors = NULL,
                                 individual_effects = matrix(sr_effects_1_rtt,nrow=2,ncol=1),
                                 dyadic_effects = NULL,
                                 fpr_effects = matrix(fpr_effects_1_rtt,nrow=3,ncol=1),
                                 rtt_effects = matrix(rtt_effects_1_rtt,nrow=3,ncol=1),
                                 theta_effects = matrix(theta_effects_1_rtt,nrow=3,ncol=1),
                                 false_positive_rate = false_positive_rate_rtt, 
                                 recall_of_true_ties = recall_of_true_ties_rtt,
                                 fpr_sigma = fpr_sigma_rtt, 
                                 rtt_sigma = rtt_sigma_rtt,
                                 theta_mean = theta_mean_rtt, 
                                 theta_sigma = theta_sigma_rtt, 
                                 N_responses = N_responses_rtt,
                                 N_periods = N_periods_rtt, 
                                 decay_curve = decay_curve_rtt,
                                 flow_rate = flow_rate_rtt
                                 )  

model_dat = make_strand_data(self_report=list(G$reporting_network[,,1],G$reporting_network[,,2]),  group_ids=factor(G$group_id), individual_covariates=data.frame(Status=predictor_1_rtt), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}

fit_rtt <- mclapply(1:PP, function(z){

  fit_latent_network_model(data=dat_list[[z]],
                           fpr_regression = ~ Status,
                           rtt_regression = ~ Status,
                           theta_regression = ~ 1,
                           focal_regression = ~ Status,
                           target_regression = ~ Status,
                           dyad_regression = ~ 1,
                           mode="mcmc",
                           return_latent_network=TRUE,
                           stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 1000,
                                                       iter_sampling = 1000, max_treedepth = NULL, adapt_delta = NULL)
                           )
  },
  mc.cores = PP)

stanfit_rtt <- vector("list", PP)
for(i in 1:PP)
stanfit_rtt[[i]] <- rstan::read_stan_csv(fit_rtt[[i]]$fit$output_files())


save.image(file = "recall-true-tie-rate.RData")

