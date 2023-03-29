####### FALSE POSITIVE RATE #######
PP <- 20
FPR_test <- seq(0.001, 0.02, length.out=PP)
RTT_test <- rev(seq(0.4, 0.999, length.out=PP))

G_list <- dat_list <- vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### FPR test pars
N_id_mix = 100
N_groups_mix = 3
gprobs_mix = c(0.2, 0.5, 0.3)
in_block_mix = 0.02
out_block_mix = 0.002
sr_mu_mix = c(0,0)  
sr_sigma_mix = c(2.1, 0.7) 
sr_rho_mix = 0.5
dr_mu_mix = c(0,0) 
dr_sigma_mix = 1.4
dr_rho_mix = 0.75
predictor_1_mix = rbern(100, 0)
sr_effects_1_mix = c(-0.3, 1.3)
fpr_effects_1_mix = c(0.1, -0.9, 0.0)
rtt_effects_1_mix = c(-0.1, 0.3, -0.3)
theta_effects_1_mix = 0.0
false_positive_rate_mix = c(FPR_test[i], FPR_test[i], FPR_test[i])
recall_of_true_ties_mix = c(RTT_test[i], RTT_test[i], RTT_test[i])
fpr_sigma_mix = c(0.7, 0.7, 0.7) 
rtt_sigma_mix = c(0.7, 0.7, 0.7)
theta_mean_mix = 0 
theta_sigma_mix = 0 
N_responses_mix = 2
N_periods_mix = 2
decay_curve_mix = 0.0*exp(-seq(1,4, length.out=2))
flow_rate_mix = c(0, 0)

G = simulate_selfreport_network( N_id = N_id_mix, 
                                 N_groups = N_groups_mix, 
                                 group_probs = gprobs_mix,
                                 in_block = in_block_mix, 
                                 out_block = out_block_mix,
                                 sr_mu = sr_mu_mix,  
                                 sr_sigma = sr_sigma_mix, 
                                 sr_rho = sr_rho_mix,
                                 dr_mu = dr_mu_mix,  
                                 dr_sigma = dr_sigma_mix, 
                                 dr_rho = dr_rho_mix,
                                 individual_predictors = data.frame(Status=predictor_1_mix),
                                 dyadic_predictors = NULL,
                                 individual_effects = matrix(sr_effects_1_mix,nrow=2,ncol=1),
                                 dyadic_effects = NULL,
                                 fpr_effects = matrix(fpr_effects_1_mix,nrow=3,ncol=1),
                                 rtt_effects = matrix(rtt_effects_1_mix,nrow=3,ncol=1),
                                 theta_effects = matrix(theta_effects_1_mix,nrow=3,ncol=1),
                                 false_positive_rate = false_positive_rate_mix, 
                                 recall_of_true_ties = recall_of_true_ties_mix,
                                 fpr_sigma = fpr_sigma_mix, 
                                 rtt_sigma = rtt_sigma_mix,
                                 theta_mean = theta_mean_mix, 
                                 theta_sigma = theta_sigma_mix, 
                                 N_responses = N_responses_mix,
                                 N_periods = N_periods_mix, 
                                 decay_curve = decay_curve_mix,
                                 flow_rate = flow_rate_mix
                                 )  

model_dat = make_strand_data(self_report=list(G$reporting_network[,,1],G$reporting_network[,,2]),  group_ids=factor(G$group_id), individual_covariates=data.frame(Status=predictor_1_mix), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_mix <- mclapply(1:PP, function(z){

  fit_latent_network_model(data=dat_list[[z]],
                           fpr_regression = ~ 1,
                           rtt_regression = ~ 1,
                           theta_regression = ~ 1,
                           focal_regression = ~ 1,
                           target_regression = ~ 1,
                           dyad_regression = ~ 1,
                           mode="mcmc",
                           return_latent_network=TRUE,
                           stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 1000,
                                                       iter_sampling = 1000, max_treedepth = NULL, adapt_delta = NULL)
                           )
  },
  mc.cores = PP)

stanfit_mix <- vector("list", PP)
for(i in 1:PP)
stanfit_mix[[i]] <- rstan::read_stan_csv(fit_mix[[i]]$fit$output_files())


save.image(file = "mixed-rate.RData")

