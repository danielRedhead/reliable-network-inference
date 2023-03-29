####### THETA RATE #######
PP = 20
TV_test <- seq(0.001, 2, length.out=PP)
G_list = dat_list = vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### FPR test pars
N_id_tv = 100
N_groups_tv = 3
gprobs_tv = c(0.2, 0.5, 0.3)
in_block_tv = 0.02
out_block_tv = 0.002
sr_mu_tv = c(0,0)  
sr_sigma_tv = c(2.1, 0.7) 
sr_rho_tv = 0.5
dr_mu_tv = c(0,0) 
dr_sigma_tv = 1.4
dr_rho_tv = 0.75
predictor_1_tv = rbern(100, 0)
sr_effects_1_tv = c(-0.3, 1.3)
fpr_effects_1_tv = c(0.1, -0.9, 0.0)
rtt_effects_1_tv = c(-0.1, 0.3, -0.3)
theta_effects_1_tv = 0.0
false_positive_rate_tv = c(0.01, 0.01, 0.01)
recall_of_true_ties_tv = c(0.95, 0.95, 0.95)
fpr_sigma_tv = c(0.7, 0.7, 0.7) 
rtt_sigma_tv = c(0.7, 0.7, 0.7)
theta_mean_tv = 0.2
theta_sigma_tv = TV_test[i]  
N_responses_tv = 2
N_periods_tv = 2
decay_curve_tv = 0.0*exp(-seq(1,4, length.out=2))
flow_rate_tv = c(0, 0)

G = simulate_selfreport_network( N_id = N_id_tv, 
                                 N_groups = N_groups_tv, 
                                 group_probs = gprobs_tv,
                                 in_block = in_block_tv, 
                                 out_block = out_block_tv,
                                 sr_mu = sr_mu_tv,  
                                 sr_sigma = sr_sigma_tv, 
                                 sr_rho = sr_rho_tv,
                                 dr_mu = dr_mu_tv,  
                                 dr_sigma = dr_sigma_tv, 
                                 dr_rho = dr_rho_tv,
                                 individual_predictors = data.frame(Status=predictor_1_tv),
                                 dyadic_predictors = NULL,
                                 individual_effects = matrix(sr_effects_1_tv,nrow=2,ncol=1),
                                 dyadic_effects = NULL,
                                 fpr_effects = matrix(fpr_effects_1_tv,nrow=3,ncol=1),
                                 rtt_effects = matrix(rtt_effects_1_tv,nrow=3,ncol=1),
                                 theta_effects = matrix(theta_effects_1_tv,nrow=3,ncol=1),
                                 false_positive_rate = false_positive_rate_tv, 
                                 recall_of_true_ties = recall_of_true_ties_tv,
                                 fpr_sigma = fpr_sigma_tv, 
                                 rtt_sigma = rtt_sigma_tv,
                                 theta_mean = theta_mean_tv, 
                                 theta_sigma = theta_sigma_tv, 
                                 N_responses = N_responses_tv,
                                 N_periods = N_periods_tv, 
                                 decay_curve = decay_curve_tv,
                                 flow_rate = flow_rate_tv
                                 )        

model_dat = make_strand_data(self_report=list(G$reporting_network[,,1],G$reporting_network[,,2]),  group_ids=factor(G$group_id), individual_covariates=data.frame(Status=predictor_1_tv), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_tv <- mclapply(1:PP, function(z){

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

stanfit_tv <- vector("list", PP)
for(i in 1:PP)
stanfit_tv[[i]] <- rstan::read_stan_csv(fit_tv[[i]]$fit$output_files())


save.image(file = "theta-var.RData")










