####### THETA RATE #######
PP = 20
TR_test <- seq(0.001, 0.66, length.out=PP)
G_list = dat_list = vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### FPR test pars
N_id_tr = 100
N_groups_tr = 3
gprobs_tr = c(0.2, 0.5, 0.3)
in_block_tr = 0.02
out_block_tr = 0.002
sr_mu_tr = c(0,0)  
sr_sigma_tr = c(2.1, 0.7) 
sr_rho_tr = 0.5
dr_mu_tr = c(0,0) 
dr_sigma_tr = 1.4
dr_rho_tr = 0.75
predictor_1_tr = rbern(100, 0)
sr_effects_1_tr = c(-0.3, 1.3)
fpr_effects_1_tr = c(0.1, -0.9, 0.0)
rtt_effects_1_tr = c(-0.1, 0.3, -0.3)
theta_effects_1_tr = 0.0
false_positive_rate_tr = c(0.01, 0.01, 0.01)
recall_of_true_ties_tr = c(0.95, 0.95, 0.95)
fpr_sigma_tr = c(0.7, 0.7, 0.7) 
rtt_sigma_tr = c(0.7, 0.7, 0.7)
theta_mean_tr = TR_test[i] 
theta_sigma_tr = 0.7 
N_responses_tr = 2
N_periods_tr = 2
decay_curve_tr = 0.0*exp(-seq(1,4, length.out=2))
flow_rate_tr = c(0, 0)

G = simulate_selfreport_network( N_id = N_id_tr, 
                                 N_groups = N_groups_tr, 
                                 group_probs = gprobs_tr,
                                 in_block = in_block_tr, 
                                 out_block = out_block_tr,
                                 sr_mu = sr_mu_tr,  
                                 sr_sigma = sr_sigma_tr, 
                                 sr_rho = sr_rho_tr,
                                 dr_mu = dr_mu_tr,  
                                 dr_sigma = dr_sigma_tr, 
                                 dr_rho = dr_rho_tr,
                                 individual_predictors = data.frame(Status=predictor_1_tr),
                                 dyadic_predictors = NULL,
                                 individual_effects = matrix(sr_effects_1_tr,nrow=2,ncol=1),
                                 dyadic_effects = NULL,
                                 fpr_effects = matrix(fpr_effects_1_tr,nrow=3,ncol=1),
                                 rtt_effects = matrix(rtt_effects_1_tr,nrow=3,ncol=1),
                                 theta_effects = matrix(theta_effects_1_tr,nrow=3,ncol=1),
                                 false_positive_rate = false_positive_rate_tr, 
                                 recall_of_true_ties = recall_of_true_ties_tr,
                                 fpr_sigma = fpr_sigma_tr, 
                                 rtt_sigma = rtt_sigma_tr,
                                 theta_mean = theta_mean_tr, 
                                 theta_sigma = theta_sigma_tr, 
                                 N_responses = N_responses_tr,
                                 N_periods = N_periods_tr, 
                                 decay_curve = decay_curve_tr,
                                 flow_rate = flow_rate_tr
                                 )        

model_dat = make_strand_data(self_report=list(G$reporting_network[,,1],G$reporting_network[,,2]),  group_ids=factor(G$group_id), individual_covariates=data.frame(Status=predictor_1_tr), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_tr <- mclapply(1:PP, function(z){

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

stanfit_tr <- vector("list", PP)
for(i in 1:PP)
stanfit_tr[[i]] <- rstan::read_stan_csv(fit_tr[[i]]$fit$output_files())


save.image(file = "theta-rate.RData")










