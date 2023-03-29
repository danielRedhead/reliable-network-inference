####### FALSE POSITIVE VAR #######
PP <- 20
FPV_test <- seq(0.001, 2, length.out=PP)
G_list <- dat_list <- vector("list",PP)

for(i in 1:PP){
  set.seed(1)

# FPV test params
N_id_fpv = 100
N_groups_fpv = 3
gprobs_fpv = c(0.2, 0.5, 0.3)
in_block_fpv = 0.02
out_block_fpv = 0.002
sr_mu_fpv = c(0,0)  
sr_sigma_fpv = c(2.1, 0.7) 
sr_rho_fpv = 0.5
dr_mu_fpv = c(0,0) 
dr_sigma_fpv = 1.4
dr_rho_fpv = 0.75
predictor_1_fpv = rbern(100, 0)
sr_effects_1_fpv = c(-0.3, 1.3)
fpr_effects_1_fpv = c(0.1, -0.9, 0.0)
rtt_effects_1_fpv = c(-0.1, 0.3, -0.3)
theta_effects_1_fpv = 0.0
false_positive_rate_fpv = c(0.01, 0.01, 0.01)
recall_of_true_ties_fpv = c(0.95, 0.95, 0.95)
fpr_sigma_fpv = c(FPV_test[i], FPV_test[i], FPV_test[i]) 
rtt_sigma_fpv = c(0.7, 0.7, 0.7)
theta_mean_fpv = 0 
theta_sigma_fpv = 0 
N_responses_fpv = 2
N_periods_fpv = 2
decay_curve_fpv = 0.0*exp(-seq(1,4, length.out=2))
flow_rate_fpv = c(0, 0)

G = simulate_selfreport_network( N_id = N_id_fpv, 
                                 N_groups = N_groups_fpv, 
                                 group_probs = gprobs_fpv,
                                 in_block = in_block_fpv, 
                                 out_block = out_block_fpv,
                                 sr_mu = sr_mu_fpv,  
                                 sr_sigma = sr_sigma_fpv, 
                                 sr_rho = sr_rho_fpv,
                                 dr_mu = dr_mu_fpv,  
                                 dr_sigma = dr_sigma_fpv, 
                                 dr_rho = dr_rho_fpv,
                                 individual_predictors = data.frame(Status=predictor_1_fpv),
                                 dyadic_predictors = NULL,
                                 individual_effects = matrix(sr_effects_1_fpv,nrow=2,ncol=1),
                                 dyadic_effects = NULL,
                                 fpr_effects = matrix(fpr_effects_1_fpv,nrow=3,ncol=1),
                                 rtt_effects = matrix(rtt_effects_1_fpv,nrow=3,ncol=1),
                                 theta_effects = matrix(theta_effects_1_fpv,nrow=3,ncol=1),
                                 false_positive_rate = false_positive_rate_fpv, 
                                 recall_of_true_ties = recall_of_true_ties_fpv,
                                 fpr_sigma = fpr_sigma_fpv, 
                                 rtt_sigma = rtt_sigma_fpv,
                                 theta_mean = theta_mean_fpv, 
                                 theta_sigma = theta_sigma_fpv, 
                                 N_responses = N_responses_fpv,
                                 N_periods = N_periods_fpv, 
                                 decay_curve = decay_curve_fpv,
                                 flow_rate = flow_rate_fpv
                                 )  

model_dat = make_strand_data(self_report=list(G$reporting_network[,,1],G$reporting_network[,,2]),  group_ids=factor(G$group_id), individual_covariates=data.frame(Status=predictor_1_fpv), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_fpv <- mclapply(1:PP, function(z){

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

stanfit_fpv <- vector("list", PP)
for(i in 1:PP)
stanfit_fpv[[i]] <- rstan::read_stan_csv(fit_fpv[[i]]$fit$output_files())


save.image(file = "false-positive-var.RData")

