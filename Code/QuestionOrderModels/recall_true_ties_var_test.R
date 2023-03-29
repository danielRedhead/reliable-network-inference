####### RECALL OF TRUE TIE VAR #######
PP <- 20
RTV_test <- seq(0.001, 2, length.out=PP)
G_list <- dat_list <- vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### RTV test pars
N_id_rtv = 100
N_groups_rtv = 3
gprobs_rtv = c(0.2, 0.5, 0.3)
in_block_rtv = 0.02
out_block_rtv = 0.002
sr_mu_rtv = c(0,0)  
sr_sigma_rtv = c(2.1, 0.7) 
sr_rho_rtv = 0.5
dr_mu_rtv = c(0,0) 
dr_sigma_rtv = 1.4
dr_rho_rtv = 0.75
predictor_1_rtv = rbern(100, 0)
sr_effects_1_rtv = c(-0.3, 1.3)
fpr_effects_1_rtv = c(0.1, -0.9, 0.0)
rtt_effects_1_rtv = c(-0.1, 0.3, -0.3)
theta_effects_1_rtv = 0.0
false_positive_rate_rtv = c(0.01, 0.01, 0.01)
recall_of_true_ties_rtv = c(0.95, 0.95, 0.95)
fpr_sigma_rtv = c(0.7, 0.7, 0.7) 
rtt_sigma_rtv = c(RTV_test[i], RTV_test[i], RTV_test[i])
theta_mean_rtv = 0.2 
theta_sigma_rtv = 0.7 
N_responses_rtv = 2
N_periods_rtv = 2
decay_curve_rtv = 0.0*exp(-seq(1,4, length.out=2))
flow_rate_rtv = c(0, 0)

G = simulate_selfreport_network( N_id = N_id_rtv, 
                                 N_groups = N_groups_rtv, 
                                 group_probs = gprobs_rtv,
                                 in_block = in_block_rtv, 
                                 out_block = out_block_rtv,
                                 sr_mu = sr_mu_rtv,  
                                 sr_sigma = sr_sigma_rtv, 
                                 sr_rho = sr_rho_rtv,
                                 dr_mu = dr_mu_rtv,  
                                 dr_sigma = dr_sigma_rtv, 
                                 dr_rho = dr_rho_rtv,
                                 individual_predictors = data.frame(Status=predictor_1_rtv),
                                 dyadic_predictors = NULL,
                                 individual_effects = matrix(sr_effects_1_rtv,nrow=2,ncol=1),
                                 dyadic_effects = NULL,
                                 fpr_effects = matrix(fpr_effects_1_rtv,nrow=3,ncol=1),
                                 rtt_effects = matrix(rtt_effects_1_rtv,nrow=3,ncol=1),
                                 theta_effects = matrix(theta_effects_1_rtv,nrow=3,ncol=1),
                                 false_positive_rate = false_positive_rate_rtv, 
                                 recall_of_true_ties = recall_of_true_ties_rtv,
                                 fpr_sigma = fpr_sigma_rtv, 
                                 rtt_sigma = rtt_sigma_rtv,
                                 theta_mean = theta_mean_rtv, 
                                 theta_sigma = theta_sigma_rtv, 
                                 N_responses = N_responses_rtv,
                                 N_periods = N_periods_rtv, 
                                 decay_curve = decay_curve_rtv,
                                 flow_rate = flow_rate_rtv
                                 )  

model_dat = make_strand_data(self_report=list(G$reporting_network[,,1],G$reporting_network[,,2]),  group_ids=factor(G$group_id), individual_covariates=data.frame(Status=predictor_1_rtv), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_rtv <- mclapply(1:PP, function(z){

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

stanfit_rtv <- vector("list", PP)
for(i in 1:PP)
stanfit_rtv[[i]] <- rstan::read_stan_csv(fit_rtv[[i]]$fit$output_files())

save.image(file = "recall-true-tie-var.RData")


