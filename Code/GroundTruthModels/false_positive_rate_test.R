####### FALSE POSITIVE RATE #######
PP = 20
FPR_test = seq(0.001, 0.05, length.out=PP)
G_list = dat_list = vector("list",PP)

for(i in 1:PP){
  set.seed(420)

####### FPR test pars
N_id_fpr = 100
N_groups_fpr = 3
gprobs_fpr = c(0.32, 0.36, 0.32)
in_block_fpr = 0.008
out_block_fpr = 0.0006
sr_mu_fpr = c(0,0)  
sr_sigma_fpr = c(1.8, 1.2) 
sr_rho_fpr = 0.5
dr_mu_fpr = c(0,0) 
dr_sigma_fpr = 2.2
dr_rho_fpr = 0.65
predictor_1_fpr = rbeta(100, 0.4*1.2, (1-0.4)*1.2) # Status is continuous, but most are low in status
sr_effects_1_fpr = c(2.1, -1.0)                       # High status are more likely to send food, and less likely to recieve
fpr_effects_1_fpr = c(-1.1, -0.9, 0.0)                # high status are less likely to exagerate ties in either direction, and no effect in layer 3
rtt_effects_1_fpr = c(-2.8, -2.5, 0.0)                # high status are less likely to recall who they gave to, but else 0
theta_effects_1_fpr = 0.0                             # No effect of status of question order
false_positive_rate_fpr = c(FPR_test[i], FPR_test[i], 0.005)
recall_of_true_ties_fpr = c(0.2, 0.2, 0.99)
fpr_sigma_fpr = c(0.8, 0.8, 0.1) 
rtt_sigma_fpr = c(1.78, 1.8, 0.1)
theta_mean_fpr = 0.2 
theta_sigma_fpr = 1.8 
N_responses_fpr = 2
N_periods_fpr = 12
effect_max = 2.2
effect_decay = 0.15
decay_curve_fpr = effect_max*exp(-effect_decay *c(1:N_periods_fpr))
flow_rate_fpr = rbeta(N_periods_fpr, 4, 30)

G = simulate_selfreport_network( N_id = N_id_fpr, 
                                 N_groups = N_groups_fpr, 
                                 group_probs = gprobs_fpr,
                                 in_block = in_block_fpr, 
                                 out_block = out_block_fpr,
                                 sr_mu = sr_mu_fpr,  
                                 sr_sigma = sr_sigma_fpr, 
                                 sr_rho = sr_rho_fpr,
                                 dr_mu = dr_mu_fpr,  
                                 dr_sigma = dr_sigma_fpr, 
                                 dr_rho = dr_rho_fpr,
                                 individual_predictors = data.frame(Status=predictor_1_fpr),
                                 dyadic_predictors = NULL,
                                 individual_effects = matrix(sr_effects_1_fpr,nrow=2,ncol=1),
                                 dyadic_effects = NULL,
                                 fpr_effects = matrix(fpr_effects_1_fpr,nrow=3,ncol=1),
                                 rtt_effects = matrix(rtt_effects_1_fpr,nrow=3,ncol=1),
                                 theta_effects = matrix(theta_effects_1_fpr,nrow=3,ncol=1),
                                 false_positive_rate = false_positive_rate_fpr, 
                                 recall_of_true_ties = recall_of_true_ties_fpr,
                                 fpr_sigma = fpr_sigma_fpr, 
                                 rtt_sigma = rtt_sigma_fpr,
                                 theta_mean = theta_mean_fpr, 
                                 theta_sigma = theta_sigma_fpr, 
                                 N_responses = N_responses_fpr,
                                 N_periods = N_periods_fpr, 
                                 decay_curve = decay_curve_fpr,
                                 flow_rate = flow_rate_fpr
                                 )    

ground_dat = vector("list",N_periods_fpr)      
for(k in 1:N_periods_fpr)
ground_dat[[k]] =  G$transfer_network[,,k]                              

model_dat = make_strand_data(self_report=list(G$reporting_network[,,1],G$reporting_network[,,2]), ground_truth=ground_dat, group_ids=factor(G$group_id), 
                             individual_covariates=data.frame(Status=predictor_1_fpr), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_fpr <- mclapply(1:PP, function(z){

  fit_latent_network_plus_flows_model(data=dat_list[[z]],
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

stanfit_fpr <- vector("list", PP)
for(i in 1:PP)
stanfit_fpr[[i]] <- rstan::read_stan_csv(fit_fpr[[i]]$fit$output_files())


save.image(file = "false-positive-rate.RData")










