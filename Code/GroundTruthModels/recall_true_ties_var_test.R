####### RECALL OF TRUE TIE VAR #######
PP <- 20
RTV_test <- seq(0.001, 2, length.out=PP)
G_list <- dat_list <- vector("list",PP)

for(i in 1:PP){
  set.seed(420)

####### RTV test pars
N_id_rtv = 100
N_groups_rtv = 3
gprobs_rtv = c(0.32, 0.36, 0.32)
in_block_rtv = 0.008
out_block_rtv = 0.0006
sr_mu_rtv = c(0,0)  
sr_sigma_rtv = c(1.8, 1.2) 
sr_rho_rtv = 0.5
dr_mu_rtv = c(0,0) 
dr_sigma_rtv = 2.2
dr_rho_rtv = 0.65
predictor_1_rtv = rbeta(100, 0.4*1.2, (1-0.4)*1.2) # Status is continuous, but most are low in status
sr_effects_1_rtv = c(2.1, -1.0)                       # High status are more likely to send food, and less likely to recieve
fpr_effects_1_rtv = c(-1.1, -0.9, 0.0)                # high status are less likely to exagerate ties in either direction, and no effect in layer 3
rtt_effects_1_rtv = c(-2.8, -2.5, 0.0)                # high status are less likely to recall who they gave to, but else 0
theta_effects_1_rtv = 0.0                             # No effect of status of question order
false_positive_rate_rtv = c(0.02, 0.02, 0.005)
recall_of_true_ties_rtv = c(0.2, 0.2, 0.99)
fpr_sigma_rtv = c(0.8, 0.8, 0.1) 
rtt_sigma_rtv = c(RTV_test[i], RTV_test[i], 0.1)
theta_mean_rtv = 0.2 
theta_sigma_rtv = 1.8 
N_responses_rtv = 2
N_periods_rtv = 12
effect_max = 2.2
effect_decay = 0.15
decay_curve_rtv = effect_max*exp(-effect_decay *c(1:N_periods_rtv))
flow_rate_rtv = rbeta(N_periods_rtv, 4, 30)

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

ground_dat = vector("list",N_periods_rtv)      
for(k in 1:N_periods_rtv)
ground_dat[[k]] =  G$transfer_network[,,k]                              

model_dat = make_strand_data(self_report=list(G$reporting_network[,,1],G$reporting_network[,,2]), ground_truth=ground_dat, group_ids=factor(G$group_id), 
                             individual_covariates=data.frame(Status=predictor_1_rtv), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_rtv <- mclapply(1:PP, function(z){

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

stanfit_rtv <- vector("list", PP)
for(i in 1:PP)
stanfit_rtv[[i]] <- rstan::read_stan_csv(fit_rtv[[i]]$fit$output_files())

save.image(file = "recall-true-tie-var.RData")


