####### FALSE POSITIVE RATE #######
PP = 20
SR_sigma_test = seq(0.1, 5, length.out=PP)
G_list = dat_list = vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### FPR test pars
N_id_sr_sigma = 100
N_groups_sr_sigma = 3
gprobs_sr_sigma = c(0.2, 0.5, 0.3)
in_block_sr_sigma = 0.02
out_block_sr_sigma = 0.002
sr_mu_sr_sigma = c(0,0)  
sr_sigma_sr_sigma = c(1.7*SR_sigma_test[i], 1*SR_sigma_test[i]) 
sr_rho_sr_sigma = 0.5
dr_mu_sr_sigma = c(0,0) 
dr_sigma_sr_sigma = 1.4
dr_rho_sr_sigma = 0.75
predictor_1_sr_sigma = rbern(100, 0)
sr_effects_1_sr_sigma = c(-0.3, 1.3)

G = simulate_sbm_plus_srm_network(N_id = N_id_sr_sigma, 
                                 N_groups = N_groups_sr_sigma, 
                                 group_probs = gprobs_sr_sigma,
                                 in_block = in_block_sr_sigma, 
                                 out_block = out_block_sr_sigma,
                                 sr_mu = sr_mu_sr_sigma,  
                                 sr_sigma = sr_sigma_sr_sigma, 
                                 sr_rho = sr_rho_sr_sigma,
                                 dr_mu = dr_mu_sr_sigma,  
                                 dr_sigma = dr_sigma_sr_sigma, 
                                 dr_rho = dr_rho_sr_sigma,
                                 individual_predictors = data.frame(Status=predictor_1_sr_sigma),
                                 dyadic_predictors = NULL,
                                 individual_effects = matrix(sr_effects_1_sr_sigma,nrow=2,ncol=1),
                                 dyadic_effects = NULL
                                 )        

model_dat = make_strand_data(self_report=list(G$network),  group_ids=factor(G$group_id), individual_covariates=data.frame(Status=predictor_1_sr_sigma), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_sr_sigma <- mclapply(1:PP, function(z){

  fit_block_plus_social_relations_model(data=dat_list[[z]],
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

stanfit_sr_sigma <- vector("list", PP)
for(i in 1:PP)
stanfit_sr_sigma[[i]] <- rstan::read_stan_csv(fit_sr_sigma[[i]]$fit$output_files())


save.image(file = "SR_sigma.RData")










