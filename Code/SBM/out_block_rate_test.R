####### SBM OUT-BLOCK RATE #######
PP = 20
out_block_rate_test = seq(0.005, 0.03, length.out=PP)
G_list = dat_list = vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### out-block test pars
N_id_out_block = 100
in_block_out_block = 0.04
out_block_out_block = out_block_rate_test[i]
N_groups_out_block = 3
group_probs_out_block = c(0.2, 0.5, 0.3)
sr_mu_out_block = c(0,0)  
sr_sigma_out_block = c(2.1,0.7) 
sr_rho_out_block = 0.5
dr_mu_out_block = c(0,0) 
dr_sigma_out_block = 1.4
dr_rho_out_block = 0.75
predictor_1_out_block = rbern(100, 0)
sr_effects_1_out_block = c(-0.3, 1.3)

G = simulate_sbm_network(N_id = N_id_out_block,     
                         N_groups = N_groups_out_block,
                         group_probs = group_probs_out_block,                         
                         in_block = in_block_out_block,
                         out_block = out_block_out_block,
                         sr_mu = sr_mu_out_block,  
                         sr_sigma = sr_sigma_out_block, 
                         sr_rho = sr_rho_out_block,
                         dr_mu = dr_mu_out_block,  
                         dr_sigma = dr_sigma_out_block, 
                         dr_rho = dr_rho_out_block,
                         individual_predictors = data.frame(Status=predictor_1_out_block),
                         dyadic_predictors = NULL,
                         individual_effects = matrix(sr_effects_1_out_block,nrow=2,ncol=1),
                         dyadic_effects = NULL
                         )        

model_dat = make_strand_data(self_report=list(G$network),  group_ids=factor(G$group_ids), individual_covariates=data.frame(Status=predictor_1_out_block), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_out_block <- mclapply(1:PP, function(z){

  fit_block_model(data=dat_list[[z]],
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

stanfit_out_block <- vector("list", PP)
for(i in 1:PP)
stanfit_out_block[[i]] <- rstan::read_stan_csv(fit_out_block[[i]]$fit$output_files())


save.image(file = "out_block_rate.RData")










