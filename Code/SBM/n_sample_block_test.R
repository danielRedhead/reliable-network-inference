####### SBM N-SAMPLE #######
PP = 20
n_sample_test = round(seq(25, 350, length.out=PP),0)
G_list = dat_list = vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### FPR test pars
N_id_n_sample = n_sample_test[i]
in_block_n_sample = 0.021
N_groups_n_sample = 3
group_probs_n_sample = c(0.2, 0.5, 0.3)
out_block_n_sample = 0.007
sr_mu_n_sample = c(0,0)  
sr_sigma_n_sample = c(2.1,0.7) 
sr_rho_n_sample = 0.5
dr_mu_n_sample = c(0,0) 
dr_sigma_n_sample = 1.4
dr_rho_n_sample = 0.75
predictor_1_n_sample = rbern(n_sample_test[i], 0)
sr_effects_1_n_sample = c(-0.3, 1.3)

G = simulate_sbm_network(N_id = N_id_n_sample,     
                         N_groups = N_groups_n_sample,
                         group_probs = group_probs_n_sample,                         
                         in_block = in_block_n_sample,
                         out_block = out_block_n_sample,
                         sr_mu = sr_mu_n_sample,  
                         sr_sigma = sr_sigma_n_sample, 
                         sr_rho = sr_rho_n_sample,
                         dr_mu = dr_mu_n_sample,  
                         dr_sigma = dr_sigma_n_sample, 
                         dr_rho = dr_rho_n_sample,
                         individual_predictors = data.frame(Status=predictor_1_n_sample),
                         dyadic_predictors = NULL,
                         individual_effects = matrix(sr_effects_1_n_sample,nrow=2,ncol=1),
                         dyadic_effects = NULL
                         )        

model_dat = make_strand_data(self_report=list(G$network),  group_ids=factor(G$group_ids), individual_covariates=data.frame(Status=predictor_1_n_sample), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_n_sample <- mclapply(1:PP, function(z){

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

stanfit_n_sample <- vector("list", PP)
for(i in 1:PP)
stanfit_n_sample[[i]] <- rstan::read_stan_csv(fit_n_sample[[i]]$fit$output_files())


save.image(file = "n_sample.RData")










