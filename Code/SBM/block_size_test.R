
library(dirmult)

####### SBM BLOCK SIZE #######
PP = 20
block_prob_rate_test = seq(0.01, 0.9, length.out=PP)
G_list = dat_list = vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### n blocks test pars
N_id_block_size = 100
in_block_block_size = 0.05
out_block_block_size = 0.01
N_groups_block_size = 5
group_probs_block_size = c(block_prob_rate_test[i], c(0.25, 0.25, 0.25, 0.25)*(1-block_prob_rate_test[i]))
sr_mu_block_size = c(0,0)  
sr_sigma_block_size = c(2.1,0.7) 
sr_rho_block_size = 0.5
dr_mu_block_size = c(0,0) 
dr_sigma_block_size = 1.4
dr_rho_block_size = 0.75
predictor_1_block_size = rbern(100, 0)
sr_effects_1_block_size = c(-0.3, 1.3)

G = simulate_sbm_network(N_id = N_id_block_size,     
                         N_groups = N_groups_block_size,
                         group_probs = group_probs_block_size,                         
                         in_block = in_block_block_size,
                         out_block = out_block_block_size,
                         sr_mu = sr_mu_block_size,  
                         sr_sigma = sr_sigma_block_size, 
                         sr_rho = sr_rho_block_size,
                         dr_mu = dr_mu_block_size,  
                         dr_sigma = dr_sigma_block_size, 
                         dr_rho = dr_rho_block_size,
                         individual_predictors = data.frame(Status=predictor_1_block_size),
                         dyadic_predictors = NULL,
                         individual_effects = matrix(sr_effects_1_block_size,nrow=2,ncol=1),
                         dyadic_effects = NULL
                         )        

model_dat = make_strand_data(self_report=list(G$network),  group_ids=factor(G$group_ids), individual_covariates=data.frame(Status=predictor_1_block_size), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_block_size <- mclapply(1:PP, function(z){

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

stanfit_block_size <- vector("list", PP)
for(i in 1:PP)
stanfit_block_size[[i]] <- rstan::read_stan_csv(fit_block_size[[i]]$fit$output_files())


save.image(file = "block_size_rate.RData")










