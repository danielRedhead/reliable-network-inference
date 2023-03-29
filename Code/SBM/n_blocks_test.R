
library(dirmult)

####### SBM N BLOCKS #######
PP = 20
n_blocks_rate_test = seq(2, 21, length.out=PP)
G_list = dat_list = vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### n blocks test pars
N_id_n_block = 100
in_block_n_block = 0.05
out_block_n_block = 0.01
N_groups_n_block = n_blocks_rate_test[i]
group_probs_n_block = rdirichlet(1, rep(1,n_blocks_rate_test[i]))
sr_mu_n_block = c(0,0)  
sr_sigma_n_block = c(2.1,0.7) 
sr_rho_n_block = 0.5
dr_mu_n_block = c(0,0) 
dr_sigma_n_block = 1.4
dr_rho_n_block = 0.75
predictor_1_n_block = rbern(100, 0)
sr_effects_1_n_block = c(-0.3, 1.3)

G = simulate_sbm_network(N_id = N_id_n_block,     
                         N_groups = N_groups_n_block,
                         group_probs = group_probs_n_block,                         
                         in_block = in_block_n_block,
                         out_block = out_block_n_block,
                         sr_mu = sr_mu_n_block,  
                         sr_sigma = sr_sigma_n_block, 
                         sr_rho = sr_rho_n_block,
                         dr_mu = dr_mu_n_block,  
                         dr_sigma = dr_sigma_n_block, 
                         dr_rho = dr_rho_n_block,
                         individual_predictors = data.frame(Status=predictor_1_n_block),
                         dyadic_predictors = NULL,
                         individual_effects = matrix(sr_effects_1_n_block,nrow=2,ncol=1),
                         dyadic_effects = NULL
                         )        

model_dat = make_strand_data(self_report=list(G$network),  group_ids=factor(G$group_ids), individual_covariates=data.frame(Status=predictor_1_n_block), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_n_block <- mclapply(1:PP, function(z){

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

stanfit_n_block <- vector("list", PP)
for(i in 1:PP)
stanfit_n_block[[i]] <- rstan::read_stan_csv(fit_n_block[[i]]$fit$output_files())


save.image(file = "n_block_rate.RData")










