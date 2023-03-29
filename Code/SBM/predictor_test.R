####### SBM OUT-BLOCK RATE #######
PP = 20
pred_test = seq(0.001, 1, length.out=PP)
G_list = dat_list = vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### out-block test pars
N_id_pred = 100
in_block_pred = 0.04
out_block_pred = 0.007 
N_groups_pred = 3
group_probs_pred = c(0.2, 0.5, 0.3)
sr_mu_pred = c(0,0)  
sr_sigma_pred = c(2.1,0.7) 
sr_rho_pred = 0.5
dr_mu_pred = c(0,0) 
dr_sigma_pred = 1.4
dr_rho_pred = 0.75
predictor_1_pred = rbern(100, 0.5)
sr_effects_1_pred = c(-3.3*pred_test[i], 2.3*pred_test[i])

G = simulate_sbm_network(N_id = N_id_pred,     
                         N_groups = N_groups_pred,
                         group_probs = group_probs_pred,                         
                         in_block = in_block_pred,
                         out_block = out_block_pred,
                         sr_mu = sr_mu_pred,  
                         sr_sigma = sr_sigma_pred, 
                         sr_rho = sr_rho_pred,
                         dr_mu = dr_mu_pred,  
                         dr_sigma = dr_sigma_pred, 
                         dr_rho = dr_rho_pred,
                         individual_predictors = data.frame(Status=predictor_1_pred),
                         dyadic_predictors = NULL,
                         individual_effects = matrix(sr_effects_1_pred,nrow=2,ncol=1),
                         dyadic_effects = NULL
                         )        

model_dat = make_strand_data(self_report=list(G$network),  group_ids=factor(G$group_ids), individual_covariates=data.frame(Status=predictor_1_pred), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_pred <- mclapply(1:PP, function(z){

  fit_block_model(data=dat_list[[z]],
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

stanfit_pred <- vector("list", PP)
for(i in 1:PP)
stanfit_pred[[i]] <- rstan::read_stan_csv(fit_pred[[i]]$fit$output_files())


save.image(file = "pred.RData")










