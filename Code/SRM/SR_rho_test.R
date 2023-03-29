####### FALSE POSITIVE RATE #######
PP = 20
SR_rho_test = seq(0.01, 0.8, length.out=PP)
G_list = dat_list = vector("list",PP)

for(i in 1:PP){
  set.seed(1)

####### FPR test pars
N_id_sr_rho = 100
B_sr_rho = 0.02
sr_mu_sr_rho = c(0,0)  
sr_sigma_sr_rho = c(2.1, 0.7) 
sr_rho_sr_rho = SR_rho_test[i]
dr_mu_sr_rho = c(0,0) 
dr_sigma_sr_rho = 1.4
dr_rho_sr_rho = 0.75
predictor_1_sr_rho = rbern(100, 0)
sr_effects_1_sr_rho = c(-0.3, 1.3)

G = simulate_srm_network(N_id = N_id_sr_rho, 
                         B = B_sr_rho, 
                         sr_mu = sr_mu_sr_rho,  
                         sr_sigma = sr_sigma_sr_rho, 
                         sr_rho = sr_rho_sr_rho,
                         dr_mu = dr_mu_sr_rho,  
                         dr_sigma = dr_sigma_sr_rho, 
                         dr_rho = dr_rho_sr_rho,
                         individual_predictors = data.frame(Status=predictor_1_sr_rho),
                         dyadic_predictors = NULL,
                         individual_effects = matrix(sr_effects_1_sr_rho,nrow=2,ncol=1),
                         dyadic_effects = NULL
                         )        

model_dat = make_strand_data(self_report=list(G$network),  group_ids=NULL, individual_covariates=data.frame(Status=predictor_1_sr_rho), dyadic_covariates=NULL)

     dat_list[[i]] <- model_dat
     G_list[[i]] <- G
}


fit_sr_rho <- mclapply(1:PP, function(z){

  fit_social_relations_model(data=dat_list[[z]],
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

stanfit_sr_rho <- vector("list", PP)
for(i in 1:PP)
stanfit_sr_rho[[i]] <- rstan::read_stan_csv(fit_sr_rho[[i]]$fit$output_files())


save.image(file = "SR_rho.RData")










