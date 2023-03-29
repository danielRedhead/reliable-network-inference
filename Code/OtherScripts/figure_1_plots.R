set.seed(1)
####### Figure 1
N_id_fpr = 100
N_groups_fpr = 3
gprobs_fpr = c(0.2, 0.5, 0.3)
in_block_fpr = 0.021
out_block_fpr = 0.0007
sr_mu_fpr = c(0,0)  
sr_sigma_fpr = c(1.7, 0.7) 
sr_rho_fpr = 0.5
dr_mu_fpr = c(0,0) 
dr_sigma_fpr = 1.4
dr_rho_fpr = 0.75
predictor_1_fpr = rbern(100, 0)
sr_effects_1_fpr = c(-0.3, 1.3)
fpr_effects_1_fpr = c(0.1, -0.9, 0.0)
rtt_effects_1_fpr = c(-0.1, 0.3, -0.3)
theta_effects_1_fpr = 0.0
false_positive_rate_fpr = c(0.01, 0.01, 0.01)
recall_of_true_ties_fpr = c(0.7, 0.7, 0.7)
fpr_sigma_fpr = c(0.7, 0.7, 0.7) 
rtt_sigma_fpr = c(0.7, 0.7, 0.7)
theta_mean_fpr = 0 
theta_sigma_fpr = 0 
N_responses_fpr = 2
N_periods_fpr = 2
decay_curve_fpr = 0.0*exp(-seq(1,4, length.out=2))

G = simulate_selfreport_network(N_id = N_id_fpr, 
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
                               flow_rate = 0.01
                                 )     


#################################################################################### Basic plots
diag(G$true_network) <- 0
true <- graph_from_adjacency_matrix(G$true_network, mode = c("directed"))

false <- G$reporting_network

false_either <- ifelse(false[,,1]==1 | t(false[,,2])==1, 1,0)
false_both <- ifelse(false[,,1]==1 & t(false[,,2])==1, 1,0)

false_either <- graph_from_adjacency_matrix(false_either, mode = c("directed"))
false_both <- graph_from_adjacency_matrix(false_both, mode = c("directed"))

coords <- layout_with_kk(true)

V(false_both)$color <- V(false_either)$color <- V(true)$color <- c("turquoise4","gray13", "goldenrod3")[G$group_id]
V(false_both)$frame.color <- V(false_either)$frame.color <- V(true)$frame.color <- c("turquoise4","gray13", "goldenrod3")[G$group_id]

deg <- degree(true, mode="all")
degf <- degree(false_either, mode="all")
degb <- degree(false_both, mode="all")

figure1_t <- plot(true, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = ((deg+4)*0.25), #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)
windows()
figure1_e <- plot(false_either, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = ((degf+4)*0.25), #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)

windows()
figure1_b <- plot(false_both, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = ((degb+4)*0.25), #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)

# plot true out network
Cairo(file="true.png", 
      type="png",
      units="in", 
      width=12, 
      height=12, 
      pointsize=36, 
      dpi=150)

par(mar=c(1,1,1,1)+0)
plot(true, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = ((deg+4)*0.25), #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)
dev.off()

# plot true out network
Cairo(file="intersection.png", 
      type="png",
      units="in", 
      width=12, 
      height=12, 
      pointsize=36, 
      dpi=150)

par(mar=c(1,1,1,1)+0)
plot(false_both, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = ((degb+4)*0.25), #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)
dev.off()


# plot true out network
Cairo(file="union.png", 
      type="png",
      units="in", 
      width=12, 
      height=12, 
      pointsize=36, 
      dpi=150)

par(mar=c(1,1,1,1)+0)
plot(false_either, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = ((degf+4)*0.25), #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)
dev.off()

######################################################################## Latent plot
# Prepare for Stan
group_ids = factor(G$group_ids)
net = list(TransferOut=G$reporting_network[,,1], TransferIn=G$reporting_network[,,2])
ind = data.frame(Status=predictor_1_fpr)

model_dat = make_strand_data(self_report=net, group_ids=group_ids, individual_covariates=ind, dyadic_covariates=NULL)


fit1 = fit_latent_network_model(data=model_dat,
                                focal_regression = ~ 1,
                                target_regression = ~ 1,
                                dyad_regression = ~ 1,
                                fpr_regression = ~ 1,
                                rtt_regression = ~ 1,
                                theta_regression = ~ 1,
                                mode="mcmc",
                                return_latent_network = TRUE,
                                stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 1000,
                                iter_sampling = 1000, max_treedepth = NULL, adapt_delta = NULL)
                                              )

stanfit <- rstan::read_stan_csv(fit1$fit$output_files())

net_post <- extract(stanfit,pars="p_tie_out")$p_tie_out

p_tie = apply(net_post,2:3,median)

latent <- graph_from_adjacency_matrix(round(p_tie), mode = c("directed"))

coords <- layout_with_kk(true)
V(latent)$color <- c("turquoise4","gray13", "goldenrod3")[G$group_id]
V(latent)$frame.color <- c("turquoise4","gray13", "goldenrod3")[G$group_id]

degl <- degree(latent, mode="all")


figure1_l <- plot(latent, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = ((degl+4)*0.25), #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)

# plot true out network
Cairo(file="latent.png", 
      type="png",
      units="in", 
      width=12, 
      height=12, 
      pointsize=36, 
      dpi=150)

par(mar=c(1,1,1,1)+0)
plot(latent, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = ((degl+4)*0.25), #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)
dev.off()

