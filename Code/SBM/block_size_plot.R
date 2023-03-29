############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=4)

recov_res_true_m[,1] <- rep(in_block_block_size, PP) 
recov_res_true_m[,2] <- rep(out_block_block_size, PP) 
recov_res_true_m[,3] <- rep(0, PP)
recov_res_true_m[,4] <- rep(0, PP)

parameter_recovery_points_plots(PP, stanfit_block_size, block_prob_rate_test, recov_res_true_m, "BRD")


