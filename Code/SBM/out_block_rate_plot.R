############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=4)

recov_res_true_m[,1] <- rep(in_block_out_block, PP) 
recov_res_true_m[,2] <- out_block_rate_test
recov_res_true_m[,3] <- rep(0, PP)
recov_res_true_m[,4] <- rep(0, PP)

parameter_recovery_points_plots(PP, stanfit_out_block, out_block_rate_test, recov_res_true_m, "OBR")


