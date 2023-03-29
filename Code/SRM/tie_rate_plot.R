############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=8)

recov_res_true_m[,1] <- in_block_rate_test
recov_res_true_m[,2] <- rep(dr_sigma_in_block, PP)
recov_res_true_m[,3] <- rep(dr_rho_in_block, PP)
recov_res_true_m[,4] <- rep(sr_sigma_in_block[1], PP)
recov_res_true_m[,5] <- rep(sr_sigma_in_block[2], PP)
recov_res_true_m[,6] <- rep(sr_rho_in_block, PP)
recov_res_true_m[,7] <- rep(0, PP)
recov_res_true_m[,8] <- rep(0, PP)

parameter_recovery_points_plots(PP, stanfit_in_block, in_block_rate_test, recov_res_true_m, "IBR")
parameter_recovery_corrs_plots(PP, stanfit_in_block, in_block_rate_test, G_list, "IBR")
parameter_recovery_MI_plots(PP, stanfit_in_block, in_block_rate_test, G_list, "IBR")

