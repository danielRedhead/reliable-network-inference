############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=9)

recov_res_true_m[,1] <- rep(in_block_out_block, PP)
recov_res_true_m[,2] <- out_block_rate_test
recov_res_true_m[,3] <- rep(dr_sigma_out_block, PP)
recov_res_true_m[,4] <- rep(dr_rho_out_block, PP)
recov_res_true_m[,5] <- rep(sr_sigma_out_block[1], PP)
recov_res_true_m[,6] <- rep(sr_sigma_out_block[2], PP)
recov_res_true_m[,7] <- rep(sr_rho_out_block, PP)
recov_res_true_m[,8] <- rep(0, PP)
recov_res_true_m[,9] <- rep(0, PP)

parameter_recovery_points_plots(PP, stanfit_out_block, out_block_rate_test, recov_res_true_m, "OBR")
parameter_recovery_corrs_plots(PP, stanfit_out_block, out_block_rate_test, G_list, "OBR")
parameter_recovery_MI_plots(PP, stanfit_out_block, out_block_rate_test, G_list, "OBR")
