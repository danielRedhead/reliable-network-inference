############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=9)

recov_res_true_m[,1] <- rep(in_block_sr_sigma, PP)
recov_res_true_m[,2] <- rep(out_block_sr_sigma, PP)
recov_res_true_m[,3] <- rep(dr_sigma_sr_sigma, PP)
recov_res_true_m[,4] <- rep(dr_rho_sr_sigma, PP)
recov_res_true_m[,5] <- 1.7*SR_sigma_test
recov_res_true_m[,6] <- SR_sigma_test
recov_res_true_m[,7] <- rep(sr_rho_sr_sigma, PP)
recov_res_true_m[,8] <- rep(0, PP)
recov_res_true_m[,9] <- rep(0, PP)

parameter_recovery_points_plots(PP, stanfit_sr_sigma, SR_sigma_test, recov_res_true_m, "SRS")
parameter_recovery_corrs_plots(PP, stanfit_sr_sigma, SR_sigma_test, G_list, "SRS")
parameter_recovery_MI_plots(PP, stanfit_sr_sigma, SR_sigma_test, G_list, "SRS")
