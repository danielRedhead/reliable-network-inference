############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=9)

recov_res_true_m[,1] <- rep(in_block_dr_rho, PP)
recov_res_true_m[,2] <- rep(out_block_dr_rho, PP)
recov_res_true_m[,3] <- rep(dr_sigma_dr_rho, PP)
recov_res_true_m[,4] <- DR_rho_test
recov_res_true_m[,5] <- rep(sr_sigma_dr_rho[1], PP)
recov_res_true_m[,6] <- rep(sr_sigma_dr_rho[2], PP)
recov_res_true_m[,7] <- rep(sr_rho_dr_rho, PP)
recov_res_true_m[,8] <- rep(0, PP)
recov_res_true_m[,9] <- rep(0, PP)

parameter_recovery_points_plots(PP, stanfit_dr_rho, DR_rho_test, recov_res_true_m, "DRR")
parameter_recovery_corrs_plots(PP, stanfit_dr_rho, DR_rho_test, G_list, "DRR")
parameter_recovery_MI_plots(PP, stanfit_dr_rho, DR_rho_test, G_list, "DRR")

