############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=8)

recov_res_true_m[,1] <- rep(B_dr_sigma, PP)
recov_res_true_m[,2] <- DR_sigma_test
recov_res_true_m[,3] <- rep(dr_rho_dr_sigma, PP)
recov_res_true_m[,4] <- rep(sr_sigma_dr_sigma[1], PP)
recov_res_true_m[,5] <- rep(sr_sigma_dr_sigma[2], PP)
recov_res_true_m[,6] <- rep(sr_rho_dr_sigma, PP)
recov_res_true_m[,7] <- rep(0, PP)
recov_res_true_m[,8] <- rep(0, PP)

parameter_recovery_points_plots(PP, stanfit_dr_sigma, DR_sigma_test, recov_res_true_m, "DRS")
parameter_recovery_corrs_plots(PP, stanfit_dr_sigma, DR_sigma_test, G_list, "DRS")
parameter_recovery_MI_plots(PP, stanfit_dr_sigma, DR_sigma_test, G_list, "DRS")


