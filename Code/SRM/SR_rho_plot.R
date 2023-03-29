############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=8)

recov_res_true_m[,1] <- rep(B_sr_rho, PP)
recov_res_true_m[,2] <- rep(dr_sigma_sr_rho, PP)
recov_res_true_m[,3] <- rep(dr_rho_sr_rho, PP)
recov_res_true_m[,4] <- rep(sr_sigma_sr_rho[1], PP)
recov_res_true_m[,5] <- rep(sr_sigma_sr_rho[2], PP)
recov_res_true_m[,6] <- SR_rho_test
recov_res_true_m[,7] <- rep(0, PP)
recov_res_true_m[,8] <- rep(0, PP)

parameter_recovery_points_plots(PP, stanfit_sr_rho, SR_rho_test, recov_res_true_m, "SRR")
parameter_recovery_corrs_plots(PP, stanfit_sr_rho, SR_rho_test, G_list, "SRR")
parameter_recovery_MI_plots(PP, stanfit_sr_rho, SR_rho_test, G_list, "SRR")

