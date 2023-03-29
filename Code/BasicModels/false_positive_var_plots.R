############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=15)

recov_res_true_m[,1] <- rep(in_block_fpv, PP)
recov_res_true_m[,2] <- rep(out_block_fpv, PP)
recov_res_true_m[,3] <- rep(dr_sigma_fpv, PP)
recov_res_true_m[,4] <- rep(dr_rho_fpv, PP)
recov_res_true_m[,5] <- rep(sr_sigma_fpv[1], PP)
recov_res_true_m[,6] <- rep(sr_sigma_fpv[2], PP)
recov_res_true_m[,7] <- rep(sr_rho_fpv, PP)
recov_res_true_m[,8] <- rep(false_positive_rate_fpv[1], PP)
recov_res_true_m[,9] <- rep(false_positive_rate_fpv[2], PP)
recov_res_true_m[,10] <- rep(recall_of_true_ties_fpv[1], PP)
recov_res_true_m[,11] <- rep(recall_of_true_ties_fpv[2], PP)
recov_res_true_m[,12] <- FPV_test
recov_res_true_m[,13] <- FPV_test
recov_res_true_m[,14] <- rep(rtt_sigma_fpv[1], PP)
recov_res_true_m[,15] <- rep(rtt_sigma_fpv[2], PP)

parameter_recovery_points_plots(PP, stanfit_fpv, FPV_test, recov_res_true_m, "FPV")
parameter_recovery_corrs_plots(PP, stanfit_fpv, FPV_test, G_list, "FPV")
parameter_recovery_MI_plots(PP, stanfit_fpv, FPV_test, G_list, "FPV")
network_level_characters_plots(PP, stanfit_fpv, FPV_test,  G_list, "FPV")
individual_level_characters_plots(PP, stanfit_fpv, FPV_test,  G_list, "FPV")
individual_level_MI_characters_plots(PP, stanfit_fpv, FPV_test,  G_list, "FPV")

F1_plots(PP, stanfit_fpv, FPV_test, G_list, "FPV")



