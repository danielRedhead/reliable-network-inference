############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=17)

recov_res_true_m[,1] <- rep(in_block_mix, PP)
recov_res_true_m[,2] <- rep(out_block_mix, PP)
recov_res_true_m[,3] <- rep(dr_sigma_mix, PP)
recov_res_true_m[,4] <- rep(dr_rho_mix, PP)
recov_res_true_m[,5] <- rep(sr_sigma_mix[1], PP)
recov_res_true_m[,6] <- rep(sr_sigma_mix[2], PP)
recov_res_true_m[,7] <- rep(sr_rho_mix, PP)
recov_res_true_m[,8] <- FPR_test
recov_res_true_m[,9] <- FPR_test
recov_res_true_m[,10] <- RTT_test
recov_res_true_m[,11] <- RTT_test
recov_res_true_m[,12] <- rep(fpr_sigma_mix[1], PP)
recov_res_true_m[,13] <- rep(fpr_sigma_mix[2], PP)
recov_res_true_m[,14] <- rep(rtt_sigma_mix[1], PP)
recov_res_true_m[,15] <- rep(rtt_sigma_mix[2], PP)
recov_res_true_m[,16] <- rep(theta_mean_mix, PP)
recov_res_true_m[,17] <- rep(theta_sigma_mix, PP)

parameter_recovery_points_plots(PP, stanfit_mix, FPR_test, recov_res_true_m, "MIX")
parameter_recovery_corrs_plots(PP, stanfit_mix, FPR_test, G_list, "MIX")
parameter_recovery_MI_plots(PP, stanfit_mix, FPR_test, G_list, "MIX")
network_level_characters_plots(PP, stanfit_mix, FPR_test,  G_list, "MIX")
individual_level_characters_plots(PP, stanfit_mix, FPR_test,  G_list, "MIX")
individual_level_MI_characters_plots(PP, stanfit_mix, FPR_test,  G_list, "MIX")

recov_slope_true_m <- matrix(NA,nrow=PP,ncol=6)

recov_slope_true_m[,1] <- rep(fpr_effects_1_mix[1], PP)
recov_slope_true_m[,2] <- rep(fpr_effects_1_mix[2], PP)
recov_slope_true_m[,3] <- rep(rtt_effects_1_mix[1], PP)
recov_slope_true_m[,4] <- rep(rtt_effects_1_mix[2], PP)
recov_slope_true_m[,5] <- rep(sr_effects_1_mix[1], PP)
recov_slope_true_m[,6] <- rep(sr_effects_1_mix[2], PP)

parameter_recovery_slopes_plots(PP, stanfit_mix, FPR_test, recov_slope_true_m, "MIX")
