############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=17)
recov_res_true_m[,1] <- rep(in_block_rtv, PP)
recov_res_true_m[,2] <- rep(out_block_rtv, PP)
recov_res_true_m[,3] <- rep(dr_sigma_rtv, PP)
recov_res_true_m[,4] <- rep(dr_rho_rtv, PP)
recov_res_true_m[,5] <- rep(sr_sigma_rtv[1], PP)
recov_res_true_m[,6] <- rep(sr_sigma_rtv[2], PP)
recov_res_true_m[,7] <- rep(sr_rho_rtv, PP)
recov_res_true_m[,8] <- rep(false_positive_rate_rtv[1],PP)
recov_res_true_m[,9] <- rep(false_positive_rate_rtv[2],PP)
recov_res_true_m[,10] <- rep(recall_of_true_ties_rtv[1],PP)
recov_res_true_m[,11] <- rep(recall_of_true_ties_rtv[2],PP)
recov_res_true_m[,12] <- rep(fpr_sigma_rtv[1], PP)
recov_res_true_m[,13] <- rep(fpr_sigma_rtv[2], PP)
recov_res_true_m[,14] <- RTV_test
recov_res_true_m[,15] <- RTV_test
recov_res_true_m[,16] <- rep(theta_mean_rtv, PP)
recov_res_true_m[,17] <- rep(theta_sigma_rtv, PP)

parameter_recovery_points_plots(PP, stanfit_rtv, RTV_test, recov_res_true_m, "RTV")
parameter_recovery_corrs_plots(PP, stanfit_rtv, RTV_test, G_list, "RTV")
parameter_recovery_MI_plots(PP, stanfit_rtv, RTV_test, G_list, "RTV")
network_level_characters_plots(PP, stanfit_rtv, RTV_test,  G_list, "RTV")
individual_level_characters_plots(PP, stanfit_rtv, RTV_test,  G_list, "RTV")
individual_level_MI_characters_plots(PP, stanfit_rtv, RTV_test,  G_list, "RTV")


recov_slope_true_m <- matrix(NA,nrow=PP,ncol=6)

recov_slope_true_m[,1] <- rep(fpr_effects_1_rtv[1], PP)
recov_slope_true_m[,2] <- rep(fpr_effects_1_rtv[2], PP)
recov_slope_true_m[,3] <- rep(rtt_effects_1_rtv[1], PP)
recov_slope_true_m[,4] <- rep(rtt_effects_1_rtv[2], PP)
recov_slope_true_m[,5] <- rep(sr_effects_1_rtv[1], PP)
recov_slope_true_m[,6] <- rep(sr_effects_1_rtv[2], PP)

parameter_recovery_slopes_plots(PP, stanfit_rtv, RTV_test, recov_slope_true_m, "RTV")


parameter_recovery_decay_plots(PP, stanfit_rtv, RTV_test, decay_curve_rtv, "RTV")
parameter_recovery_flow_plots(PP, stanfit_rtv, RTV_test, flow_rate_rtv, "RTV")


