############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=8)

recov_res_true_m[,1] <- rep(B_sr_predictor, PP)
recov_res_true_m[,2] <- rep(dr_sigma_sr_predictor, PP)
recov_res_true_m[,3] <- rep(dr_rho_sr_predictor, PP)
recov_res_true_m[,4] <- rep(sr_sigma_sr_predictor[1], PP)
recov_res_true_m[,5] <- rep(sr_sigma_sr_predictor[2], PP)
recov_res_true_m[,6] <- rep(sr_rho_sr_predictor, PP)
recov_res_true_m[,7] <- -0.3*SR_predictor_test
recov_res_true_m[,8] <- SR_predictor_test

parameter_recovery_points_plots(PP, stanfit_sr_predictor, SR_predictor_test, recov_res_true_m, "SRP")
parameter_recovery_corrs_plots(PP, stanfit_sr_predictor, SR_predictor_test, G_list, "SRP")
parameter_recovery_MI_plots(PP, stanfit_sr_predictor, SR_predictor_test, G_list, "SRP")


