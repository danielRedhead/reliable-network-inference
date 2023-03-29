############################################## Parameter recovery
recov_res_true_m <- matrix(NA,nrow=PP,ncol=4)

recov_res_true_m[,1] <- rep(in_block_pred, PP) 
recov_res_true_m[,2] <- rep(out_block_pred, PP) 
recov_res_true_m[,3] <- -3.3*pred_test
recov_res_true_m[,4] <- 2.3*pred_test

parameter_recovery_points_plots(PP, stanfit_pred, pred_test, recov_res_true_m, "PPP")


