#####################################################################################################################
parameter_recovery_points_plots <- function(PP, fit, Input,  recov_res_true_m, ID){
 recov_res_latent_m <- recov_res_latent_l <- recov_res_latent_h <- matrix(NA, nrow=PP, ncol=15)
 
# First extract sample from Stan    
for(k in 1:PP){
 sr_L = extract(fit[[k]],pars="sr_L")$sr_L
 dr_L = extract(fit[[k]],pars="dr_L")$dr_L
 sr_sigma = extract(fit[[k]],pars="sr_sigma")$sr_sigma
 dr_sigma = extract(fit[[k]],pars="dr_sigma")$dr_sigma

 false_positive_rate = extract(fit[[k]],pars="false_positive_rate")$false_positive_rate
 recall_of_true_ties = extract(fit[[k]],pars="recall_of_true_ties")$recall_of_true_ties

 fpr_sigma = extract(fit[[k]],pars="fpr_sigma")$fpr_sigma
 rtt_sigma = extract(fit[[k]],pars="rtt_sigma")$rtt_sigma
 fpr_raw = extract(fit[[k]],pars="fpr_raw")$fpr_raw
 rtt_raw = extract(fit[[k]],pars="rtt_raw")$rtt_raw

 B = extract(fit[[k]],pars="B")$B

# Now organize estimates of interest
 latent_res <- matrix(nrow=dim(fpr_raw)[1], ncol=15)  

 for( q in 1:dim(fpr_raw)[1]){
  latent_res[q,1] <- logistic(B[q,1,1])
  latent_res[q,2] <- logistic(B[q,1,2])

  latent_res[q,3] <- dr_sigma[q]
  latent_res[q,4] <- dr_L[q,2,1]

  latent_res[q,5] <- sr_sigma[q,1]
  latent_res[q,6] <- sr_sigma[q,2]
  latent_res[q,7] <- sr_L[q,2,1]

  latent_res[q,8] <- false_positive_rate[q,1]
  latent_res[q,9] <- false_positive_rate[q,2]

  latent_res[q,10] <- recall_of_true_ties[q,1]
  latent_res[q,11] <- recall_of_true_ties[q,2]

  latent_res[q,12] <- fpr_sigma[q,1]
  latent_res[q,13] <- fpr_sigma[q,2]

  latent_res[q,14] <- rtt_sigma[q,1]
  latent_res[q,15] <- rtt_sigma[q,2]
  }
                                              
  x <- apply(latent_res, 2, HPDI)

  recov_res_latent_m[k,] <- apply(latent_res, 2, mean)
  recov_res_latent_l[k,] <- x[1,] 
  recov_res_latent_h[k,] <- x[2,]
 }

Estimate <- c("In block", "Out Block", "Dyadic sigma", "Dyadic rho", 
               "Focal sigma", "Target sigma", "Focal-Target rho",
               "Outgoing FPR", "Incoming FPR", "Outgoing RTT", "Incoming RTT", 
               "Outgoing FPR sigma", "Incoming FPR sigma", "Outgoing RTT sigma", "Incoming RTT sigma")

EstimateOrder <- c("Outgoing FPR", "Incoming FPR", "Outgoing RTT", "Incoming RTT", 
               "Outgoing FPR sigma", "Incoming FPR sigma", "Outgoing RTT sigma", "Incoming RTT sigma",
               "In block", "Out Block", "Dyadic rho",  "Dyadic sigma", 
               "Focal-Target rho", "Focal sigma", "Target sigma")

 mp <- vector("list",15)
 for(k in 1:15)
 mp[[k]] <- data.frame(Value=c(recov_res_true_m[,k], recov_res_latent_m[,k]), 
                  L=c(rep(NA,PP),recov_res_latent_l[,k]),
                  H=c(rep(NA,PP),recov_res_latent_h[,k]), 
                  Input=rep(Input,2),
                  Mode=rep( c("True","Estimated"), each=PP),
                  Estimate=Estimate[k]
                  )

mp <- do.call(rbind, mp) 
cols <- c( "True" = "gray13", "Estimated" = "goldenrod3")
mp$Estimate <- factor(mp$Estimate)

mp$Estimate <- factor(mp$Estimate, EstimateOrder)

if(ID=="MIX" | ID == "FPR"){
  title= "False positive rate"
}
if(ID=="FPV"){
  title= "Dispersion in false positive rate"
}
if(ID=="RTT"){
  title= "True tie recall rate"
}
if(ID=="RTV"){
  title= "Dispersion in true tie recall rate"
}

###############
(p4 <- ggplot(mp, aes(x=Input, y=Value, color=Mode, fill=Mode))+
    geom_line(size=1) +
    geom_ribbon(aes(ymin=L, ymax=H), alpha=0.3, colour = NA) +
    facet_wrap(. ~ Estimate, scales = "free_y", nrow=4) +
    scale_fill_manual(name = "Mode:",values = cols) +
    scale_colour_manual(name = "Mode:",values = cols) + xlab(title) + ylab("Parameter value") + theme(legend.position = c(0.89, 0.1))+
    theme(legend.title=element_text(size=14), 
          legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size = 14))+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91", size = 0.5, linetype = "solid")) +
     theme(plot.margin = unit(c(0.1,0.4,0.1,0.1), "cm"))
    )

ggsave(paste0("ParamRecov_points_",ID,".pdf"), p4, height=8,width=12)
}


n_mut_inf = function(x,y){
 dat = discretize(cbind(x,y))
 res = mutinformation(dat[,1], dat[,2]) / sqrt(entropy(dat[,1]) * entropy(dat[,2]))
 return(res)
}


parameter_recovery_corrs_plots <- function(PP, fit, Input, G_list, ID){
####### Parameter-level properties for correlations of effects
param_res_latent_m <- param_res_latent_l <- param_res_latent_h <- matrix(NA,nrow=PP,ncol=7)
    
for(k in 1:PP){
 sr = extract(fit[[k]],pars="sr")$sr
 dr = extract(fit[[k]],pars="dr")$dr

 fpr = extract(fit[[k]],pars="fpr")$fpr
 rtt = extract(fit[[k]],pars="rtt")$rtt

 latent_res <- matrix(nrow=dim(dr)[1], ncol=7)  

  for( q in 1:dim(dr)[1]){
    scrap = dr[q,,]
    diag(scrap) = NA

    latent_res[q,1] <- cor(sr[q,,1], G_list[[k]]$sr[,1])
    latent_res[q,2] <- cor(sr[q,,2], G_list[[k]]$sr[,2])
    latent_res[q,3] <- cor(c(scrap), c(G_list[[k]]$dr),use="pairwise.complete")
    latent_res[q,4] <- cor(fpr[q,,1], G_list[[k]]$fpr[,1])
    latent_res[q,5] <- cor(fpr[q,,2], G_list[[k]]$fpr[,2])
    latent_res[q,6] <- cor(rtt[q,,1], G_list[[k]]$rtt[,1])
    latent_res[q,7] <- cor(rtt[q,,2], G_list[[k]]$rtt[,2])
  }

  x <- apply(latent_res, 2, HPDI)

  param_res_latent_m[k,] <- apply(latent_res, 2, mean)
  param_res_latent_l[k,] <- x[1,] 
  param_res_latent_h[k,] <- x[2,]

   print(paste("Set",k, "of", PP))
 }

 Estimate <- c("Focal offsets", "Target offsets", "Dyadic offsets", "Outgoing FPR", "Incoming FPR", "Outgoing RTT", "Incoming RTT")

 mp <- vector("list",7)
 for(k in 1:7)
 mp[[k]] <- data.frame(Correlation=c(param_res_latent_m[,k]), 
                  L=c(param_res_latent_l[,k]),
                  H=c(param_res_latent_h[,k]), 
                  Input=rep(Input,1),
                  Mode=rep(c("Estimated"), each=PP),
                  Estimate=Estimate[k]
                  )

mp <- do.call(rbind, mp) 
cols <- c("Reported" = "turquoise4", "True" = "gray13", "Estimated" = "goldenrod3")
mp$Estimate <- factor(mp$Estimate)
mp$Estimate <- factor(mp$Estimate, Estimate)

###############
if(ID=="MIX" | ID == "FPR"){
  title= "False positive rate"
}
if(ID=="FPV"){
  title= "Dispersion in false positive rate"
}
if(ID=="RTT"){
  title= "True tie recall rate"
}
if(ID=="RTV"){
  title= "Dispersion in true tie recall rate"
}

(p3 <- ggplot(mp, aes(Input, Correlation, color=Mode, fill=Mode))+
    geom_line(size=1) +
    geom_ribbon(aes(ymin=L, ymax=H), alpha=0.3, colour = NA) +
    facet_wrap(. ~ Estimate,  nrow=1) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) + xlab(title) + theme(legend.position="none") + 
    theme(legend.title=element_text(size=14), 
          legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size = 14))+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
    theme(panel.background = element_rect(fill = "gray91",
                                colour = "gray91",
                                size = 0.5, linetype = "solid"))
    )

ggsave(paste0("ParamRecov_corrs_",ID,".pdf"),p3, height=3,width=15)
}


parameter_recovery_MI_plots <- function(PP, fit, Input, G_list, ID){
####### Parameter-level properties for correlations of effects
param_res_latent_m <- param_res_latent_l <- param_res_latent_h <- matrix(NA,nrow=PP,ncol=7)
    
for(k in 1:PP){
 sr = extract(fit[[k]],pars="sr")$sr
 dr = extract(fit[[k]],pars="dr")$dr

 fpr = extract(fit[[k]],pars="fpr")$fpr
 rtt = extract(fit[[k]],pars="rtt")$rtt

 latent_res <- matrix(nrow=dim(dr)[1], ncol=7)  

  for( q in 1:dim(dr)[1]){
    scrap = dr[q,,]
    diag(scrap) = NA

    latent_res[q,1] <- n_mut_inf(sr[q,,1], G_list[[k]]$sr[,1])
    latent_res[q,2] <- n_mut_inf(sr[q,,2], G_list[[k]]$sr[,2])
    latent_res[q,3] <- n_mut_inf(c(scrap), c(G_list[[k]]$dr))
    latent_res[q,4] <- n_mut_inf(fpr[q,,1], G_list[[k]]$fpr[,1])
    latent_res[q,5] <- n_mut_inf(fpr[q,,2], G_list[[k]]$fpr[,2])
    latent_res[q,6] <- n_mut_inf(rtt[q,,1], G_list[[k]]$rtt[,1])
    latent_res[q,7] <- n_mut_inf(rtt[q,,2], G_list[[k]]$rtt[,2])
  }

  x <- apply(latent_res, 2, HPDI)

  param_res_latent_m[k,] <- apply(latent_res, 2, mean)
  param_res_latent_l[k,] <- x[1,] 
  param_res_latent_h[k,] <- x[2,]

   print(paste("Set",k, "of", PP))
 }

 Estimate <- c("Focal offsets", "Target offsets", "Dyadic offsets", "Outgoing FPR", "Incoming FPR", "Outgoing RTT", "Incoming RTT")

 mp <- vector("list",7)
 for(k in 1:7)
 mp[[k]] <- data.frame(MutualInformation=c(param_res_latent_m[,k]), 
                  L=c(param_res_latent_l[,k]),
                  H=c(param_res_latent_h[,k]), 
                  Input=rep(Input,1),
                  Mode=rep(c("Estimated"), each=PP),
                  Estimate=Estimate[k]
                  )

mp <- do.call(rbind, mp) 
cols <- c("Reported" = "turquoise4", "True" = "gray13", "Estimated" = "goldenrod3")
mp$Estimate <- factor(mp$Estimate)
mp$Estimate <- factor(mp$Estimate, Estimate)

###############
if(ID=="MIX" | ID == "FPR"){
  title= "False positive rate"
}
if(ID=="FPV"){
  title= "Dispersion in false positive rate"
}
if(ID=="RTT"){
  title= "True tie recall rate"
}
if(ID=="RTV"){
  title= "Dispersion in true tie recall rate"
}

(p3 <- ggplot(mp, aes(Input, MutualInformation, color=Mode, fill=Mode))+
    geom_line(size=1) +
    geom_ribbon(aes(ymin=L, ymax=H), alpha=0.3, colour = NA) +
    facet_wrap(. ~ Estimate,  nrow=1) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) + xlab(title) + theme(legend.position="none") + 
    theme(legend.title=element_text(size=14), 
          legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size = 14))+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
    theme(panel.background = element_rect(fill = "gray91",
                                colour = "gray91",
                                size = 0.5, linetype = "solid"))
    )

ggsave(paste0("ParamRecov_MI_",ID,".pdf"),p3, height=3,width=15)
}

network_level_characters_plots <- function(PP, fit, Input,  G_list, ID){
####### Network-level properties #######
 net_res_true <- net_res_claimed2 <- net_res_claimed1 <- net_res_latent_m <- net_res_latent_l <- net_res_latent_h <- matrix(NA, nrow=PP, ncol=5)
  
for(i in 1:PP){
  True_Net <- G_list[[i]]$true_network
  diag(True_Net) <- 0

  Claimed_Net2 <- Claimed_Net1 <- G_list[[i]]$reporting_network
  Claimed_Net1 <- ifelse((Claimed_Net1[,,1] + t(Claimed_Net1[,,2]))>0,1,0)
  Claimed_Net2 <- ifelse((Claimed_Net2[,,1] + t(Claimed_Net2[,,2]))==2,1,0)

  net_res_true[i,] <- get_network_level_properties(True_Net)
  net_res_claimed1[i,] <- get_network_level_properties(Claimed_Net1)
  net_res_claimed2[i,] <- get_network_level_properties(Claimed_Net2)


  net_post <- extract(fit[[i]],pars="p_tie_out")$p_tie_out
 
   latent_res <- matrix(nrow=dim(net_post)[1], ncol=5)  
  for( j in 1:dim(net_post)[1]){
    result <- net_post[j,,]

    for(m in 1:nrow(result)){
      for(n in 1:ncol(result)){
      result[m,n] <- rbern(1,result[m,n])
      }
    }
    latent_res[j,] <- get_network_level_properties(result)
  }

  x <- apply(latent_res, 2, HPDI)

  net_res_latent_m[i,] <- apply(latent_res, 2, mean)
  net_res_latent_l[i,] <- x[1,] 
  net_res_latent_h[i,] <- x[2,]

  print(paste("Set",i, "of", PP))
 }


 Estimate <- c("Density", "Reciprocity", "Transitivity", "Betw. Centra.", "Eigen. Centra." )

 mp <- vector("list",5)
 for(k in 1:5)
 mp[[k]] <- data.frame(Value=c(net_res_true[,k], net_res_claimed1[,k], net_res_claimed2[,k], net_res_latent_m[,k]), 
                  L=c(rep(NA,PP), rep(NA,PP), rep(NA,PP), net_res_latent_l[,k]),
                  H=c(rep(NA,PP), rep(NA,PP), rep(NA,PP), net_res_latent_h[,k]), 
                  Input=rep(Input,4),
                  Mode=rep( c("True","Reported (union)","Reported (intersection)","Estimated"), each=PP),
                  Estimate=Estimate[k]
                  )

mp <- do.call(rbind, mp) 
cols <- c("Reported (union)" = "#b54a31", "Reported (intersection)" = "turquoise4", "Estimated" = "goldenrod3", "True" = "gray13")
mp$Estimate <- factor(mp$Estimate)
mp$Estimate <- factor(mp$Estimate, levels(mp$Estimate)[c(2,4,5,1,3)])

if(ID=="MIX" | ID == "FPR"){
  title= "False positive rate"
}
if(ID=="FPV"){
  title= "Dispersion in false positive rate"
}
if(ID=="RTT"){
  title= "True tie recall rate"
}
if(ID=="RTV"){
  title= "Dispersion in true tie recall rate"
}

###############
(p1 <- ggplot(mp, aes(Input, Value, color=Mode, fill=Mode))+
    geom_line(size=1.25) +
    geom_ribbon(aes(ymin=L, ymax=H), alpha=0.3, colour = NA) +
    facet_wrap(. ~ Estimate, scales = "free_y", nrow=1) +
    scale_fill_manual(name = "Mode:",values = cols) +
    scale_colour_manual(name = "Mode:",values = cols) + xlab(title) + theme(legend.position="bottom") + 
    theme(legend.title=element_text(size=14), 
          legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size = 14))+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91", size = 0.5, linetype = "solid")) +
     theme(plot.margin = unit(c(0.1,0.4,0.1,0.1), "cm"))
    )

ggsave(paste0("NetLevel_",ID,".pdf"),p1, height=4,width=15)
}


individual_level_characters_plots <- function(PP, fit, Input,  G_list, ID){
####### Individual-level properties #######
  ind_res_true <- ind_res_claimed2 <- ind_res_claimed1 <- ind_res_latent_m <- ind_res_latent_l <- ind_res_latent_h <- matrix(NA,nrow=PP,ncol=5)
    
  for(i in 1:PP){
  True_Net <- G_list[[i]]$true_network
  diag(True_Net) <- 0

  Claimed_Net2 <- Claimed_Net1 <- G_list[[i]]$reporting_network
  Claimed_Net1 <- ifelse((Claimed_Net1[,,1] + t(Claimed_Net1[,,2]))>0,1,0)
  Claimed_Net2 <- ifelse((Claimed_Net2[,,1] + t(Claimed_Net2[,,2]))==2,1,0)

  ind_res_true[i,] <- get_node_level_properties(True_Net, True_Net)
  ind_res_claimed1[i,] <- get_node_level_properties(True_Net, Claimed_Net1)
  ind_res_claimed2[i,] <- get_node_level_properties(True_Net, Claimed_Net2)

  net_post <- extract(fit[[i]],pars="p_tie_out")$p_tie_out
 
   latent_res <- matrix(nrow=dim(net_post)[1], ncol=5)  
  for( j in 1:dim(net_post)[1]){
    result <- net_post[j,,]

    for(m in 1:nrow(result)){
      for(n in 1:ncol(result)){
      result[m,n] <- rbern(1,result[m,n])
      }
    }
    latent_res[j,] <- get_node_level_properties(True_Net,result)
  }

  x <- apply(latent_res, 2, HPDI)

  ind_res_latent_m[i,] <- apply(latent_res, 2, mean)
  ind_res_latent_l[i,] <- x[1,] 
  ind_res_latent_h[i,] <- x[2,]
   print(paste("Set",i, "of", PP))
 }

 Estimate <- c("In degree", "Out degree", "Betw. Centra.", "Eigen. Centra.", "P.R. Centra.")

 mp <- vector("list",5)
 for(k in 1:5)
 mp[[k]] <- data.frame(Correlation=c(ind_res_true[,k], ind_res_claimed1[,k],ind_res_claimed2[,k], ind_res_latent_m[,k]), 
                  L=c(rep(NA,PP), rep(NA,PP), rep(NA,PP), ind_res_latent_l[,k]),
                  H=c(rep(NA,PP), rep(NA,PP), rep(NA,PP), ind_res_latent_h[,k]), 
                  Input=rep(Input,4),
                  Mode=rep( c("True","Reported (union)","Reported (intersection)","Estimated"), each=PP),
                  Estimate=Estimate[k]
                  )

mp <- do.call(rbind, mp) 
cols <- c("Reported (union)" = "#b54a31", "Reported (intersection)" = "turquoise4", "Estimated" = "goldenrod3", "True" = "gray13")
mp$Estimate <- factor(mp$Estimate)
mp$Estimate <- factor(mp$Estimate, levels(mp$Estimate)[c(3,4,1,2,5)])

mp$Mode = factor(mp$Mode)

###############
if(ID=="MIX" | ID == "FPR"){
  title= "False positive rate"
}
if(ID=="FPV"){
  title= "Dispersion in false positive rate"
}
if(ID=="RTT"){
  title= "True tie recall rate"
}
if(ID=="RTV"){
  title= "Dispersion in true tie recall rate"
}
(p2 <- ggplot(mp, aes(Input, Correlation, color=Mode, fill=Mode))+
    geom_line(size=1.25) +
    geom_ribbon(aes(ymin=L, ymax=H), alpha=0.3, colour = NA) +
    facet_wrap(. ~ Estimate, scales = "free_y", nrow=1) +
    scale_fill_manual(name = "Mode:", values = cols) +
    scale_colour_manual(name = "Mode:", values = cols) + xlab(title) + theme(legend.position="bottom") + 
    theme(legend.title=element_text(size=14), 
          legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size = 14))+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
    theme(panel.background = element_rect(fill = "gray91",
                                colour = "gray91",
                                size = 0.5, linetype = "solid"))

    )

ggsave(paste0("IndLevel_",ID,".pdf"),p2, height=4,width=15)
}


individual_level_MI_characters_plots <- function(PP, fit, Input,  G_list, ID){
####### Individual-level properties #######
  ind_res_true <- ind_res_claimed2 <- ind_res_claimed1 <- ind_res_latent_m <- ind_res_latent_l <- ind_res_latent_h <- matrix(NA,nrow=PP,ncol=5)
    
  for(i in 1:PP){
  True_Net <- G_list[[i]]$true_network
  diag(True_Net) <- 0

  Claimed_Net2 <- Claimed_Net1 <- G_list[[i]]$reporting_network
  Claimed_Net1 <- ifelse((Claimed_Net1[,,1] + t(Claimed_Net1[,,2]))>0,1,0)
  Claimed_Net2 <- ifelse((Claimed_Net2[,,1] + t(Claimed_Net2[,,2]))==2,1,0)

  ind_res_true[i,] <- get_node_level_properties_MI(True_Net, True_Net)
  ind_res_claimed1[i,] <- get_node_level_properties_MI(True_Net, Claimed_Net1)
  ind_res_claimed2[i,] <- get_node_level_properties_MI(True_Net, Claimed_Net2)

  net_post <- extract(fit[[i]],pars="p_tie_out")$p_tie_out
 
   latent_res <- matrix(nrow=dim(net_post)[1], ncol=5)  
  for( j in 1:dim(net_post)[1]){
    result <- net_post[j,,]

    for(m in 1:nrow(result)){
      for(n in 1:ncol(result)){
      result[m,n] <- rbern(1,result[m,n])
      }
    }
    latent_res[j,] <- get_node_level_properties_MI(True_Net,result)
  }

  x <- apply(latent_res, 2, HPDI)

  ind_res_latent_m[i,] <- apply(latent_res, 2, mean)
  ind_res_latent_l[i,] <- x[1,] 
  ind_res_latent_h[i,] <- x[2,]
   print(paste("Set",i, "of", PP))
 }

 Estimate <- c("In degree", "Out degree", "Betw. Centra.", "Eigen. Centra.", "P.R. Centra.")

 mp <- vector("list",5)
 for(k in 1:5)
 mp[[k]] <- data.frame(MutualInformation=c(ind_res_true[,k], ind_res_claimed1[,k],ind_res_claimed2[,k], ind_res_latent_m[,k]), 
                  L=c(rep(NA,PP), rep(NA,PP), rep(NA,PP), ind_res_latent_l[,k]),
                  H=c(rep(NA,PP), rep(NA,PP), rep(NA,PP), ind_res_latent_h[,k]), 
                  Input=rep(Input,4),
                  Mode=rep( c("True","Reported (union)","Reported (intersection)","Estimated"), each=PP),
                  Estimate=Estimate[k]
                  )

mp <- do.call(rbind, mp) 
cols <- c("Reported (union)" = "#b54a31", "Reported (intersection)" = "turquoise4", "Estimated" = "goldenrod3", "True" = "gray13")
mp$Estimate <- factor(mp$Estimate)
mp$Estimate <- factor(mp$Estimate, levels(mp$Estimate)[c(3,4,1,2,5)])

mp$Mode = factor(mp$Mode)

###############
if(ID=="MIX" | ID == "FPR"){
  title= "False positive rate"
}
if(ID=="FPV"){
  title= "Dispersion in false positive rate"
}
if(ID=="RTT"){
  title= "True tie recall rate"
}
if(ID=="RTV"){
  title= "Dispersion in true tie recall rate"
}
(p2 <- ggplot(mp, aes(Input, MutualInformation, color=Mode, fill=Mode))+
    geom_line(size=1.25) +
    geom_ribbon(aes(ymin=L, ymax=H), alpha=0.3, colour = NA) +
    facet_wrap(. ~ Estimate, scales = "free_y", nrow=1) +
    scale_fill_manual(name = "Mode:", values = cols) +
    scale_colour_manual(name = "Mode:", values = cols) + xlab(title) + theme(legend.position="bottom") + 
    theme(legend.title=element_text(size=14), 
          legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size = 14))+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
    theme(panel.background = element_rect(fill = "gray91",
                                colour = "gray91",
                                size = 0.5, linetype = "solid"))

    )

ggsave(paste0("IndLevel_MI_",ID,".pdf"),p2, height=4,width=15)
}



#####################################################################################################################
F1_plots <- function(PP, fit, Input, G_list, ID){
 recov_res_old_i <- recov_res_old_u <- recov_res_latent_l <- recov_res_latent_h  <- recov_res_latent_m <- matrix(NA, nrow=PP, ncol=3)

for(k in 1:PP){
   Q_real = G_list[[k]]$true_network
   diag(Q_real) = 0

   Claimed_Net1 = G_list[[k]]$reporting_network[,,1]
   Claimed_Net2 = G_list[[k]]$reporting_network[,,2]

   diag(Claimed_Net1) = 0
   diag(Claimed_Net2) = 0

   Claimed_Net_U = ifelse((Claimed_Net1 + t(Claimed_Net2))>0,1,0)
   Claimed_Net_I = ifelse((Claimed_Net1 + t(Claimed_Net2))==2,1,0)

   ##########
   Q_true = c(ifelse(Q_real==1,"tie","no_tie"))
   Z = table(c(Claimed_Net_U),c(Q_true))

   tp = Z[2,2]
   tn = Z[1,1]

   fp = Z[2,1]
   fn = Z[1,2]


   precision_u = tp / (tp + fp)
   recall_u = tp / (tp + fn)
   F1_u = tp / (tp + 0.5*(fp + fn))

   ##########
   Q_true = c(ifelse(Q_real==1,"tie","no_tie"))
   Z = table(c(Claimed_Net_I),c(Q_true))

   tp = Z[2,2]
   tn = Z[1,1]

   fp = Z[2,1]
   fn = Z[1,2]


   precision_i = tp / (tp + fp)
   recall_i = tp / (tp + fn)
   F1_i = tp / (tp + 0.5*(fp + fn))

   recov_res_old_u[k,] = c(precision_u, recall_u, F1_u)
   recov_res_old_i[k,] = c(precision_i, recall_i, F1_i)

} 



Estimate <- c("Precision", "Recall", "F1 score")

 mp1 <- vector("list",3)
 for(k in 1:3)
 mp1[[k]] <- data.frame(Value=c(recov_res_old_u[,k]), 
                  L=c(rep(NA,PP)),
                  H=c(rep(NA,PP)),
                  Input=Input,
                  Estimate=Estimate[k],
                  Mode="Reported (union)"
                  )

mp1 <- do.call(rbind, mp1) 

 mp2 <- vector("list",3)
 for(k in 1:3)
 mp2[[k]] <- data.frame(Value=c(recov_res_old_i[,k]), 
                  L=c(rep(NA,PP)),
                  H=c(rep(NA,PP)),
                  Input=Input,
                  Estimate=Estimate[k],
                  Mode="Reported (intersection)"
                  )

mp2 <- do.call(rbind, mp2) 


 
# First extract sample from Stan    
for(k in 1:PP){
 Q_real = G_list[[k]]$true_network
 diag(Q_real) = 0
 Q = rstan::extract(fit[[k]], pars="p_tie_out")$p_tie_out

# Now organize estimates of interest
 latent_res <- matrix(nrow=dim(Q)[1], ncol=3)  

 for( q in 1:dim(Q)[1]){

 #####
  Q_true = c(ifelse(Q_real==1,"tie","no_tie"))
  Z = table(c(round(Q[q,,],0)),c(Q_true))

  tp = Z[2,2]
  tn = Z[1,1]

  fp = Z[2,1]
  fn = Z[1,2]


  precision = tp / (tp + fp)
  recall = tp / (tp + fn)
  F1 = tp / (tp + 0.5*(fp + fn))

  latent_res[q,1] <- precision
  latent_res[q,2] <- recall
  latent_res[q,3] <- F1
  }
                                              
  x <- apply(latent_res, 2, HPDI)

  recov_res_latent_m[k,] <- apply(latent_res, 2, mean)
  recov_res_latent_l[k,] <- x[1,] 
  recov_res_latent_h[k,] <- x[2,]
 }

Estimate <- c("Precision", "Recall", "F1 score")

 mp <- vector("list",3)
 for(k in 1:3)
 mp[[k]] <- data.frame(Value=c(recov_res_latent_m[,k]), 
                  L=c(recov_res_latent_l[,k]),
                  H=c(recov_res_latent_h[,k]), 
                  Input=Input,
                  Estimate=Estimate[k],
                  Mode = "Estimated"
                  )

mp <- do.call(rbind, mp) 


mp <- rbind(mp,rbind(mp1,mp2))

cols <- c("Reported (union)" = "#b54a31", "Reported (intersection)" = "turquoise4", "Estimated" = "goldenrod3")
mp$Estimate <- factor(mp$Estimate)

#mp$Estimate <- factor(mp$Estimate, EstimateOrder)

if(ID=="MIX" | ID == "FPR"){
  title= "False positive rate"
}
if(ID=="FPV"){
  title= "Dispersion in false positive rate"
}
if(ID=="RTT"){
  title= "True tie recall rate"
}
if(ID=="RTV"){
  title= "Dispersion in true tie recall rate"
}

###############
(p4 <- ggplot(mp, aes(x=Input, y=Value, color=Mode, fill=Mode))+
    geom_line(size=1) +
    geom_ribbon(aes(ymin=L, ymax=H), alpha=0.3, colour = NA) +
    facet_wrap(. ~ Estimate, scales = "free_y", nrow=1) +
    scale_fill_manual(name = "Mode:",values = cols) +
    scale_colour_manual(name = "Mode:",values = cols) + xlab(title) + ylab("Performance metric") + theme(legend.position = "bottom")+
    theme(strip.text.x = element_text(size = 14))+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91", size = 0.5, linetype = "solid")) +
     theme(plot.margin = unit(c(0.1,0.4,0.1,0.1), "cm")) + ylim(0, 1)
    )

ggsave(paste0("F1_Prec_Rec_",ID,".pdf"), p4, height=4,width=12)
}

