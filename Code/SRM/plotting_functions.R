#####################################################################################################################
parameter_recovery_points_plots <- function(PP, fit, Input,  recov_res_true_m, ID){
 recov_res_latent_m <- recov_res_latent_l <- recov_res_latent_h <- matrix(NA, nrow=PP, ncol=8)
 
# First extract sample from Stan    
for(k in 1:PP){
 sr_L = extract(fit[[k]],pars="sr_L")$sr_L
 dr_L = extract(fit[[k]],pars="dr_L")$dr_L
 sr_sigma = extract(fit[[k]],pars="sr_sigma")$sr_sigma
 dr_sigma = extract(fit[[k]],pars="dr_sigma")$dr_sigma

 B = extract(fit[[k]],pars="B")$B

# Now organize estimates of interest
 latent_res <- matrix(nrow=dim(sr_sigma)[1], ncol=8)  

 for( q in 1:dim(sr_sigma)[1]){
  latent_res[q,1] <- logistic(B[q,1,1])

  latent_res[q,2] <- dr_sigma[q]
  latent_res[q,3] <- dr_L[q,2,1]

  latent_res[q,4] <- sr_sigma[q,1]
  latent_res[q,5] <- sr_sigma[q,2]
  latent_res[q,6] <- sr_L[q,2,1]

  latent_res[q,7] <- 0
  latent_res[q,8] <- 0

  if("focal_effects" %in%  fit[[k]]@model_pars){
     focal_effects = extract(fit[[k]],pars="focal_effects")$focal_effects
   latent_res[q,7] <- focal_effects[q,1]
   } 

  if("focal_effects" %in%  fit[[k]]@model_pars){
     target_effects = extract(fit[[k]],pars="target_effects")$target_effects
   latent_res[q,8] <- target_effects[q,1]
   }

  }
                                              
  x <- apply(latent_res, 2, HPDI)

  recov_res_latent_m[k,] <- apply(latent_res, 2, mean)
  recov_res_latent_l[k,] <- x[1,] 
  recov_res_latent_h[k,] <- x[2,]
 }

Estimate <- c("Intercept", "Dyadic sigma", "Dyadic rho", 
               "Focal sigma", "Target sigma", "Focal-Target rho", "Focal predictor", "Target predictor")

EstimateOrder <- c("Intercept", "Dyadic rho",  "Dyadic sigma", 
               "Focal-Target rho", "Focal sigma", "Target sigma", "Focal predictor", "Target predictor")

 mp <- vector("list",8)
 for(k in 1:8)
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

if(ID=="DRR"){
  title= "Dyadic reciprocity"
}
if(ID=="SRR"){
  title= "Generalized reciprocity"
}
if(ID=="DRS"){
  title= "Dyadic sigma"
}
if(ID=="SRS"){
  title= "Sender receiver sigma"
}
if(ID=="SRP"){
  title= "Sender receiver predictor"
}
if(ID=="IBR"){
  title= "Within block rate"
}
if(ID=="OBR"){
  title= "Cross block rate"
}


      

###############
(p4 <- ggplot(mp, aes(x=Input, y=Value, color=Mode, fill=Mode))+
    geom_line(size=1) +
    geom_ribbon(aes(ymin=L, ymax=H), alpha=0.3, colour = NA) +
    facet_wrap(. ~ Estimate, scales = "free_y", nrow=4) +
    scale_fill_manual(name = "Mode:",values = cols) +
    scale_colour_manual(name = "Mode:",values = cols) + xlab(title) + ylab("Parameter value") + theme(legend.position = "bottom") +
    theme(legend.title=element_text(size=14), 
          legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size = 14))+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91", size = 0.5, linetype = "solid")) +
     theme(plot.margin = unit(c(0.1,0.4,0.1,0.1), "cm"))
    )

ggsave(paste0("ParamRecov_points_",ID,".pdf"), p4, height=8,width=8)
}




parameter_recovery_corrs_plots <- function(PP, fit, Input, G_list, ID){
####### Parameter-level properties for correlations of effects
param_res_latent_m <- param_res_latent_l <- param_res_latent_h <- matrix(NA,nrow=PP,ncol=3)
    
for(k in 1:PP){
 sr = extract(fit[[k]],pars="sr")$sr
 dr = extract(fit[[k]],pars="dr")$dr
 
 latent_res <- matrix(nrow=dim(dr)[1], ncol=3)  

  for( q in 1:dim(dr)[1]){
    scrap = dr[q,,]
    diag(scrap) = NA

    latent_res[q,1] <- cor(sr[q,,1], G_list[[k]]$sr[,1])
    latent_res[q,2] <- cor(sr[q,,2], G_list[[k]]$sr[,2])
    latent_res[q,3] <- cor(c(scrap), c(G_list[[k]]$dr),use="pairwise.complete")
  }

  x <- apply(latent_res, 2, HPDI)

  param_res_latent_m[k,] <- apply(latent_res, 2, mean)
  param_res_latent_l[k,] <- x[1,] 
  param_res_latent_h[k,] <- x[2,]

   print(paste("Set",k, "of", PP))
 }

 Estimate <- c("Focal offsets", "Target offsets", "Dyadic offsets")

 mp <- vector("list",3)
 for(k in 1:3)
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
if(ID=="DRR"){
  title= "Dyadic reciprocity"
}
if(ID=="SRR"){
  title= "Generalized reciprocity"
}
if(ID=="DRS"){
  title= "Dyadic sigma"
}
if(ID=="SRS"){
  title= "Sender receiver sigma"
}
if(ID=="SRP"){
  title= "Sender receiver predictor"
}
if(ID=="IBR"){
  title= "Within block rate"
}
if(ID=="OBR"){
  title= "Cross block rate"
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

ggsave(paste0("ParamRecov_corrs_",ID,".pdf"),p3, height=3,width=9)
}


n_mut_inf = function(x,y){
 dat = discretize(cbind(x,y))
 res = mutinformation(dat[,1], dat[,2]) / sqrt(entropy(dat[,1]) * entropy(dat[,2]))
 return(res)
}


parameter_recovery_MI_plots <- function(PP, fit, Input, G_list, ID){
####### Parameter-level properties for correlations of effects
param_res_latent_m <- param_res_latent_l <- param_res_latent_h <- matrix(NA,nrow=PP,ncol=3)
    
for(k in 1:PP){
 sr = extract(fit[[k]],pars="sr")$sr
 dr = extract(fit[[k]],pars="dr")$dr
 
 latent_res <- matrix(nrow=dim(dr)[1], ncol=3)  

  for( q in 1:dim(dr)[1]){
    scrap = dr[q,,]
    diag(scrap) = NA

    latent_res[q,1] <- n_mut_inf(sr[q,,1], G_list[[k]]$sr[,1])
    latent_res[q,2] <- n_mut_inf(sr[q,,2], G_list[[k]]$sr[,2])
    latent_res[q,3] <- n_mut_inf(c(scrap), c(G_list[[k]]$dr))
  }

  x <- apply(latent_res, 2, HPDI)

  param_res_latent_m[k,] <- apply(latent_res, 2, mean)
  param_res_latent_l[k,] <- x[1,] 
  param_res_latent_h[k,] <- x[2,]

   print(paste("Set",k, "of", PP))
 }

 Estimate <- c("Focal offsets", "Target offsets", "Dyadic offsets")

 mp <- vector("list",3)
 for(k in 1:3)
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
if(ID=="DRR"){
  title= "Dyadic reciprocity"
}
if(ID=="SRR"){
  title= "Generalized reciprocity"
}
if(ID=="DRS"){
  title= "Dyadic sigma"
}
if(ID=="SRS"){
  title= "Sender receiver sigma"
}
if(ID=="SRP"){
  title= "Sender receiver predictor"
}
if(ID=="IBR"){
  title= "Within block rate"
}
if(ID=="OBR"){
  title= "Cross block rate"
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

ggsave(paste0("ParamRecov_MI_",ID,".pdf"),p3, height=3,width=9)
}

