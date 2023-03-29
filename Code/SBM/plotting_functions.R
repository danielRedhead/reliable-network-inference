#####################################################################################################################
parameter_recovery_points_plots <- function(PP, fit, Input,  recov_res_true_m, ID, mode){
 recov_res_latent_m <- recov_res_latent_l <- recov_res_latent_h <- matrix(NA, nrow=PP, ncol=4)
 
# First extract sample from Stan    
for(k in 1:PP){
 B = extract(fit[[k]],pars="block_effects")$block_effects

# Now organize estimates of interest
 latent_res <- matrix(nrow=dim(B)[1], ncol=4)  

 for( q in 1:dim(B)[1]){
  latent_res[q,1] <- B[q,1] + B[q,2] + B[q,14]
  latent_res[q,2] <- B[q,1] + B[q,3] + B[q,14]

  latent_res[q,3] <- 0
  latent_res[q,4] <- 0

  if("focal_effects" %in%  fit[[k]]@model_pars){
     focal_effects = extract(fit[[k]],pars="focal_effects")$focal_effects
   latent_res[q,3] <- focal_effects[q,1]
   } 

  if("focal_effects" %in%  fit[[k]]@model_pars){
     target_effects = extract(fit[[k]],pars="target_effects")$target_effects
   latent_res[q,4] <- target_effects[q,1]
   }



  }
                                              
  x <- apply(latent_res, 2, HPDI)

  recov_res_latent_m[k,] <- apply(latent_res, 2, mean)
  recov_res_latent_l[k,] <- x[1,] 
  recov_res_latent_h[k,] <- x[2,]
 }

Estimate <- c("In block", "Out Block", "Focal predictor", "Target predictor")

EstimateOrder <- c("In block", "Out Block", "Focal predictor", "Target predictor")

 mp <- vector("list",4)
 for(k in 1:4)
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


if(ID=="SST"){
  title= "Sample size"
}
if(ID=="BRD"){
  title= "Block relative density"
}
if(ID=="NBR"){
  title= "Number of blocks"
}
if(ID=="IBR"){
  title= "Within block rate"
}
if(ID=="OBR"){
  title= "Cross block rate"
}
if(ID=="PPP"){
  title= "Predictor strength"
}

      

###############
(p4 <- ggplot(mp, aes(x=Input, y=Value, color=Mode, fill=Mode))+
    geom_line(size=1) +
    geom_ribbon(aes(ymin=L, ymax=H), alpha=0.3, colour = NA) +
    facet_wrap(. ~ Estimate, scales = "free_y", nrow=2) +
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

ggsave(paste0("ParamRecov_points_",ID,mode,".pdf"), p4, height=8,width=8)
}


