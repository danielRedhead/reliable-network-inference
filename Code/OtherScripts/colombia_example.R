# Load data
load("ColombianDataWithImputations.RData")
set.seed(1)

# Prepare for Stan
data(FoodSharing_Data)
g = FoodSharing_Data
g$Male = ifelse(g$Sex=="Male",1,0)

Outgoing_Transfers = g$TransferOut 
Incoming_Transfers = g$TransferIn 

net = list(
 TransferOut=g$TransferOut,
 TransferIn=g$TransferIn
 )

ind = data.frame(
 Age=center(g$Age), 
 GoodsValues=center(g$GoodsValues), 
 Male=g$Male,
 CantWork=g$CantWork, 
 Grip=center(g$GripStrength), 
 NoFood=g$NoFood, 
 Depressed=g$Depressed
 )

dyad = list(
 Relatedness=g$Relatedness, 
 Friends=g$Friends
 )

groups = data.frame(Ethnicity = 
 as.factor(g$Ethnicity)
 )

model_dat = make_strand_data(
 self_report=net, 
 block_covariates=groups, 
 individual_covariates=ind, 
 dyadic_covariates=dyad
 )

fit1 = fit_latent_network_model(
 data=model_dat,
 block_regression = ~ Ethnicity,
 focal_regression = ~ Age + GoodsValues + Male + NoFood + CantWork + Grip + Depressed,
 target_regression = ~ Age + GoodsValues + Male + NoFood + CantWork + Grip + Depressed,
 dyad_regression = ~ Relatedness + Friends,
 fpr_regression = ~ Age + GoodsValues + Depressed,
 rtt_regression = ~ Age + GoodsValues + Depressed,
 theta_regression = ~ 1,
 mode="mcmc",
 return_latent_network = TRUE
 )

res1 = summarize_strand_results(fit1)

###################################### Plotting
diag(Outgoing_Transfers) = 0
diag(Incoming_Transfers) = 0

group_ids= as.factor(g$Ethnicity)

outgoing = graph_from_adjacency_matrix(Outgoing_Transfers, mode = c("directed"))
incoming = graph_from_adjacency_matrix(Incoming_Transfers, mode = c("directed"))

reported_either <- ifelse(Outgoing_Transfers==1 | t(Incoming_Transfers)==1, 1,0)
reported_both <- ifelse(Outgoing_Transfers==1 & t(Incoming_Transfers)==1, 1,0)

union <- graph_from_adjacency_matrix(reported_either, mode = c("directed"))
intersection <- graph_from_adjacency_matrix(reported_both, mode = c("directed"))

latent_res = graph_from_adjacency_matrix(round(apply(res1$samples$latent_network_sample,2:3,median )), mode = c("directed"))

coords <- layout_with_kk(union)

V(latent_res)$color <-V(intersection)$color <- V(union)$color <- V(outgoing)$color <- V(incoming)$color <- c("gray13","turquoise4",  "goldenrod3")[group_ids]
V(latent_res)$frame.color <-V(intersection)$frame.color <- V(union)$frame.color <- V(outgoing)$frame.color <- V(incoming)$frame.color <- c("gray13", "turquoise4", "goldenrod3")[group_ids]


# plot true out network
Cairo(file="colombia_outgoing.png", 
      type="png",
      units="in", 
      width=6, 
      height=6, 
      pointsize=36, 
      dpi=150)

par(mar=c(1,1,1,1)+0)
plot(outgoing, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = 3, #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)

dev.off()

# plot true out network
Cairo(file="colombia_incoming.png", 
      type="png",
      units="in", 
      width=6, 
      height=6, 
      pointsize=36, 
      dpi=150)

par(mar=c(1,1,1,1)+0)
plot(incoming, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = 3, #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)

dev.off()

# plot true out network
Cairo(file="colombia_intersection.png", 
      type="png",
      units="in", 
      width=6, 
      height=6, 
      pointsize=36, 
      dpi=150)

par(mar=c(1,1,1,1)+0)
plot(intersection, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = 3, #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)
dev.off()


# plot true out network
Cairo(file="colombia_union.png", 
      type="png",
      units="in", 
      width=6, 
      height=6, 
      pointsize=36, 
      dpi=150)

par(mar=c(1,1,1,1)+0)
plot(union, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = 3, #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)
dev.off()

# plot true out network
Cairo(file="colombia_latent.png", 
      type="png",
      units="in", 
      width=6, 
      height=6, 
      pointsize=36, 
      dpi=150)

par(mar=c(1,1,1,1)+0)
plot(latent_res, edge.arrow.size =0.1, edge.curved = 0.3,
                  vertex.label=NA, vertex.size = 3, #vertex.color="#252525", vertex.frame.color="#252525",
                  layout = coords, seed = 3)
dev.off()




  vis1 = strand_caterpillar_plot(res1, submodel=c("Focal efffects: Out-degree","Target effects: In-degree","Dyadic effects"), normalized=TRUE)
  ggsave("Colombia_slopes_latent.pdf", vis1, width=8,height=8)

  vis2 = strand_caterpillar_plot(res1, submodel=c("False positive rate", "Recall of true ties","Theta: question-order effects"), normalized=TRUE)
  ggsave("Colombia_slopes_measurement.pdf", vis2, width=8,height=8)

  vis3 = strand_caterpillar_plot(res1, submodel=c("False positive rate", "Recall of true ties","Theta: question-order effects"), only_slopes=FALSE, 
                                  normalized=FALSE, only_technicals=TRUE)
  ggsave("Colombia_intercepts_measurement.pdf", vis3, width=8,height=4)

            


strand_caterpillar_plot2=function (
      results = res1
      submodels = NULL
      normalized = FALSE
      only_slopes = FALSE
    only_technicals = TRUE
    site = "BOB"
    export_as_table = FALSE

    dat = vector("list", length(results$summary_list))
    for (k in 1:length(results$summary_list)) {
        dat[[k]] = data.frame(results$summary_list[[k]])
        dat[[k]]$SubModel = names(results$summary_list)[k]
        colnames(dat[[k]]) = c("Variable", "Median", "LI", "HI", 
            "Mean", "SD", "SubModel")
        for (j in 2:6) dat[[k]][, j] = as.numeric(dat[[k]][, 
            j])
    }
    df = do.call(rbind, dat)
    colnames(df) = c("Variable", "Median", "LI", "HI", "Mean", 
        "SD", "SubModel")
    df$Submodel = factor(df$SubModel)
    df$Submodel = factor(df$SubModel, levels = c("False positive rate", 
        "Recall of true ties", "Theta: question-order effects", 
        "Focal efffects: Out-degree", "Target effects: In-degree", 
        "Dyadic effects", "Other estimates"))
    if (only_slopes == TRUE) {
        exclude = c("false positive rate intercept, layer 1", 
            "false positive rate intercept, layer 2", "false positive rate sd, layer 1", 
            "false positive rate sd, layer 2", "recall rate of true ties intercept, layer 1", 
            "recall rate of true ties intercept, layer 2", "recall rate of true ties sd, layer 1", 
            "recall rate of true ties sd, layer 2", "theta intercept, layer 1 to 2", 
            "theta sd, layer 1 to 2", "focal effects sd", "target effects sd", 
            "dyadic effects sd", "focal-target effects rho (generalized recipocity)", 
            "dyadic effects rho (dyadic recipocity)")
        df = df[which(!df$Variable %in% exclude), ]
    }
    if (only_technicals == TRUE) {
        include = c("false positive rate intercept, layer 1", 
            "false positive rate intercept, layer 2", "false positive rate sd, layer 1", 
            "false positive rate sd, layer 2", "recall rate of true ties intercept, layer 1", 
            "recall rate of true ties intercept, layer 2", "recall rate of true ties sd, layer 1", 
            "recall rate of true ties sd, layer 2", "theta intercept, layer 1 to 2", 
            "theta sd, layer 1 to 2", "focal effects sd", "target effects sd", 
            "dyadic effects sd", "focal-target effects rho (generalized recipocity)", 
            "dyadic effects rho (dyadic recipocity)")
        unit = c("false positive rate intercept, layer 1", "false positive rate intercept, layer 2", 
            "recall rate of true ties intercept, layer 1", "recall rate of true ties intercept, layer 2", 
            "theta intercept, layer 1 to 2", "focal-target effects rho (generalized recipocity)", 
            "dyadic effects rho (dyadic recipocity)")
        df = df[which(df$Variable %in% include), ]
        df$Scaling = ifelse(df$Variable %in% unit, "Rates", "Dispersion")
        df$SubModel2 = ifelse(df$SubModel == "Other estimates", 
            "Correlation", "Dispersion")
    }
    if (!is.null(submodels)) 
        df = df[which(df$SubModel %in% submodels), ]
    df$Diff = df$HI - df$LI
    if (normalized == TRUE) {
        df$Median = df$Median/df$Diff
        df$LI = df$LI/df$Diff
        df$HI = df$HI/df$Diff
    }
    df$Site = site
    df$Variable[which(df$Variable == "focal-target effects rho (generalized recipocity)")] = "focal-target effects rho"
    df$Variable[which(df$Variable == "dyadic effects rho (dyadic recipocity)")] = "dyadic effects rho"
    df$Submodel = factor(df$SubModel)
    df$Submodel = factor(df$SubModel, levels = c("False positive rate", 
        "Recall of true ties", "Theta: question-order effects", 
        "Focal efffects: Out-degree", "Target effects: In-degree", 
        "Dyadic effects", "Other estimates"))

    df=df[which(df$SubModel %in% c("False positive rate", 
        "Recall of true ties", "Theta: question-order effects")),]
    p <- ggplot(df, aes(x = Variable, y = Median, ymin = LI, 
        ymax = HI)) + geom_linerange(size = 1) + geom_point(size = 2) + 
        facet_grid(Submodel ~ ., scales = "free", space = "free") + 
        geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") + 
        labs(y = "Regression parameters", x = "") + theme(strip.text.x = element_text(size = 12, 
        face = "bold"), strip.text.y = element_text(size = 12, 
        face = "bold"), axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 14, face = "bold"), 
        axis.title.x = element_blank()) + theme(strip.text.y = element_text(angle = 360)) + 
        coord_flip() + theme(panel.spacing = unit(1, "lines"))
    p2 <- ggplot(df, aes(x = Variable, y = Median, ymin = LI, 
        ymax = HI)) + geom_linerange(size = 1) + geom_point(size = 2) + 
        facet_grid(SubModel ~ Scaling, scales = "free") + geom_hline(aes(yintercept = 0), 
        color = "black", linetype = "dashed") + labs(y = "Regression parameters", 
        x = "") + theme(strip.text.x = element_text(size = 12, 
        face = "bold"), strip.text.y = element_text(size = 12, 
        face = "bold"), axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 14, face = "bold"), 
        axis.title.x = element_blank()) + theme(strip.text.y = element_text(angle = 360)) + 
        coord_flip() + theme(panel.spacing = unit(1, "lines"))
    if (export_as_table == FALSE) {
        if (only_technicals == TRUE) {
            return(p2)
        }
        else {
            return(p)
        }
    }
    if (export_as_table == TRUE) {
        df
    }
}

  vis3 = strand_caterpillar_plot2(res1, submodel=c("False positive rate", "Recall of true ties","Theta: question-order effects"), only_slopes=FALSE, 
                                  normalized=FALSE, only_technicals=TRUE)