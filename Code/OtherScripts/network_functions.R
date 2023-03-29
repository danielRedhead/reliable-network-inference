######################################################################## For plots

get_network_level_properties <- function(adjmatrix){
  res <- rep(NA, length=5)
  graph <- graph_from_adjacency_matrix(adjmatrix, mode = "directed", weighted = NULL, diag = FALSE)

  res[1] <- edge_density(graph)
  res[2] <- reciprocity(graph, ignore.loops = TRUE)
  res[3] <- transitivity(graph) 
  res[4] <- centr_betw(graph, directed = FALSE)$centralization
  res[5] <- centr_eigen(graph, directed = FALSE)$centralization

return(res)
}


get_node_level_properties <- function(adjmatrix1,adjmatrix2){
  res <- rep(NA, length=5)
  graph1 <- graph_from_adjacency_matrix(adjmatrix1, mode = "directed", weighted = NULL, diag = FALSE)
  graph2 <- graph_from_adjacency_matrix(adjmatrix2, mode = "directed", weighted = NULL, diag = FALSE)

  res[1] <- cor(degree(graph1, mode = "in", loops = FALSE, normalized = FALSE), degree(graph2, mode = "in", loops = FALSE, normalized = FALSE))
  res[2] <- cor(degree(graph1, mode = "out", loops = FALSE, normalized = FALSE), degree(graph2, mode = "out", loops = FALSE, normalized = FALSE))
  res[3] <- cor(centr_betw(graph1, directed = FALSE)$res, centr_betw(graph2, directed = FALSE)$res ) 
  res[4] <- cor(centr_eigen(graph1, directed = FALSE)$vector, centr_eigen(graph2, directed = FALSE)$vector)
  res[5] <- cor(page_rank(graph1, directed = FALSE)$vector, page_rank(graph2, directed = FALSE)$vector)

  return(res)
}

get_node_level_properties_MI <- function(adjmatrix1,adjmatrix2){
  res <- rep(NA, length=5)
  graph1 <- graph_from_adjacency_matrix(adjmatrix1, mode = "directed", weighted = NULL, diag = FALSE)
  graph2 <- graph_from_adjacency_matrix(adjmatrix2, mode = "directed", weighted = NULL, diag = FALSE)

  res[1] <- n_mut_inf(degree(graph1, mode = "in", loops = FALSE, normalized = FALSE), degree(graph2, mode = "in", loops = FALSE, normalized = FALSE))
  res[2] <- n_mut_inf(degree(graph1, mode = "out", loops = FALSE, normalized = FALSE), degree(graph2, mode = "out", loops = FALSE, normalized = FALSE))
  res[3] <- n_mut_inf(centr_betw(graph1, directed = FALSE)$res, centr_betw(graph2, directed = FALSE)$res ) 
  res[4] <- n_mut_inf(centr_eigen(graph1, directed = FALSE)$vector, centr_eigen(graph2, directed = FALSE)$vector)
  res[5] <- n_mut_inf(page_rank(graph1, directed = FALSE)$vector, page_rank(graph2, directed = FALSE)$vector)

  return(res)
}

strand_caterpillar_plot = function(results, submodels=NULL, normalized=FALSE, only_slopes=TRUE, only_technicals=FALSE){
  dat = vector("list",length(results$summary_list))

  for(k in 1:length(results$summary_list)){
   dat[[k]] = data.frame(results$summary_list[[k]])
   dat[[k]]$SubModel = names(results$summary_list)[k]
   colnames(dat[[k]]) = c("Variable", "Median", "LI", "HI", "Mean","SD", "SubModel")
   for(j in 2:6)
   dat[[k]][,j] = as.numeric(dat[[k]][,j])
  }


df = do.call(rbind, dat)

colnames(df) = c("Variable", "Median", "LI", "HI", "Mean","SD", "SubModel")


df$Submodel = factor(df$SubModel)
df$Submodel = factor(df$SubModel, levels=c("False positive rate", "Recall of true ties","Theta: question-order effects",
                                           "Focal efffects: Out-degree","Target effects: In-degree","Dyadic effects", "Other estimates" ))

if(only_slopes==TRUE){
exclude=c("false positive rate intercept, layer 1",               
"false positive rate intercept, layer 2",               
"false positive rate sd, layer 1",                     
"false positive rate sd, layer 2",                            
"recall rate of true ties intercept, layer 1",          
"recall rate of true ties intercept, layer 2",          
"recall rate of true ties sd, layer 1",                 
"recall rate of true ties sd, layer 2",                 
"theta intercept, layer 1 to 2",                        
"theta sd, layer 1 to 2",                               
"focal effects sd",                                            
"target effects sd",                                            
"dyadic effects sd",                                                         
"focal-target effects rho (generalized recipocity)",    
"dyadic effects rho (dyadic recipocity)")   

df = df[which(!df$Variable %in% exclude),]
}

if(only_technicals==TRUE){
include=c("false positive rate intercept, layer 1",               
"false positive rate intercept, layer 2",               
"false positive rate sd, layer 1",                     
"false positive rate sd, layer 2",                            
"recall rate of true ties intercept, layer 1",          
"recall rate of true ties intercept, layer 2",          
"recall rate of true ties sd, layer 1",                 
"recall rate of true ties sd, layer 2",                 
"theta intercept, layer 1 to 2",                        
"theta sd, layer 1 to 2",                               
"focal effects sd",                                            
"target effects sd",                                            
"dyadic effects sd",                                                         
"focal-target effects rho (generalized recipocity)",    
"dyadic effects rho (dyadic recipocity)")   

unit=c("false positive rate intercept, layer 1",               
"false positive rate intercept, layer 2",                                  
"recall rate of true ties intercept, layer 1",          
"recall rate of true ties intercept, layer 2",                          
"theta intercept, layer 1 to 2",                                                                                                             
"focal-target effects rho (generalized recipocity)",    
"dyadic effects rho (dyadic recipocity)")

df = df[which(df$Variable %in% include),]

df$Scaling = ifelse(df$Variable %in% unit, "Rates", "Dispersion")
}


if(!is.null(submodels))
df = df[which(df$SubModel %in% submodels),]

df$Diff = df$HI-df$LI   

if(normalized==TRUE) {
  df$Median = df$Median/df$Diff
  df$LI = df$LI/df$Diff
  df$HI =  df$HI/df$Diff
}
 

p <- ggplot(df,aes(x=Variable,y=Median,ymin=LI,ymax=HI))+ 
     geom_linerange(size=1)+
     geom_point(size=2)+
     facet_grid( SubModel~., scales = "free", space='free')+
       #facet_wrap(vars(SubModel), ncol=1,scales = "free")+
       geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
     labs(y="Regression parameters", x="") + theme(strip.text.x = element_text(size=12,face="bold"), 
     strip.text.y = element_text(size=12,face="bold"),axis.text=element_text(size=12),axis.title.y=element_text(size=14,
     face="bold"), axis.title.x=element_blank())+theme(strip.text.y = element_text(angle = 360)) + coord_flip() + theme(panel.spacing = unit(1, "lines")) 

p2 <- ggplot(df,aes(x=Variable,y=Median,ymin=LI,ymax=HI))+ 
     geom_linerange(size=1)+
     geom_point(size=2)+
     facet_grid( SubModel~Scaling, scales = "free")+
       #facet_wrap(vars(SubModel), ncol=1,scales = "free")+
       geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
     labs(y="Regression parameters", x="") + theme(strip.text.x = element_text(size=8,face="bold"), 
     strip.text.y = element_text(size=8,face="bold"),axis.text=element_text(size=8),axis.title.y=element_text(size=14,
     face="bold"), axis.title.x=element_blank())+theme(strip.text.y = element_text(angle = 360)) + coord_flip() + theme(panel.spacing = unit(1, "lines")) 


if(only_technicals==TRUE){
 return(p2)} else{
    return(p)
 }
 }