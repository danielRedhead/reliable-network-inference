################################################################# Load packages
library(igraph)
library(qgraph)
library(rethinking)
library(Rlab)
library(parallel)
library(cmdstanr)
library(reshape2)
library(plyr)
library(kinship2)
library(geosphere)
library(GGally)
library(network)
library(ggplot2)
library(rethinking)
library(colorspace)
library(parallel)
library(Cairo)
library(infotheo)
library(dirmult)
library(STRAND)

# Set your working directory, e.g.,
#setwd("./reliable-network-inference/Code/")

################################################################# Figure and analysis replication
source("OtherScripts/network_functions.R")
source("OtherScripts/figure_1_plots.R")
source("OtherScripts/colombia_example.R")

################################################################ Full robustness checks
# This code is also reproducible from here on.
# However, it will take an eternity to run on a normal computer. We recomend using at least a 100 core server, and then running
# each sub-block of code below in a separate R session.

################################################################# Basic Models
# Simulate and analyse fake data

# False positive mean rate
source("OtherScripts/network_functions.R")
source("BasicModels/plotting_functions.R")
source("BasicModels/false_positive_rate_test.R")
source("BasicModels/false_positive_rate_plots.R")

rm(list=ls(all=TRUE))

# False positive variance
source("OtherScripts/network_functions.R")
source("BasicModels/plotting_functions.R")
source("BasicModels/false_positive_var_test.R")
source("BasicModels/false_positive_var_plots.R")

rm(list=ls(all=TRUE))

# Recall of true ties mean rate
source("OtherScripts/network_functions.R")
source("BasicModels/plotting_functions.R")
source("BasicModels/recall_true_ties_rate_test.R")
source("BasicModels/recall_true_ties_rate_plots.R")

rm(list=ls(all=TRUE))

# Recall of true ties variance
source("OtherScripts/network_functions.R")
source("BasicModels/plotting_functions.R")
source("BasicModels/recall_true_ties_var_test.R")
source("BasicModels/recall_true_ties_var_plots.R")

rm(list=ls(all=TRUE))

# Mixed rate
source("OtherScripts/network_functions.R")
source("BasicModels/plotting_functions.R")
source("BasicModels/mixed_rate_test.R")
source("BasicModels/mixed_rate_plots.R")

rm(list=ls(all=TRUE))


################################################################# Question Order Models
# Simulate and analyse fake data

# Theta mean
source("OtherScripts/network_functions.R")
source("QuestionOrderModels/plotting_functions.R")
source("QuestionOrderModels/theta_rate_test.R")
source("QuestionOrderModels/theta_rate_plots.R")
rm(list=ls(all=TRUE))

# Theta variance
source("OtherScripts/network_functions.R")
source("QuestionOrderModels/plotting_functions.R")
source("QuestionOrderModels/theta_var_test.R")
source("QuestionOrderModels/theta_var_plots.R")
rm(list=ls(all=TRUE))


# False positive mean rate
source("OtherScripts/network_functions.R")
source("QuestionOrderModels/plotting_functions.R")
source("QuestionOrderModels/false_positive_rate_test.R")
source("QuestionOrderModels/false_positive_rate_plots.R")
rm(list=ls(all=TRUE))

# False positive variance
source("OtherScripts/network_functions.R")
source("QuestionOrderModels/plotting_functions.R")
source("QuestionOrderModels/false_positive_var_test.R")
source("QuestionOrderModels/false_positive_var_plots.R")
rm(list=ls(all=TRUE))

# Recall of true ties mean rate
source("OtherScripts/network_functions.R")
source("QuestionOrderModels/plotting_functions.R")
source("QuestionOrderModels/recall_true_ties_rate_test.R")
source("QuestionOrderModels/recall_true_ties_rate_plots.R")
rm(list=ls(all=TRUE))

# Recall of true ties variance
source("OtherScripts/network_functions.R")
source("QuestionOrderModels/plotting_functions.R")
source("QuestionOrderModels/recall_true_ties_var_test.R")
source("QuestionOrderModels/recall_true_ties_var_plots.R")
rm(list=ls(all=TRUE))

# Mixed rate
source("OtherScripts/network_functions.R")
source("QuestionOrderModels/plotting_functions.R")
source("QuestionOrderModels/mixed_rate_test.R")
source("QuestionOrderModels/mixed_rate_plots.R")
rm(list=ls(all=TRUE))

################################################################# Status Models
# Simulate and analyse fake data

# False positive mean rate
source("OtherScripts/network_functions.R")
source("StatusBiasModels/plotting_functions.R")
#source("StatusBiasModels/false_positive_rate_test.R")
source("StatusBiasModels/false_positive_rate_test_new.R")
source("StatusBiasModels/false_positive_rate_plots.R")
rm(list=ls(all=TRUE))

# False positive var
source("OtherScripts/network_functions.R")
source("StatusBiasModels/plotting_functions.R")
source("StatusBiasModels/false_positive_var_test.R")
source("StatusBiasModels/false_positive_var_plots.R")
rm(list=ls(all=TRUE))

# Recall of true ties mean rate
source("OtherScripts/network_functions.R")
source("StatusBiasModels/plotting_functions.R")
source("StatusBiasModels/recall_true_ties_rate_test.R")
source("StatusBiasModels/recall_true_ties_rate_plots.R")
rm(list=ls(all=TRUE))

# Recall of true ties var
source("OtherScripts/network_functions.R")
source("StatusBiasModels/plotting_functions.R")
source("StatusBiasModels/recall_true_ties_var_test.R")
source("StatusBiasModels/recall_true_ties_var_plots.R")
rm(list=ls(all=TRUE))

# mixed
source("OtherScripts/network_functions.R")
source("StatusBiasModels/plotting_functions.R")
source("StatusBiasModels/mixed_rate_test.R")
source("StatusBiasModels/mixed_rate_plots.R")
rm(list=ls(all=TRUE))

################################################################# Status Models with extra nets
# Simulate and analyse fake data

# False positive mean rate
source("OtherScripts/network_functions.R")
source("GroundTruthModels/plotting_functions.R")
source("GroundTruthModels/false_positive_rate_test.R")
source("GroundTruthModels/false_positive_rate_plots.R")
rm(list=ls(all=TRUE))

# False positive var
source("OtherScripts/network_functions.R")
source("GroundTruthModels/plotting_functions.R")
source("GroundTruthModels/false_positive_var_test.R")
source("GroundTruthModels/false_positive_var_plots.R")
rm(list=ls(all=TRUE))

# Recall of true ties mean rate
source("OtherScripts/network_functions.R")
source("GroundTruthModels/plotting_functions.R")
source("GroundTruthModels/recall_true_ties_rate_test.R")
source("GroundTruthModels/recall_true_ties_rate_plots.R")
rm(list=ls(all=TRUE))

# Recall of true ties var
source("OtherScripts/network_functions.R")
source("GroundTruthModels/plotting_functions.R")
source("GroundTruthModels/recall_true_ties_var_test.R")
source("GroundTruthModels/recall_true_ties_var_plots.R")
rm(list=ls(all=TRUE))

# mixed
source("OtherScripts/network_functions.R")
source("GroundTruthModels/plotting_functions.R")
source("GroundTruthModels/mixed_rate_test.R")
source("GroundTruthModels/mixed_rate_plots.R")
rm(list=ls(all=TRUE))


################################################################# SRM + SBM Models
# Simulate and analyse fake data

# DR rho
source("OtherScripts/network_functions.R")
source("SRMandSBM/plotting_functions.R")
source("SRMandSBM/DR_rho_test.R")
source("SRMandSBM/DR_rho_plot.R")
rm(list=ls(all=TRUE))

# DR sigma
source("OtherScripts/network_functions.R")
source("SRMandSBM/plotting_functions.R")
source("SRMandSBM/DR_sigma_test.R")
source("SRMandSBM/DR_sigma_plot.R")
rm(list=ls(all=TRUE))

# in-block
source("OtherScripts/network_functions.R")
source("SRMandSBM/plotting_functions.R")
source("SRMandSBM/in_block_rate_test.R")
source("SRMandSBM/in_block_rate_plot.R")
rm(list=ls(all=TRUE))

# out-block
source("OtherScripts/network_functions.R")
source("SRMandSBM/plotting_functions.R")
source("SRMandSBM/out_block_rate_test.R")
source("SRMandSBM/out_block_rate_plot.R")
rm(list=ls(all=TRUE))

# SR rho
source("OtherScripts/network_functions.R")
source("SRMandSBM/plotting_functions.R")
source("SRMandSBM/SR_rho_test.R")
source("SRMandSBM/SR_rho_plot.R")
rm(list=ls(all=TRUE))

# SR sigma
source("OtherScripts/network_functions.R")
source("SRMandSBM/plotting_functions.R")
source("SRMandSBM/SR_sigma_test.R")
source("SRMandSBM/SR_sigma_plot.R")
rm(list=ls(all=TRUE))

# SR predictor
source("OtherScripts/network_functions.R")
source("SRMandSBM/plotting_functions.R")
source("SRMandSBM/SR_predictor_test.R")
source("SRMandSBM/SR_predictor_plot.R")
rm(list=ls(all=TRUE))



################################################################# SRM 
# Simulate and analyse fake data

# DR rho
source("OtherScripts/network_functions.R")
source("SRM/plotting_functions.R")
source("SRM/DR_rho_test.R")
source("SRM/DR_rho_plot.R")
rm(list=ls(all=TRUE))

# DR sigma
source("OtherScripts/network_functions.R")
source("SRM/plotting_functions.R")
source("SRM/DR_sigma_test.R")
source("SRM/DR_sigma_plot.R")
rm(list=ls(all=TRUE))

# SR rho
source("OtherScripts/network_functions.R")
source("SRM/plotting_functions.R")
source("SRM/SR_rho_test.R")
source("SRM/SR_rho_plot.R")
rm(list=ls(all=TRUE))

# SR sigma
source("OtherScripts/network_functions.R")
source("SRM/plotting_functions.R")
source("SRM/SR_sigma_test.R")
source("SRM/SR_sigma_plot.R")
rm(list=ls(all=TRUE))

# tie rate
source("OtherScripts/network_functions.R")
source("SRM/plotting_functions.R")
source("SRM/tie_rate_test.R")
source("SRM/tie_rate_plot.R")
rm(list=ls(all=TRUE))

# SR predictor
source("OtherScripts/network_functions.R")
source("SRM/plotting_functions.R")
source("SRM/SR_predictor_test.R")
source("SRM/SR_predictor_plot.R")
rm(list=ls(all=TRUE))


################################################################# SBM 
# Simulate and analyse fake data

# in block
source("OtherScripts/network_functions.R")
source("SBM/plotting_functions.R")
source("SBM/in_block_rate_test.R")
source("SBM/in_block_rate_plot.R")
rm(list=ls(all=TRUE))

# out block
source("OtherScripts/network_functions.R")
source("SBM/plotting_functions.R")
source("SBM/out_block_rate_test.R")
source("SBM/out_block_rate_plot.R")
rm(list=ls(all=TRUE))

# n blocks
source("OtherScripts/network_functions.R")
source("SBM/plotting_functions.R")
source("SBM/n_blocks_test.R")
source("SBM/n_blocks_plot.R")
rm(list=ls(all=TRUE))

# n_sample blocks
source("OtherScripts/network_functions.R")
source("SBM/plotting_functions.R")
source("SBM/n_sample_block_test.R")
source("SBM/n_sample_block_plot.R")
rm(list=ls(all=TRUE))

# block size
source("OtherScripts/network_functions.R")
source("SBM/plotting_functions.R")
source("SBM/block_size_test.R")
source("SBM/block_size_plot.R")
rm(list=ls(all=TRUE))

# predictor
source("OtherScripts/network_functions.R")
source("SBM/plotting_functions.R")
source("SBM/predictor_test.R")
source("SBM/predictor_plot.R")
rm(list=ls(all=TRUE))



############################################## Colombia example
