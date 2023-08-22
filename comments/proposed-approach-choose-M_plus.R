# -----------------------------------------------------------------------------
# Script for "A note on the number of components retained for multivariate
# functional principal components analysis" for the estimation of a given
# number of multivariate eigencomponents
# -----------------------------------------------------------------------------

# Implements an idea:
# Step 1 = Use a high PVE to choose M_j in the univariate FPCA each dimensions.
# Step 2 = Choose M based on the eigenvalues of the mv-FPCA to achieve a given $\alpha%. Then retain only M mv-FCPS.

# -----------------------------------------------------------------------------
# Load packages
library(foreach)
library(funData)
library(MFPCA)
library(optparse)
library(tidyverse)


setwd("/Users/edwardgunning/Dropbox/Happ-Greven-Comment-Steven")
# -----------------------------------------------------------------------------
# Get parameters
option_list = list(
  make_option(
    c("-S", "--n_simu"), type = "integer", default = 10, 
    help = "Number of simulations", metavar = "integer"
  ),
  make_option(
    c("-N", "--n_curves"), type = "integer", default = 50, 
    help = "Number of curves", metavar = "integer"
  ),
  make_option(
    c("-P", "--n_features"), type = "integer", default = 5, 
    help = "Number of features", metavar = "integer"
  ),
  make_option(
    c("-M", "--n_points"), type = "integer", default = 101, 
    help = "Number of sampling points", metavar = "integer"
  ),
  make_option(
    c("-K", "--n_components"), type = "integer", default = 50, 
    help = "Number of components", metavar = "integer"
  ),
  make_option(
    c("-k", "--npc_univ"), type = "double", default = 0.9,
    help = "Percentage of variance explained by the univariate components",
    metavar = "integer"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# -----------------------------------------------------------------------------
# Register parallel back-end
n_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl = cl)

# -----------------------------------------------------------------------------
# Parameters
N_sim <- opt$n_simu  # Number of simulations

N <- opt$n_curves  # Number of curves
P <- opt$n_features  # Number of features
M <- opt$n_points  # Number of sampling points per curves
K <- opt$n_components  # Number of components to simulate the curves
pct_univ <- rep(opt$npc_univ, P)  # Number of univariate components

argvals <- seq(0, 1, length.out = M)
argvals_list <- rep(list(argvals), P)


# -----------------------------------------------------------------------------
# Define univariate expansion
uni_expansion_list <- vector("list", length = P)
for (idx in seq_len(P)) {
  uni_expansion <- list(type = "uFPCA", pve = pct_univ[idx], nbasis = 10)
  uni_expansion_list[idx] <- list(uni_expansion)
}

uni_expansion_lossless_list <- vector("list", length = P)
for (idx in seq_len(P)) {
  uni_expansion <- list(type = "uFPCA", pve = 0.9999, nbasis = 10)
  uni_expansion_lossless_list[idx] <- list(uni_expansion)
}



# Run simulations
results_list <- foreach(
  idx = seq_len(N_sim),
  .combine = 'c',
  .packages = c('funData', 'MFPCA')
) %dopar% {
  print(paste("Simulation", idx))
  
  # Simulate the data
  sim <- simMultiFunData(
    type = "split",
    argvals = argvals_list,
    M = K,
    eFunType = "Fourier",
    eValType = "exponential",
    N = N
  )
  
  # Univariate FPCA
  npc_univ <- rep(0, P)
  for (idx in 1:P) {
    pace <- MFPCA::PACE(
      funDataObject = sim$simData[[idx]], pve = pct_univ[idx]
    )
    npc_univ[idx] <- pace$npc
  }
  npc <- sum(npc_univ)
  
  # MFPCA
  MFPCA_est <- MFPCA(
    sim$simData,
    M = min(npc, K),
    uniExpansions = uni_expansion_list,
    fit = TRUE
  )
  
  # do MFPCA with lossless univariate expansions:
  MFPCA_lossless_init <- MFPCA(
    sim$simData,
    M = min(npc, K),
    uniExpansions = uni_expansion_lossless_list,
    fit = TRUE
  )
  
  # Choose the cutoff based on the MULTIVARIATE EIGENVALUES
  K_retain <- min(which(cumsum(MFPCA_lossless_init$values / sum(MFPCA_lossless_init$values)) > pct_univ[idx]))
  
  ## Retain only the first K_retain multivariate components.
  # (not sure how to do this perfectly, this is just an idea...)
  MFPCA_lossless_final <- MFPCA_lossless_init
  MFPCA_lossless_final$values <- MFPCA_lossless_final$values[seq_len(K_retain)]
  for (idx in 1:P) {
    MFPCA_lossless_final$functions[[idx]] <- MFPCA_lossless_final$functions[[idx]][seq_len(K_retain)]
  }
  MFPCA_lossless_final$scores <- MFPCA_lossless_final$scores[, seq_len(K_retain)]
  
  # return both ways...
  list(way_1 = MFPCA_est,
       way_2 = MFPCA_lossless_final)
}

# Stop the cluster
parallel::stopCluster(cl = cl)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Save the results
file_name <- paste0(
  'results_pve_N_', N, '_P_', P, '_M_', M, '_K_', K,
  '_univ_', opt$npc_univ, '_', N_sim, '.rds'
)
saveRDS(results_list, file = file_name)



