# -----------------------------------------------------------------------------
# Script for "A note on the number of components retained for multivariate
# functional principal components analysis" for the estimation of a given
# number of multivariate eigencomponents
# -----------------------------------------------------------------------------

# Load packages
library(funData)
library(MFPCA)
library(tidyverse)

# -----------------------------------------------------------------------------
# Parameters
N_sim <- 50  # Number of simulations

N <- 100  # Number of curves
P <- 5  # Number of features
M <- 100  # Number of sampling points per curves
K <- 50  # Number of components to simulate the curves
npc <- 5  # Number of (multivariate) components to estimate
npc_univ <- seq(npc, K, by = 5)  # Number of univariate components to estimate
# Maybe it is the other way round.
# We fix the number of univariate components to be npc and we estimate the 
# number of multivarirate components to be seq(npc, sum_p npc).

argvals <- seq(0, 1, length.out = M)
argvals_list <- rep(list(argvals), P)

# -----------------------------------------------------------------------------
# Exponential
true_vals <- eVal(K, 'exponential')
true_vals_df <- tibble(
    variable = 1:length(true_vals),
    value = true_vals
)

simu_list <- vector("list", length = N_sim)
for (idx in seq_len(N_sim)) {
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

    evals_list <- vector("list", length = length(npc_univ))
    for (i in 1:length(npc_univ)) {
        print(paste("Univariate NPC:", npc_univ[i]))

        # Univariate expansion
        uni_expansion <- list(type = "uFPCA", npc = npc_univ[i], nbasis = 10)
        uni_expansion_list <- rep(list(uni_expansion), P)
        
        # MFPCA
        MFPCA_est <- MFPCA(
            sim$simData,
            M = npc, 
            uniExpansions = uni_expansion_list,
            fit = TRUE
        )
        
        evals_list[[i]] <- MFPCA_est
    }
    
    simu_list[[idx]] <- evals_list

}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Extract results from one simulation
# -----------------------------------------------------------------------------
eigenvalues <- evals_list |> 
    lapply(function(x) x$values)

eigenvalues_tbl <- tibble(
    K = rep(npc_univ, each = npc),
    N = rep(1:npc, length(npc_univ)),
    value = unlist(eigenvalues)
)

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------

idx <- 5

par(mfrow = c(1, 2))
plot(
    sim$trueFuns[[1]]@argvals[[1]],
    sim$trueFuns[[1]]@X[idx,],
    col = 'red', type = 'l'
)
lines(
    MFPCA_est$functions[[1]]@argvals[[1]],
    MFPCA_est$functions[[1]]@X[idx,],
    col = 'blue'
)

plot(
    sim$trueFuns[[2]]@argvals[[1]],
    sim$trueFuns[[2]]@X[idx,],
    col = 'red', type = 'l'
)
lines(
    MFPCA_est$functions[[2]]@argvals[[1]],
    MFPCA_est$functions[[2]]@X[idx,],
    col = 'blue'
)
