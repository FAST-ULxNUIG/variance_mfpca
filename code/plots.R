#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Script for "A note on the number of components retained for multivariate
# functional principal components analysis" for the estimation of a given
# number of multivariate eigencomponents
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Load packages
library(funData)
library(MFPCA)
library(tidyverse)

# -----------------------------------------------------------------------------
# Exponential
# -----------------------------------------------------------------------------
true_vals <- eVal(K, 'exponential')
true_vals_df <- tibble(
    variable = 1:length(true_vals),
    value = true_vals
)


# -----------------------------------------------------------------------------
# Extract results for eigenvalues
# -----------------------------------------------------------------------------
eigenvalues <- results_list |>
    lapply(function(x) x$values)

eigenvalues_compare <- simu_list |>
    lapply(function(x) abs(true_vals[1:npc] - x$values) / true_vals[1:npc])

eigenvalues_tbl <- tibble(
    n_simu = rep(1:N_sim, each = npc),
    N = rep(1:npc, times = N_sim),
    value = unlist(eigenvalues),
    values_compare = unlist(eigenvalues_compare)
)

# -----------------------------------------------------------------------------
# Extract results for eigenfunctions
# -----------------------------------------------------------------------------
compare_mfd <- function(true_curve, estim_curve) {
    results <- rep(0, nObs(estim_curve))
    for (idx in 1:nObs(estim_curve)) {
        results[idx] = min(
            norm(true_curve[idx] - estim_curve[idx])**2,
            norm(true_curve[idx] + estim_curve[idx])**2
        )
    }
    return(results)
}

eigenfunctions <- simu_list |>
    lapply(function(x) x$functions)

eigenfunctions_compare <- simu_list |>
    lapply(function(x) compare_mfd(sim$trueFuns[1:npc], x$functions))

eigenfunctions_tbl <- tibble(
    n_simu = rep(1:N_sim, each = npc),
    N = rep(1:npc, times = N_sim),
    values_compare = unlist(eigenfunctions_compare)
)


# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------

ggplot(eigenvalues_tbl) +
    geom_boxplot(aes(x = N, y = values_compare, group = N)) +
    xlab('Eigenvalues') +
    ylab('Error') +
    see::theme_modern()

ggplot(eigenfunctions_tbl) +
    geom_boxplot(aes(x = N, y = values_compare, group = N)) +
    xlab('Eigenvalues') +
    ylab('Error') +
    see::theme_modern()
