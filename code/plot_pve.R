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
library(tikzDevice)

# -----------------------------------------------------------------------------
# Paths
PATH <- './results/pve'
GRAPHS <- './graphs/pve'

# -----------------------------------------------------------------------------
# Variables
K <- 50
true_eigenvalues <- eVal(K, 'exponential')
true_percentage <- cumsum(true_eigenvalues) / sum(true_eigenvalues)

# -----------------------------------------------------------------------------
# Functions
extract_eigenvalues <- function(filename) {
    result_list <- readRDS(paste0(PATH, '/', filename))
    
    # Get the different parameters
    N <- as.integer(
        str_extract(str_extract(filename, "N_[:digit:]+"), "[:digit:]+")
    )
    P <- as.integer(
        str_extract(str_extract(filename, "P_[:digit:]+"), "[:digit:]+")
    )
    M <- as.integer(
        str_extract(str_extract(filename, "M_[:digit:]+"), "[:digit:]+")
    )
    K <- as.integer(
        str_extract(str_extract(filename, "K_[:digit:]+"), "[:digit:]+")
    )
    NPC <- as.double(
        str_extract(str_extract(filename, "univ_0.[:digit:]+"), "0.[:digit:]+")
    )
    N_SIMU <- as.integer(
        str_extract(str_extract(filename, "[:digit:]+.rds"), "[:digit:]+")
    )
    
    eigenvalues <- result_list |>
        lapply(function(x) x$values)
    
    rm(result_list)
    gc()
    
    list(
        'N' = N, 'P' = P, 'M' = M, 'K' = K, 'NPC' = NPC, 'N_SIMU' = N_SIMU,
        eigenvalues = eigenvalues
    )
}

error_eigenvalues <- function(eigenvalues, eigenvalues_estim) {
    (eigenvalues - eigenvalues_estim)**2 / eigenvalues**2
}

plot_eigenvalues <- function(idx) {
    # Extract eigenvalues
    eigenvalues <- extract_eigenvalues(results_fls[idx])
    
    pct_explained <- eigenvalues$eigenvalues |> 
        sapply(function(x) sum(x) / sum(true_eigenvalues))
    NPC_estim <- eigenvalues$eigenvalues |> 
        sapply(function(x) sum((cumsum(x) / sum(x)) <= eigenvalues$NPC) + 1)

    NPC <- min(
        eigenvalues$K,
        min(eigenvalues$eigenvalues |> sapply(function(x) length(x)))
    )
    
    eigenvalues$errors <- eigenvalues$eigenvalues |> lapply(
        function(x) {
            vals <- true_eigenvalues[1:NPC]
            x <- x[1:NPC]
            error_eigenvalues(vals, x)
        }
    )
    
    errors_tbl <- tibble(
        number = rep(1:NPC, eigenvalues$N_SIMU),
        value = unlist(eigenvalues$errors)
    )
    
    title <- paste0(
        '$N = ', eigenvalues$N, ' - M = ', eigenvalues$M,
        ' - K_p = ', eigenvalues$NPC, '$'
    )
    gg <- ggplot(errors_tbl) +
        geom_boxplot(aes(x = number, y = value, group = number)) +
        labs(
            title = title,
            x = "Eigenvalues",
            y = "Errors"
        ) +
        see::theme_modern()
    
    name <- paste0(
        '/pve_eigenvalues_N_', eigenvalues$N, '_M_', eigenvalues$M,
        '_univ_', eigenvalues$NPC, '.tex'
    )
    tikzDevice::tikz(
        filename = paste0(GRAPHS, name), 
        width = 10, height = 6.18047, 
        standAlone = TRUE, sanitize = FALSE
    )
    plot(gg)
    dev.off()
    
    name <- paste0(
        '/pve_eigenvalues_N_', eigenvalues$N, '_M_', eigenvalues$M,
        '_univ_', eigenvalues$NPC, '.rds'
    )
    saveRDS(pct_explained, file = paste0('./results/pct_explained/', name))
    
    name <- paste0(
        '/npc_estim_N_', eigenvalues$N, '_M_', eigenvalues$M,
        '_univ_', eigenvalues$NPC, '.rds'
    )
    saveRDS(NPC_estim, file = paste0('./results/npc_estim/', name))
}

# -----------------------------------------------------------------------------
results_fls <- list.files(PATH, full.names = FALSE, pattern = '.rds')

for (idx in 1:length(results_fls)) {
    plot_eigenvalues(idx)
}
