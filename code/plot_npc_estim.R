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
PATH <- './results/npc_estim'
GRAPHS <- './graphs/npc_estim'

# -----------------------------------------------------------------------------
# Variables
K <- 50
true_eigenvalues <- eVal(K, 'exponential')
true_percentage <- cumsum(true_eigenvalues) / sum(true_eigenvalues)


# -----------------------------------------------------------------------------
# Functions
extract_npc <- function(filename) {
    result_list <- readRDS(paste0(PATH, '/', filename))
    
    # Get the different parameters
    N <- as.integer(
        str_extract(str_extract(filename, "N_[:digit:]+"), "[:digit:]+")
    )
    M <- as.integer(
        str_extract(str_extract(filename, "M_[:digit:]+"), "[:digit:]+")
    )
    NPC <- as.double(
        str_extract(str_extract(filename, "univ_0.[:digit:]+"), "0.[:digit:]+")
    )
    
    tibble(
        'N' = N, 'M' = M, 'NPC' = NPC,
        npc_estim = result_list
    )
}

# -----------------------------------------------------------------------------
results_fls <- list.files(PATH, full.names = FALSE, pattern = '.rds')

results <- results_fls |> 
    lapply(function(x) extract_npc(x)) |> 
    bind_rows()


results_unique <- results |> select(N, M, NPC) |> unique()

ggplot(results) +
    geom_count(aes(y = npc_estim, x = as.factor(NPC))) +
    #geom_line(aes(y = NPC, x = NPC), color = 'red') +
    facet_grid(rows = vars(N), cols = vars(M)) + 
    labs(
        x = "Percentage of variance explained",
        y = "NPC"
    ) +
    see::theme_modern()


