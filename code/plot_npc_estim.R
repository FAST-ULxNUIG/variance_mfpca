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
compute_npc <- function(NPC, true_percentage) {
    sum(true_percentage < NPC) + 1
}

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
results$N_lab <- factor(
    results$N,
    labels = c("$N = 25$", "$N = 50$", "$N = 100$")
)
results$M_lab <- factor(
    results$M,
    labels = c("$S = 25$", "$S = 50$", "$S = 100$")
)
results$lab <- interaction(results$N_lab, results$M_lab, sep = ' and ')

results_unique <- results |> 
    select(N, M, NPC) |> 
    unique()
results_unique <- results_unique |> 
    mutate(
        NPC_int = sapply(
            results_unique$NPC, compute_npc, true_percentage = true_percentage
        )
    )

gg <- ggplot(results) +
    geom_count(aes(x = as.factor(NPC), y = npc_estim)) +
    scale_size_area() +
    geom_point(
        aes(x = as.factor(NPC), y = NPC_int), data = results_unique, 
        color = 'red'
    ) +
    scale_y_continuous(breaks = seq(1, 10, by = 1)) +
    facet_wrap(vars(lab)) +
    labs(
        x = "Proportion of variance explained",
        y = "Number of multivariate components"
    ) +
    guides(size = guide_legend(title = 'Number of simulations')) +
    see::theme_modern() +
    theme(
        legend.position = "bottom", strip.text = element_text(size = 14)
    )

tikzDevice::tikz(
    filename = paste0(GRAPHS, '/npc_estim.tex'), 
    width = 10, height = 10, 
    standAlone = TRUE, sanitize = FALSE
)
plot(gg)
dev.off()
