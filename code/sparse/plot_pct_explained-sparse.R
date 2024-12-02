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
PATH <- './sparse/results/pct_explained'
GRAPHS <- './graphs/sparse/pct_explained'

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
        pct_explained = result_list
    )
}

# -----------------------------------------------------------------------------
results_fls <- list.files(PATH, full.names = FALSE, pattern = '.rds')

results <- results_fls |> 
    lapply(function(x) extract_npc(x)) |> 
    bind_rows() |> 
    arrange(N, M)
results$N_lab <- factor(
    results$N,
    labels = c("$N = 25$", "$N = 50$", "$N = 100$")
)
results$M_lab <- factor(
    results$M,
    labels = c("$M = 25$", "$M = 50$", "$M = 100$")
)
results$lab <- interaction(results$N_lab, results$M_lab, sep = ' and ')

results_unique <- results |> 
    select(N, M, NPC) |> 
    unique() |> 
    mutate(
        NPC_int = as.numeric(as.factor(NPC))
    )

gg <- ggplot(results) +
    geom_boxplot(aes(y = pct_explained, x = as.factor(NPC))) +
    geom_segment(
        aes(x = NPC_int - 0.5, xend = NPC_int + 0.5,
            y = NPC, yend = NPC), data = results_unique,
        color = 'red'
    ) +
    facet_wrap(vars(lab)) + 
    labs(
        x = "True percentage explained",
        y = "Percentage explained estimed"
    ) +
    see::theme_modern() +
    theme(
        legend.position = "bottom", strip.text = element_text(size = 14)
    )

tikzDevice::tikz(
    filename = paste0(GRAPHS, '/pct_estim.tex'), 
    width = 10, height = 10, 
    standAlone = TRUE, sanitize = FALSE
)
plot(gg)
dev.off()
