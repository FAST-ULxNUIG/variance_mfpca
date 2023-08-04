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
PATH <- './results/ncomp'
GRAPHS <- './graphs/ncomp'


# -----------------------------------------------------------------------------
# Variables
K <- 50
true_eigenvalues <- eVal(K, 'exponential')


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
    NPC <- as.integer(
        str_extract(str_extract(filename, "univ_[:digit:]+"), "[:digit:]+")
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

create_dataframes <- function(idx) {
    # Extract eigenvalues
    eigenvalues <- extract_eigenvalues(results_fls[idx])
    
    NPC <- min(eigenvalues$K, eigenvalues$P * eigenvalues$NPC)
    
    eigenvalues$errors <- eigenvalues$eigenvalues |> lapply(
        function(x) {
            vals <- true_eigenvalues[1:NPC]
            error_eigenvalues(vals, x)
        }
    )
    
    return(eigenvalues)
}

create_errors <- function(idx) {
    
    eigenvalues <- create_dataframes(idx)
    NPC <- min(eigenvalues$K, eigenvalues$P * eigenvalues$NPC)
    errors_tbl <- tibble(
        N = rep(eigenvalues$N, eigenvalues$N_SIMU * NPC),
        M = rep(eigenvalues$M, eigenvalues$N_SIMU * NPC),
        P = rep(eigenvalues$P, eigenvalues$N_SIMU * NPC),
        number = rep(1:NPC, eigenvalues$N_SIMU),
        NPC = rep(eigenvalues$NPC, eigenvalues$N_SIMU * NPC),
        value = unlist(eigenvalues$errors)
    )
    
    return(errors_tbl)
}

plot_eigenvalues <- function(idx) {

    errors_tbl <- create_errors(idx)

    title = paste0(
        '$N = ', errors_tbl$N, ' - M = ', errors_tbl$M,
        ' - K_p = ', errors_tbl$NPC, '$'
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
        '/ncomp_eigenvalues_N_', errors_tbl$N, '_M_', errors_tbl$M,
        '_univ_', errors_tbl$NPC, '.tex'
    )
    tikzDevice::tikz(
        filename = paste0(GRAPHS, name), 
        width = 10, height = 6.18047, 
        standAlone = TRUE, sanitize = FALSE
    )
    plot(gg)
    dev.off()
}


# -----------------------------------------------------------------------------
results_fls <- list.files(PATH, full.names = FALSE, pattern = '.rds')

for (idx in 1:length(results_fls)) {
    plot_eigenvalues(idx)
}


# -----------------------------------------------------------------------------
errors_concat <- 1:length(results_fls) |>
    lapply(create_errors) |> 
    bind_rows(.id = "column_label")

errors_concat$N_lab <- factor(
    errors_concat$N,
    labels = c("$N = 25$", "$N = 50$", "$N = 100$")
)
errors_concat$M_lab <- factor(
    errors_concat$M,
    labels = c("$S = 25$", "$S = 50$", "$S = 100$")
)
errors_concat$lab <- errors_concat$lab <- interaction(errors_concat$N_lab, errors_concat$M_lab, sep = ' and ')

gg <- errors_concat |>
    filter((NPC == 5) | (NPC == 10)) |>
    filter(number <= 25) |> 
    mutate(NPC = as.factor(NPC)) |> 
    ggplot() +
    geom_boxplot(
        aes(x = number, y = value, colour = NPC, group = interaction(number, NPC))
    ) +
    facet_wrap(vars(lab)) +
    labs(
        x = "Eigenvalues",
        y = "Errors",
        colour = "Number of univariate\ncomponents estimated"
    ) +
    ylim(0, 1.1) + 
    ylab("Err$(\\widehat{\\nu}_m)$") +
    see::theme_modern() +
    theme(
        legend.position = "bottom",
        legend.key.size = unit(1.5, 'cm'),
        strip.text = element_text(size = 14)
    )

tikzDevice::tikz(
    filename = paste0(GRAPHS, '/ncomp.tex'), 
    width = 10, height = 10, 
    standAlone = TRUE, sanitize = FALSE
)
plot(gg)
dev.off()
