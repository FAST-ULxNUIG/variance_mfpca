#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Script for "A note on the number of components retained for multivariate
# functional principal components analysis" for the application.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Load packages
library(fda)
library(MFPCA)
library(tidyverse)

# -----------------------------------------------------------------------------
# Get the data
data(CanadianWeather)

temperature <- CanadianWeather$dailyAv[,,'Temperature.C']
precipitation <- CanadianWeather$dailyAv[,,'Precipitation.mm']


# -----------------------------------------------------------------------------
# Create the functional data objects
temperature <- funData(argvals = 1:365, X = t(temperature))
precipitation <- funData(argvals = 1:365, X = t(precipitation))

weather <- multiFunData(temperature, precipitation)


# -----------------------------------------------------------------------------
# Univariate expansions

# Temperature
ufpca_temp_medium <- PACE(temperature, nbasis = 10, pve = 0.9)
print(ufpca_temp_medium$npc)
ufpca_temp_large <- PACE(temperature, nbasis = 10, pve = 0.999)
print(ufpca_temp_large$npc)

# Precipitation
ufpca_prec_medium <- PACE(precipitation, nbasis = 10, pve = 0.9)
print(ufpca_prec_medium$npc)
ufpca_prec_large <- PACE(precipitation, nbasis = 10, pve = 0.999)
print(ufpca_prec_large$npc)


# -----------------------------------------------------------------------------
# Compute MFPCA
npc <- 4
uniExpansions <- list(
    list(type = 'uFPCA', pve = 0.9, nbasis = 10),
    list(type = 'uFPCA', pve = 0.9, nbasis = 10)
)
results_medium <- MFPCA(weather, M = npc, uniExpansions = uniExpansions)


uniExpansions <- list(
    list(type = 'uFPCA', pve = 0.999, nbasis = 10),
    list(type = 'uFPCA', pve = 0.999, nbasis = 10)
)
results_large <- MFPCA(weather, M = npc, uniExpansions = uniExpansions)

# -----------------------------------------------------------------------------
# Results
summary(results_medium)
summary(results_large)


functions_medium <- results_medium$functions
functions_large <- flipFuns(functions_medium, results_large$functions)


# For temperature
df <- data.frame(
    G = as.factor(rep(1:npc, each = 365)),
    X = rep(1:365, npc),
    Y_m = as.vector(t(functions_medium[[1]]@X)),
    Y_l = as.vector(t(functions_large[[1]]@X))
) |> reshape2::melt(c('G', 'X'))
gg <- ggplot(df, aes(x = X, y = value, group = interaction(G, variable))) +
    geom_line(aes(color = G, linetype = variable), linewidth = 2) +
    scale_color_manual(
        values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40"),
        name = "Eigenfunctions",
        labels = c("1st", "2nd", "3rd", "4th")
    ) +
    scale_linetype_manual(
        values = c("solid", "dotted"),
        name = "Univariate expansions",
        labels = c("0.9\\%", "0.999\\%")
    ) +
    labs(x = "Day of year", y = '') +
    see::theme_modern() +
    theme(
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.justification = 'center',
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.key.size = unit(3, "line"),
        strip.text = element_text(size = 20)
    ) +
    guides(color = guide_legend(order = 2), linetype = guide_legend(order = 1))

name <- './graphs/canadian_weather/temperature_eigen.tex'
tikzDevice::tikz(
    filename = name, 
    width = 10, height = 6.18047, 
    standAlone = TRUE, sanitize = FALSE
)
plot(gg)
dev.off()

# For precipitation
df <- data.frame(
    G = as.factor(rep(1:npc, each = 365)),
    X = rep(1:365, npc),
    Y_m = as.vector(t(functions_medium[[2]]@X)),
    Y_l = as.vector(t(functions_large[[2]]@X))
) |> reshape2::melt(c('G', 'X'))
gg <- ggplot(df, aes(x = X, y = value, group = interaction(G, variable))) +
    geom_line(aes(color = G, linetype = variable), linewidth = 2) +
    scale_color_manual(
        values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40"),
        name = "Eigenfunctions",
        labels = c("1st", "2nd", "3rd", "4th")
    ) +
    scale_linetype_manual(
        values = c("solid", "dotted"),
        name = "Univariate expansions",
        labels = c("0.9\\%", "0.999\\%")
    ) +
    labs(x = "Day of year", y = '') +
    see::theme_modern() +
    theme(
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.justification = 'center',
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.key.size = unit(3, "line"),
        strip.text = element_text(size = 20)
    ) +
    guides(color = guide_legend(order = 2), linetype = guide_legend(order = 1))

name <- './graphs/canadian_weather/precipitation_eigen.tex'
tikzDevice::tikz(
    filename = name, 
    width = 10, height = 6.18047, 
    standAlone = TRUE, sanitize = FALSE
)
plot(gg)
dev.off()
