#!/bin/sh 

for n_curves in 25 50 100; do
    for n_points in 25 50 100; do
        for npc_univ in 1 5 10 25; do
            Rscript not-so-simple-simulation-sparse.R --n_simu 500 --n_curves $n_curves --n_features 5 --n_points $n_points --n_components 50 --npc_univ $npc_univ
        done
    done
done
