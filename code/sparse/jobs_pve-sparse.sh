#!/bin/sh 

for n_curves in 25 50 100; do
    for n_points in 25 50 100; do
        for npc_univ in 0.5 0.7 0.9 0.95 0.99; do
            Rscript for_pve-sparse.R --n_simu 500 --n_curves $n_curves --n_features 5 --n_points $n_points --n_components 50 --npc_univ $npc_univ
        done
    done
done
