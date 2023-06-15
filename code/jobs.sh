#!/bin/sh 
#SBATCH --time=12:00:00
#SBATCH --nodes=2
#SBATCH -A ngmat008c
#SBATCH -p DevQ

for n_curves in 25 50 100; do
    for n_points in 25 50 100; do
        for npc_univ in 1 5 10 25; do
            Rscript not-so-simple-simulation.R --n_simu 4 --n_curves $n_curves --n_features 5 --n_points $n_points --n_components 50 --npc_univ $npc_univ
        done
    done
done
