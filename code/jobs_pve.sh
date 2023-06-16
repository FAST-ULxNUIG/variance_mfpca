#!/bin/sh 
#SBATCH --time=12:00:00
#SBATCH --nodes=2
#SBATCH -A ngmat008c
#SBATCH -p ProdQ
#SBATCH -o output.log
#SBATCH --mail-user=steven.golovkine@ul.ie
#SBATCH --mail-type=BEGIN,END

module load r/4.1.2

for n_curves in 25 50 100; do
    for n_points in 25 50 100; do
        for npc_univ in 0.5 0.7 0.9 0.95 0.99; do
            Rscript for_pve.R --n_simu 4 --n_curves $n_curves --n_features 5 --n_points $n_points --n_components 50 --npc_univ $npc_univ
        done
    done
done
