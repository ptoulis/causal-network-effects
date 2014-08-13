#!/bin/bash
#SBATCH -J parallel-testing # name for job array
#SBATCH -o odyssey/testing_%a.out #Standard output
#SBATCH -e odyssey/testing_%a.err #Standard error
#SBATCH -p general #Partition
#SBATCH -t 00:20:00 #Running time of 20 mins.
#SBATCH --mem-per-cpu 3000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine

Rscript main-odyssey.R 200
