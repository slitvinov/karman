#!/bin/sh
#SBATCH --constraint gpu
#SBATCH --ntasks 2
#SBATCH --time 1

. /etc/profile
module load daint-mc
module load cray-mpich
make MPICC=cc

srun ./cylinder -v -r 220 -l 7 -m 11 -p 10 -e 2600 -f force.da
