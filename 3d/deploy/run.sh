#!/bin/sh
#SBATCH --constraint gpu
#SBATCH --account s1160

. /etc/profile
module load daint-mc
module load cray-mpich

make MPICC=cc 'MPICCFLAGS = -O2 -g' && {
    r=2000
    mkdir -p $r
    srun ./cylinder -v -r $r -l 8 -m 11 -p 100 -e 2600 -f $r/force.dat -o $r/h
}
