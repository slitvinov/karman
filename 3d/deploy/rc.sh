#!/bin/sh
module load gcc openmpi
set -x
r=$1
d=$SCRATCH/koumoutsakos_lab/slitvinov/$r
make MPICC=mpicc 'MPICCFLAGS = -O2 -g' && {
    mkdir -p $d
    exec srun -n $SLURM_NTASKS --mpi=pmix ./cylinder -v -r $r -l 9 -m 12 -p 10 -e 2600 -f $d/force.dat -o $d/h
}
