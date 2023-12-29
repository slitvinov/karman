#!/bin/sh
module load gcc openmpi
set -x
r=$1
d=$SCRATCH/koumoutsakos_lab/slitvinov/$r
make MPICC=mpicc 'MPICCFLAGS = -O2 -g' && {
    mkdir -p $d
    exec srun -n $SLURM_NTASKS --mpi=pmix ./cylinder -v -r $r -l 7 -m 9 -p 10 -e 10 -f $d/force.dat -S cylinder -o $d/h -z 2.5
}
