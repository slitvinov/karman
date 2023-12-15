#!/bin/sh
. /etc/profile
module load intel intelmpi

set -x
r=$1
d=$SCRATCH/koumoutsakos_lab/slitvinov/$r
make MPICC=mpiicx 'MPICCFLAGS = -O2 -g -Wno-deprecated-non-prototype' && {
    mkdir -p $d
    exec srun -n $SLURM_NTASKS --mpi=pmi2 ./cylinder -v -r $r -l 8 -m 12 -p 100 -e 2600 -f $d/force.dat -o $d/h
}
