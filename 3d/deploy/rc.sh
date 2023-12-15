#!/bin/sh
. /etc/profile
module load gcc mpich

set -x
make MPICC=cc 'MPICCFLAGS = -O2 -g' && {
    r=$1
    mkdir -p $r
    srun ./cylinder -v -r $r -l 8 -m 12 -p 100 -e 2600 -f $r/force.dat -o $r/h
}
