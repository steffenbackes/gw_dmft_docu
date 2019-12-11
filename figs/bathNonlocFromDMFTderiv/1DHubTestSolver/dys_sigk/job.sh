#!/bin/bash
#MSUB -r 1d_dyson                # jobname
#MSUB -n 96                      # No. tasks
###MSUB -N 1                       # No. nodes 
###MSUB -c 1                           # Cores per task
###MSUB -M 4000                    # Memory per core in Mbyte
####MSUB -Q long                     # this allows 3 days max, otherwise 24h
#MSUB -T 1000                # Max. runtime in sec.
#MSUB -o job.out
#MSUB -e job.err
#MSUB -q skylake
####MSUB -q hybrid
####MSUB -q xlarge
#MSUB -A gen1393
#MSUB -m scratch

set -x
cd ${BRIDGE_MSUB_PWD}
export OMP_NUM_THREADS=1
#
#


module load hdf5/1.8.20
module load cmake/3.13.3
module load fftw3/mkl/17.0.6.256
module load boost/1.63.0

ccc_mprun alps_cthyb alps_parm.dat
