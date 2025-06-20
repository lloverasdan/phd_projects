## Specify shell
#!/bin/bash

## Required PBS directives
#PBS -A ONRDC43015529
#PBS -q debug
#PBS -l select=32:ncpus=128:mpiprocs=128
#PBS -l walltime=00:30:00

## Optional PBS directives
#PBS -N 300
#PBS -j eo

## Set up the environment
module swap PrgEnv-cray PrgEnv-intel
module load cray-parallel-netcdf

## Change to working directory
cd ${WORKDIR}/pl_shear/WRF_pl/WRF_300/test/em_real

## Execute
aprun -n 4096 ./wrf.exe

## Clean-up
mv rsl.error.0000 rsl.error
mv rsl.out.0000 rsl.out
rm -rf rsl.error.*
rm -rf rsl.out.*
