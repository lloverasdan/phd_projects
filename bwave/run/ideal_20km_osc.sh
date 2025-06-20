## Specify shell
#!/bin/bash

## Required PBS directives
#PBS -A ONRDC43015529
#PBS -q debug
#PBS -l select=4:ncpus=128:mpiprocs=128
#PBS -l walltime=00:15:00

## Optional PBS directives
#PBS -N osc_wrf
#PBS -j eo

## Set up the environment
module swap PrgEnv-cray PrgEnv-intel
module load cray-netcdf-hdf5parallel
module load cray-parallel-netcdf

## Change to working directory
cd ${WORKDIR}/bwave_nov/WRF_surface/WRFV3_dry/run

## Execute
aprun -n 1024 ./wrf.exe

## Clean-up
mv rsl.error.0000 rsl.error
mv rsl.out.0000 rsl.out
rm -rf rsl.error.*
rm -rf rsl.out.*
