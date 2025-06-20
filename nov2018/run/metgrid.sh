## Specify shell
#!/bin/bash

## Required PBS directives
#PBS -A ONRDC43015529
#PBS -q debug
#PBS -l select=2:ncpus=128:mpiprocs=128
#PBS -l walltime=00:30:00

## Optional PBS directives
#PBS -N 4km_metgrid
#PBS -j eo

## Set up the environment
module swap PrgEnv-cray PrgEnv-intel
module load cray-parallel-netcdf

## Change to working directory
cd ${WORKDIR}/pl_shear/WRF_pl/WPS

## Execute
aprun -n 256 ./metgrid.exe
rm -rf metgrid.log.0*
