## Specify shell
#!/bin/bash

## Required PBS directives
#PBS -A ONRDC43015529
#PBS -q transfer
#PBS -l select=1:ncpus=1
#PBS -l walltime=48:00:00

## Optional PBS directives
#PBS -N arch
#PBS -j eo

## Set up the environment
#PBS -V

## Execution

cd ${WORKDIR}/smallstuff
archive get -C smallstuff wrfin_ctl
archive get -C smallstuff wrfout_ctl

cd ${WORKDIR}/smallstuff/long_runs/ctl_48h/
archive get -C smallstuff/long_runs/ctl_48h wrfout_d01_2021-01-03_00_00_00

cd ${WORKDIR}/smallstuff/long_runs/adj_48h/
archive get -C smallstuff/long_runs/adj_48h wrfout_d01_2021-01-03_00_00_00

cd ${WORKDIR}/smallstuff/long_runs/adj_48h_10/
archive get -C smallstuff/long_runs/adj_48h_10 wrfout_d01_2021-01-03_00_00_00

cd ${WORKDIR}/smallstuff/long_runs/adj_48h_100/
archive get -C smallstuff/long_runs/adj_48h_100 wrfout_d01_2021-01-03_00_00_00

cd ${WORKDIR}/smallstuff/long_runs/wave_48h/
archive get -C smallstuff/long_runs/wave_48h wrfout_d01_2021-01-03_00_00_00

cd ${WORKDIR}/smallstuff/long_runs/wave_48h_10/
archive get -C smallstuff/long_runs/wave_48h_10 wrfout_d01_2021-01-03_00_00_00

cd ${WORKDIR}/smallstuff/long_runs/wave_48h_100/
archive get -C smallstuff/long_runs/wave_48h_100 wrfout_d01_2021-01-03_00_00_00
