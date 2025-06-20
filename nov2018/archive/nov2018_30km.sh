## Specify shell
#!/bin/bash

## Required PBS directives
#PBS -A NRLMR03795024
#PBS -q transfer
#PBS -l select=1:ncpus=1
#PBS -l walltime=48:00:00

## Optional PBS directives
#PBS -N 30km_nov2018
#PBS -j eo

## Set up the environment
#PBS -V

## Execution
# cd ${WORKDIR}/nov2018/30km_files
# archive get -C nov2018/30km_files wrfbdy_d01

cd ${WORKDIR}/nov2018/30km_files/ctl
archive get -C nov2018/30km_files/ctl wrfout_d01_2018-11-13_12_00_00

# cd ${WORKDIR}/nov2018/30km_files/gfs
# archive get -C nov2018/30km_files/gfs wrfin_d01_2018-11-13_12_00_00
# archive get -C nov2018/30km_files/gfs wrfin_d01_2018-11-14_12_00_00
# archive get -C nov2018/30km_files/gfs wrfin_d01_2018-11-15_12_00_00
# archive get -C nov2018/30km_files/gfs wrfin_d01_2018-11-16_12_00_00
# archive get -C nov2018/30km_files/gfs wrfin_d01_2018-11-17_12_00_00

cd ${WORKDIR}/nov2018/30km_files/adj_full
archive get -C nov2018/30km_files/adj_full wrfout_d01_2018-11-13_12_00_00

# cd ${WORKDIR}/nov2018/30km_files/adj_large
# archive get -C nov2018/30km_files/adj_large wrfout_d01_2018-11-13_12_00_00

# cd ${WORKDIR}/nov2018/30km_files/adj_small
# archive get -C nov2018/30km_files/adj_small wrfout_d01_2018-11-13_12_00_00

# cd ${WORKDIR}/nov2018/30km_files/adj_box
# archive get -C nov2018/30km_files/adj_box wrfout_d01_2018-11-13_12_00_00

# cd ${WORKDIR}/nov2018/30km_files/adj_hole
# archive get -C nov2018/30km_files/adj_hole wrfout_d01_2018-11-13_12_00_00

# cd ${WORKDIR}/nov2018/30km_files/adj_cape_box
# archive get -C nov2018/30km_files/adj_cape_box wrfout_d01_2018-11-13_12_00_00

# cd ${WORKDIR}/nov2018/30km_files/adj_cape_hole
# archive get -C nov2018/30km_files/adj_cape_hole wrfout_d01_2018-11-13_12_00_00
