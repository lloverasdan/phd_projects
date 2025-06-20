## Specify shell
#!/bin/bash

## Required PBS directives
#PBS -A ONRDC43015529
#PBS -q transfer
#PBS -l select=1:ncpus=1
#PBS -l walltime=24:00:00

## Optional PBS directives
#PBS -N arch
#PBS -j eo

## Set up the environment
#PBS -V

## Execution

cd ${WORKDIR}/bwave
tar -czvf 20km_files_mar.tar.gz 20km_files
archive put -C bwave 20km_files_mar.tar.gz

tar -czvf proc_mar.tar.gz processed
archive put -C bwave proc_mar.tar.gz

# cd ${WORKDIR}/bwave
# archive get -C bwave proc_nov.tar.gz
# tar -xzvf proc_nov.tar.gz

# cd ${WORKDIR}/bwave
# archive get -C bwave 20km_files.tar.gz
# tar -xzvf 20km_files.tar.gz

# cd ${WORKDIR}/bwave/4km_files/input
# archive get -C bwave/4km_files/input wrfin_standard
# archive get -C bwave/4km_files/input wrfin_surface
# archive get -C bwave/4km_files/input wrfin_barotropic

# cd ${WORKDIR}/bwave/4km_files/output
# archive get -C bwave/4km_files/output wrfout_standard
# archive get -C bwave/4km_files/output wrfout_surface
# archive get -C bwave/4km_files/output wrfout_barotropic

# cd ${WORKDIR}/bwave/4km_files/restart_standard
# archive get -C bwave/4km_files/restart_standard wrfrst_d01_2021-01-02_00_00_00
# archive get -C bwave/4km_files/restart_standard wrfrst_d01_2021-01-03_00_00_00
# archive get -C bwave/4km_files/restart_standard wrfrst_d01_2021-01-04_00_00_00
# archive get -C bwave/4km_files/restart_standard wrfrst_d01_2021-01-05_00_00_00
# archive get -C bwave/4km_files/restart_standard wrfrst_d01_2021-01-06_00_00_00
# archive get -C bwave/4km_files/restart_standard wrfrst_d01_2021-01-07_00_00_00
# archive get -C bwave/4km_files/restart_standard wrfrst_d01_2021-01-08_00_00_00
# archive get -C bwave/4km_files/restart_standard wrfrst_d01_2021-01-09_00_00_00

# cd ${WORKDIR}/bwave/4km_files/restart_surface
# archive get -C bwave/4km_files/restart_surface wrfrst_d01_2021-01-02_00_00_00
# archive get -C bwave/4km_files/restart_surface wrfrst_d01_2021-01-03_00_00_00
# archive get -C bwave/4km_files/restart_surface wrfrst_d01_2021-01-04_00_00_00
# archive get -C bwave/4km_files/restart_surface wrfrst_d01_2021-01-05_00_00_00
# archive get -C bwave/4km_files/restart_surface wrfrst_d01_2021-01-06_00_00_00
# archive get -C bwave/4km_files/restart_surface wrfrst_d01_2021-01-07_00_00_00
# archive get -C bwave/4km_files/restart_surface wrfrst_d01_2021-01-08_00_00_00
# archive get -C bwave/4km_files/restart_surface wrfrst_d01_2021-01-09_00_00_00

# cd ${WORKDIR}/bwave/4km_files/restart_barotropic
# archive get -C bwave/4km_files/restart_barotropic wrfrst_d01_2021-01-02_00_00_00
# archive get -C bwave/4km_files/restart_barotropic wrfrst_d01_2021-01-03_00_00_00
# archive get -C bwave/4km_files/restart_barotropic wrfrst_d01_2021-01-04_00_00_00
# archive get -C bwave/4km_files/restart_barotropic wrfrst_d01_2021-01-05_00_00_00
# archive get -C bwave/4km_files/restart_barotropic wrfrst_d01_2021-01-06_00_00_00
# archive get -C bwave/4km_files/restart_barotropic wrfrst_d01_2021-01-07_00_00_00
# archive get -C bwave/4km_files/restart_barotropic wrfrst_d01_2021-01-08_00_00_00
# archive get -C bwave/4km_files/restart_barotropic wrfrst_d01_2021-01-09_00_00_00
