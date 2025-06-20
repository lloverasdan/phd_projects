## Specify shell
#!/bin/bash

## Required PBS directives
#PBS -A NRLMR03795024
#PBS -q transfer
#PBS -l select=1:ncpus=1
#PBS -l walltime=48:00:00

## Optional PBS directives
#PBS -N tar_nov2018
#PBS -j eo

## Set up the environment
#PBS -V

## Execution
cd ${WORKDIR}/nov2018/

# tar -czvf coamps_files_may21.tar.gz coamps_files
# archive put -C nov2018 coamps_files_may21.tar.gz

# tar -czvf gfs_files_may21.tar.gz gfs_files
# archive put -C nov2018 gfs_files_may21.tar.gz

# tar -czvf gefs_files_may21.tar.gz gefs_data
# archive put -C nov2018 gefs_files_may21.tar.gz

# tar -czvf gfs_nav_may21.tar.gz gfs_nav
# archive put -C nov2018 gfs_nav_may21.tar.gz

# archive get -C nov2018 proc_oct.tar.gz
# tar -xvf proc_oct.tar.gz

# archive get -C nov2018 gfs_files_may21.tar.gz
# tar -xvf gfs_files_may21.tar.gz

archive get -C nov2018 coamps_files_may21.tar.gz
tar -xvf coamps_files_may21.tar.gz

# archive get -C nov2018 gefs_files_may21.tar.gz
# tar -xvf gefs_files_may21.tar.gz

# archive get -C nov2018 gfs_nav_may21.tar.gz
# tar -xvf gfs_nav_may21.tar.gz