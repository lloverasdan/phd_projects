## Specify shell
#!/bin/bash

## Required PBS directives
#PBS -A ONRDC43015529
#PBS -q debug
#PBS -l select=1:ncpus=1
#PBS -l walltime=00:30:00

## Optional PBS directives
#PBS -N gefs_anl
#PBS -j eo

## Set up the environment
#PBS -V

## Execute
cd ${WORKDIR}/nov2018/gefs_data
${HOME}/aws-cli/bin/aws s3 --no-sign-request cp s3://noaa-gefs-pds/gefs.20181115/12/pgrb2a/ . --recursive --exclude "*" --include "*aanl*"
${HOME}/aws-cli/bin/aws s3 --no-sign-request cp s3://noaa-gefs-pds/gefs.20181116/12/pgrb2a/ . --recursive --exclude "*" --include "*aanl*"