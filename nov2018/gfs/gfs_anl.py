#!/usr/bin/env python
""" 
Python script to download selected files from rda.ucar.edu.
After you save the file, don't forget to make it executable
i.e. - "chmod 755 <name_of_script>"
"""
import sys, os
from urllib.request import build_opener

opener = build_opener()

filelist = [
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111300.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111306.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f000.grib2',    
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111318.f000.grib2',
    
  'https://data.rda.ucar.edu/ds084.1/2018/20181114/gfs.0p25.2018111400.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181114/gfs.0p25.2018111406.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181114/gfs.0p25.2018111412.f000.grib2',    
  'https://data.rda.ucar.edu/ds084.1/2018/20181114/gfs.0p25.2018111418.f000.grib2',    
    
  'https://data.rda.ucar.edu/ds084.1/2018/20181115/gfs.0p25.2018111500.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181115/gfs.0p25.2018111506.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181115/gfs.0p25.2018111512.f000.grib2',    
  'https://data.rda.ucar.edu/ds084.1/2018/20181115/gfs.0p25.2018111518.f000.grib2',
    
  'https://data.rda.ucar.edu/ds084.1/2018/20181116/gfs.0p25.2018111600.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181116/gfs.0p25.2018111606.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181116/gfs.0p25.2018111612.f000.grib2',    
  'https://data.rda.ucar.edu/ds084.1/2018/20181116/gfs.0p25.2018111618.f000.grib2',     

  'https://data.rda.ucar.edu/ds084.1/2018/20181117/gfs.0p25.2018111700.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181117/gfs.0p25.2018111706.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181117/gfs.0p25.2018111712.f000.grib2',    
  'https://data.rda.ucar.edu/ds084.1/2018/20181117/gfs.0p25.2018111718.f000.grib2', 
    
]

for file in filelist:
    ofile = os.path.basename(file)
    sys.stdout.write("downloading " + ofile + " ... ")
    sys.stdout.flush()
    infile = opener.open(file)
    outfile = open(ofile, "wb")
    outfile.write(infile.read())
    outfile.close()
    sys.stdout.write("done\n")
