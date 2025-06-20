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
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111306.f003.grib2',        
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f000.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f003.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f006.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f009.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f012.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f015.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f018.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f021.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f024.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f027.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f030.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f033.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f036.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f039.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f042.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f045.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f048.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f051.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f054.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f057.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f060.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f063.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f066.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f069.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f072.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f075.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f078.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f081.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f084.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f087.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f090.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f093.grib2',
  'https://data.rda.ucar.edu/ds084.1/2018/20181113/gfs.0p25.2018111312.f096.grib2'
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
