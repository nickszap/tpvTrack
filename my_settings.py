#module for storing the files/parameters/... for the tracking code.
#Is it necessary to reload the import to re-call all functions?

import glob
import os
import numpy as np

rEarth = 6370.e3 #radius of spherical Earth (m)
dFilter = 300.e3 #radius for whether extremum is extremum

latThresh = 45.*np.pi/180. #segment N of this latitude
trackMinMaxBoth = 0 #0-min, 1-max, 2-both
mesh = 'll' #'ll','wrf','mpas'

fDirData = '/data02/cases/summer2006/eraI/pv/'
filesData = sorted(glob.glob(fDirData+'eraI_theta-u-v_2pvu_2006-07-20*.nc'), key=os.path.getmtime)
deltaT = 6.*60.*60. #time frequency (s)

fDirSave = '/data02/cases/test_segment/testUnified/'
info = 'eraI_45N_test'
