#module for storing the files/parameters/... for the tracking code.
#Is it necessary to reload the import to re-call all functions?

import glob
import os
import numpy as np

rEarth = 6370.e3 #radius of spherical Earth (m)
dFilter = 300.e3 #radius for whether nbr extremum is regional extremum
areaOverlap = .2 #fraction of tpv area overlap for determining correspondence

latThresh = 45.*np.pi/180. #segment N of this latitude
trackMinMaxBoth = 0 #0-min, 1-max, 2-both
info = 'eraI_45N_test'

fDirData = '/data02/cases/summer2006/eraI/pv/'
filesData = sorted(glob.glob(fDirData+'eraI_theta-u-v_2pvu_2006-07-20*.nc'), key=os.path.getmtime)
deltaT = 6.*60.*60. #timestep (s)

#iTimeStart = 0
#iTimeEnd = 19

fDirSave = '/data02/cases/test_segment/testUnified/'
fMesh = filesData[0]  
fMetr = fDirSave+'fields_debug.nc'
fSeg = fDirSave+'seg_debug.nc'
fCorr = fDirSave+'correspond_debug.txt'
fTrack = fDirSave+'tracks_debug.txt'
fMetrics = fDirSave+'metrics_debug.nc'
