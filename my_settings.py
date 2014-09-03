#module for storing the files/parameters/... for the tracking code.
#Is it necessary to reload the import to re-call all functions?

import glob
import os

rEarth = 6370.e3 #radius of spherical Earth (m)
dFilter = 300.e3 #radius for whether extremum is extremum

fDirSave = '/data02/cases/test_segment/stevenCase/'  

fDirData = '/data02/cases/2006/eraI/pv/'
filesData = sorted(glob.glob(fDirData+'eraI_theta-u-v_2pvu_2006-07-20*.nc'), key=os.path.getmtime)

filesMetr = sorted(glob.glob(fDir+'/seg/fields_*.npz'), key=os.path.getmtime)
filesSeg = sorted(glob.glob(fDir+'/seg/seg_*.npz'), key=os.path.getmtime)
filesMetrics = sorted(glob.glob(fDir+'/seg/metrics_*.npz'), key=os.path.getmtime)
filesTrack = sorted(glob.glob(fDir+'/track/track_seg*.npz'), key=os.path.getmtime)
