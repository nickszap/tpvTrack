#module for storing the files/parameters/... for the tracking code.
#Is it necessary to reload the import to re-call all functions?

import glob
import os, errno
import numpy as np

rEarth = 6370.e3 #radius of spherical Earth (m)
dFilter = 300.e3 #radius for whether nbr extremum is regional extremum
areaOverlap = .1 #fraction of tpv area overlap for determining correspondence

latThresh = 45.*np.pi/180. #segment N of this latitude
trackMinMaxBoth = 0 #0-min, 1-max, 2-both
info = 'eraI_45N_test'

fDirData = '/data02/cases/summer2006/eraI/pv/'
filesData = sorted(glob.glob(fDirData+'eraI_theta-u-v_2pvu_2006-07-20*.nc'), key=os.path.getmtime)
deltaT = 6.*60.*60. #timestep (s)
#select time intervals within filesData[iFile]...end[-1] means use all times
iTimeStart_fData = [48]
iTimeEnd_fData = [88]
if (True):
  nFiles = len(filesData)
  if (len(iTimeStart_fData) != nFiles or len(iTimeEnd_fData) != nFiles):
    print "Uhoh, wrong iTime*_data settings in my_settings.py"

fDirSave = '/data02/cases/test_segment/testUnified/200608/'
if not os.path.exists(fDirSave):
    os.makedirs(fDirSave)

fMesh = filesData[0]  
fMetr = fDirSave+'fields_debug.nc'
fSeg = fDirSave+'seg_debug.nc'
fCorr = fDirSave+'correspond_debug.pkl'
fTrack = fDirSave+'tracks_debug.txt'
fMetrics = fDirSave+'metrics_debug.nc'

doPreProc = True
doSeg = False
doMetrics = False
doCorr = False
doTracks = False

def silentremove(filename):
  #from http://stackoverflow.com/questions/10840533/most-pythonic-way-to-delete-a-file-which-may-not-exist
  print "Removing file (if it exists): ",filename
  try:
      os.remove(filename)
  except OSError as e: # this would be "except OSError, e:" before Python 2.6
      if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
          raise # re-raise exception if a different error occured
